// main.cpp
#include <mpi.h>
#include <vector>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>   // memcpy
#include <algorithm> // min
#include <functional>
#include <omp.h>
#include <chrono>
#include <iostream> 
#include <string>
#include <unordered_map>
#include <mutex>
#include "page_ops.hpp"
#include "message.hpp"

#include "utilities.hpp"  // owner_of_shape(...)
#include "othello6.hpp"   // swap to othello4.hpp when running 4x4

// Throughput knobs (correctness does not depend on specific values)
static constexpr int BATCH_MAX  = 4096*2; // send/recv chunk size
static constexpr int SEND_DEPTH = 8;    // in-flight slots per destination
static constexpr int TAG_HITS   = 1;
static constexpr int TAG_MODULE_TO_MAIN = 2;
static constexpr int TAG_MAIN_TO_MODULE = 3;
static constexpr int TAG_SEND_TOKEN = 4;

// How often to (re)kick a non-blocking global check if nothing else triggers it.
static constexpr int REDUCE_INTERVAL_ITERS = 128;
// If we've had this many consecutive idle loops, kick a non-blocking reduce even sooner.
static constexpr int IDLE_TRIGGER_ITERS    = 8;

struct Outgoing {
    // In-flight send slots (never modify a buf while its req != MPI_REQUEST_NULL)
    std::array<std::vector<HitMsg>, SEND_DEPTH> bufs;
    std::array<MPI_Request, SEND_DEPTH> reqs;
    // Small staging per destination; non-blocking pumps keep it small
    std::vector<HitMsg> staging;
    std::size_t head = 0;
};

// --- timing utility ---
struct TimerAggregator {
    std::unordered_map<std::string, double> totals;
    std::unordered_map<std::string, std::size_t> counts;
    std::mutex mtx;

    static TimerAggregator& instance() {
        static TimerAggregator agg;
        return agg;
    }

    void add(const std::string& name, double sec) {
        std::lock_guard<std::mutex> lock(mtx);
        totals[name] += sec;
        counts[name] += 1;
    }

    void report() {
        std::cout << "\n=== Timing summary ===\n";
        for (auto& kv : totals) {
            const std::string& name = kv.first;
            double total = kv.second;
            std::size_t n = counts[name];
            std::cout << name << ": "
                    << total << " s total, "
                    << (total / n) << " s avg over "
                    << n << " calls\n";
        }
    }
};

// wrapper for timing blocks
template <typename Func>
auto time_block(const std::string& name, Func&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    auto result = f();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    TimerAggregator::instance().add(name, diff.count());
    return result;
}

// ---- MPI datatype for HitMsg (portable) ----
static inline MPI_Datatype make_HitMsg_type() {
    MPI_Datatype T;
    const int N = 3;
    int          bl[N]  = {1,1,1};
    MPI_Aint     dis[N] = {};
    MPI_Datatype ty[N]  = {MPI_UINT64_T, MPI_UINT32_T, MPI_UINT32_T};

    HitMsg probe{};
    MPI_Aint base=0, a0=0, a1=0, a2=0;
    MPI_Get_address(&probe,        &base);
    MPI_Get_address(&probe.shape,  &a0);
    MPI_Get_address(&probe.page,   &a1);
    MPI_Get_address(&probe.bit,    &a2);
    dis[0] = a0 - base; dis[1] = a1 - base; dis[2] = a2 - base;

    MPI_Type_create_struct(N, bl, dis, ty, &T);
    MPI_Type_commit(&T);
    return T;
}

// ---- Canonicalize -> HitMsg ----
static inline HitMsg make_HitMsg_from_pos(const game::Encoding& position) noexcept {
    HitMsg m{};
    const auto c  = game::canonical(position);
    const auto sh = game::shape(c);
    const auto h  = game::hash(c);
    m.shape = sh;
    m.page  = static_cast<std::uint32_t>(h >> pageops::PAGE_BITS);
    m.bit   = static_cast<std::uint32_t>(h & ((1u << pageops::PAGE_BITS) - 1u));
    return m;
}

static inline int reap_one_free_slot_nb(Outgoing& o) {
    for (int k = 0; k < SEND_DEPTH; ++k)
        if (o.reqs[k] == MPI_REQUEST_NULL && o.bufs[k].empty()) return k;
    int index = MPI_UNDEFINED, flag = 0;
    MPI_Testany(SEND_DEPTH, o.reqs.data(), &index, &flag, MPI_STATUS_IGNORE);
    if (flag && index != MPI_UNDEFINED) {
        o.reqs[index] = MPI_REQUEST_NULL;
        o.bufs[index].clear();
        return index;
    }
    return -1;
}

static inline void try_send_dest_nb(Outgoing& o, int dest,
                                    MPI_Datatype HIT_T, int tag,
                                    unsigned long long& sent_hits)
{
    // nothing staged or fully drained
    if (o.head >= o.staging.size()) {
        o.staging.clear();
        o.head = 0;
        return;
    }

    int slot = reap_one_free_slot_nb(o);
    if (slot < 0) return;

    const std::size_t avail = o.staging.size() - o.head;
    const int n = (int)std::min<std::size_t>(avail, (std::size_t)BATCH_MAX);

    o.bufs[slot].resize(n);
    std::memcpy(o.bufs[slot].data(), o.staging.data() + o.head, n * sizeof(HitMsg));
    o.head += (std::size_t)n;

    MPI_Isend(o.bufs[slot].data(), n, HIT_T, dest, tag, MPI_COMM_WORLD, &o.reqs[slot]);
    sent_hits += (unsigned long long)n;
}

// ---- Pump all destinations once (non-blocking) ----
static inline void pump_all_nb(std::vector<Outgoing>& outs, MPI_Datatype HIT_T, int tag,
                               unsigned long long& sent_hits)
{
    const int W = (int)outs.size();
    for (int d = 0; d < W; ++d) try_send_dest_nb(outs[d], d, HIT_T, tag, sent_hits);
}

// ---- Reap completed sends (free slots) ----
static inline int progress_sends(std::vector<Outgoing>& outs) {
    int progressed = 0;
    const int W = (int)outs.size();
    for (int d = 0; d < W; ++d) {
        auto& o = outs[d];
        for (int k = 0; k < SEND_DEPTH; ++k) {
            if (o.reqs[k] != MPI_REQUEST_NULL) {
                int done = 0;
                MPI_Test(&o.reqs[k], &done, MPI_STATUS_IGNORE);
                if (done) {
                    o.reqs[k] = MPI_REQUEST_NULL;
                    o.bufs[k].clear();
                    progressed = 1;
                }
            }
        }
    }
    return progressed;
}

// ---- Any local work pending (remote side only; self work is processed immediately) ----
static inline int have_work(const std::vector<Outgoing>& outs) {
    for (const auto& o : outs) {
        if (!o.staging.empty()) return 1;
        for (int k = 0; k < SEND_DEPTH; ++k) {
            if (o.reqs[k] != MPI_REQUEST_NULL) return 1;
            if (!o.bufs[k].empty())           return 1;
        }
    }
    return 0;
}

// ---- Drain completed receives ----
static inline int progress_receives(std::vector<HitMsg>& rxbuf, MPI_Request& rxreq,
                                    MPI_Datatype HIT_T, int tag,
                                    const std::function<void(const HitMsg*, int)>& handle_batch,
                                    unsigned long long& recv_hits)
{
    int progressed = 0;
    for (;;) {
        int done = 0; MPI_Status st;
        MPI_Test(&rxreq, &done, &st);
        if (!done) break;

        int count = 0; MPI_Get_count(&st, HIT_T, &count);
        if (count > 0) {
            handle_batch(rxbuf.data(), count);
            recv_hits += (unsigned long long)count;
            progressed = 1;
        }
        MPI_Irecv(rxbuf.data(), (int)rxbuf.size(), HIT_T, MPI_ANY_SOURCE, tag,
                  MPI_COMM_WORLD, &rxreq);
    }
    return progressed;
}

static inline int progress_receives_tier(std::vector<TierMsg>& rxbuf, MPI_Request& rxreq,
                                    MPI_Datatype TIER_T,
                                    const std::function<void(const TierMsg*, int, int)>& handle_batch)
{
    int progressed = 0;
    for (;;) {
        int done = 0; MPI_Status st;
        MPI_Test(&rxreq, &done, &st);
        if (!done) break;

        int count = 0; MPI_Get_count(&st, TIER_T, &count);
        if (count > 0) {
            handle_batch(rxbuf.data(), count, st.MPI_SOURCE);
            progressed = 1;
        }
        MPI_Irecv(rxbuf.data(), (int)rxbuf.size(), TIER_T, MPI_ANY_SOURCE, TAG_MODULE_TO_MAIN,
                  MPI_COMM_WORLD, &rxreq);
    }
    return progressed;
}

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int W = 0, R = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &W);
    MPI_Comm_rank(MPI_COMM_WORLD, &R);

    MPI_Datatype HIT_T = make_HitMsg_type();

    // Outgoing per-destination state
    std::vector<Outgoing> out(W);
    for (int d = 0; d < W; ++d) {
        for (int k = 0; k < SEND_DEPTH; ++k) {
            out[d].bufs[k].reserve(BATCH_MAX);
            out[d].reqs[k] = MPI_REQUEST_NULL;
        }
        out[d].staging.reserve(BATCH_MAX);
    }

    // Per-tier page contexts (lazy init)
    pageops::Context ctx_by_tier[game::MAX_TIER + 1] = {};
    auto get_ctx = [&](std::uint8_t tier) -> pageops::Context& {
        if (tier > game::MAX_TIER) {
            std::fprintf(stderr, "Error: tier %d exceeds MAX_TIER %d\n", tier, game::MAX_TIER);
            std::exit(1);
        }
        auto& ctx = ctx_by_tier[tier];
        if (!ctx.inited) {
            std::printf("Rank %d: Initializing context for tier %d\n", R, tier);
            pageops::pages_init(ctx, R, tier, "data");
        }
        return ctx;
    };    

    // One wildcard Irecv buffer
    std::vector<HitMsg> rxbuf(BATCH_MAX);
    MPI_Request rxreq = MPI_REQUEST_NULL;
    MPI_Irecv(rxbuf.data(), (int)rxbuf.size(), HIT_T, MPI_ANY_SOURCE, TAG_HITS,
              MPI_COMM_WORLD, &rxreq);

    // Global counters (termination)
    unsigned long long sent_hits = 0ULL;
    unsigned long long recv_hits = 0ULL;

    // Forward declare batch handler; self-processing is iterative (no recursion)
    std::function<void(const HitMsg*, int)> handle_batch;
    std::function<void(const TierMsg*, int, int)> handle_batch_tier;
    unsigned long long how_many_pos = 0;

    // Process a single self-dest hit and expand iteratively (DFS-like).
    // (No MPI here; only stage for remote destinations.)
    auto process_self_iter = [&](const HitMsg& seed) {
        std::vector<HitMsg> stack;
        stack.reserve(1024);
        stack.push_back(seed);

        while (!stack.empty()) {
            HitMsg msg = stack.back();
            stack.pop_back();

            const std::uint8_t tier = game::tier_of(msg.shape);
            pageops::Context& ctx   = get_ctx(tier);

            if (!pageops::set_bit(ctx, msg.shape, msg.page, msg.bit))
                continue; // already seen

            ++how_many_pos;
            if ((how_many_pos % 100000) == 0) {
                // printf("Rank %d has found %llu\n", R, how_many_pos);
                fflush(stdout);
                TimerAggregator::instance().report();
            }

            const std::uint64_t h =
                ((std::uint64_t)msg.page << pageops::PAGE_BITS) |
                 (std::uint64_t)msg.bit;
            game::Encoding pos = time_block("unhash", [&] {
                return game::unhash(msg.shape, h);
            });

            std::uint64_t moves = time_block("legal_moves", [&] {
                return game::legal_moves(pos);
            });
            if (moves) {
                while (moves) {
                    std::uint64_t mv = moves & -moves;
                    moves ^= mv;
                    const game::Encoding nxt = time_block("do_move+flip", [&] {
                        return game::flip(game::do_move(pos, mv));
                    });
                    const HitMsg m2 = make_HitMsg_from_pos(nxt);
                    const int dest = owner_of_shape(m2.shape, W);
                    if (dest == R) {
                        stack.emplace_back(m2); // keep processing locally
                    } else {
                        out[dest].staging.emplace_back(m2); // stage only
                    }
                }
            } else {
                // pass move
                const HitMsg m2 = make_HitMsg_from_pos(game::flip(pos));
                const int dest = owner_of_shape(m2.shape, W);
                if (dest == R) {
                    stack.emplace_back(m2);
                } else {
                    out[dest].staging.emplace_back(m2); // stage only
                }
            }
        }
    };

    // Batch handler: route each message
    handle_batch = [&](const HitMsg* msgs, int n) {
        for (int i = 0; i < n; ++i) {
            const HitMsg& m = msgs[i];
            const int dest = owner_of_shape(m.shape, W);
            if (dest == R) {
                process_self_iter(m);
            } else {
                out[dest].staging.emplace_back(m); // stage onlyx
            }
        }
    };

    // Near-terminal position: 
    game::Encoding test = {
        0b0000000000000000000000000000111111100001100001100001100001111110,  // Player (P)
        0b0000000000000000000000000000000000011110011110011110011110000000  // player (O)
    };

    // game::Encoding test = {
    //     0b0000000000000000000000000000010100100001000001100001000001010100,  // Player (P)
    //     0b0000000000000000000000000000000000010010001010011010010100000000  // player (O)
    // };



    // Seed from rank 0
    if (R == 0) {
        const HitMsg m0 = make_HitMsg_from_pos(test); // instead of starting position
        printf("actual result: %llu, %u, %u \n", m0.shape, m0.page, m0.bit);
        const int dest0 = owner_of_shape(m0.shape, W);
        if (dest0 == R) {
            process_self_iter(m0);
        } else {
            out[dest0].staging.emplace_back(m0); // stage only
        }
    }

    // ---- Non-blocking global check plumbing ----
    // We reduce three values at once via MPI_SUM: [Σsent, Σrecv, Σlocal_active]
    // local vector and global result
    unsigned long long local_vec[3]  = {0ULL, 0ULL, 0ULL};
    unsigned long long global_vec[3] = {0ULL, 0ULL, 0ULL};
    MPI_Request        red_req       = MPI_REQUEST_NULL;
    bool               reduce_in_flight = false;

    auto kick_reduce_nb = [&] (int local_active) {
        if (reduce_in_flight) return; // don't start a new one until previous completes
        local_vec[0] = sent_hits;
        local_vec[1] = recv_hits;
        local_vec[2] = (unsigned long long)local_active; // 0 or 1
        MPI_Iallreduce(local_vec, global_vec, 3, MPI_UNSIGNED_LONG_LONG,
                       MPI_SUM, MPI_COMM_WORLD, &red_req);
        reduce_in_flight = true;
    };

    auto poll_reduce_nb = [&] (unsigned long long& g_sent, unsigned long long& g_recv,
                               int& any_active, int& completed) {
        completed = 0;
        if (!reduce_in_flight) return;
        int done = 0;
        MPI_Test(&red_req, &done, MPI_STATUS_IGNORE);
        if (done) {
            reduce_in_flight = false;
            g_sent = global_vec[0];
            g_recv = global_vec[1];
            any_active = (global_vec[2] > 0ULL) ? 1 : 0;
            completed = 1;
        }
    };

    // Main progress loop
    int global_continue = 1;
    int iter = 0;
    int idle_spins = 0;

    do {
        // Drain MPI receives (expands self work immediately, stages remote work)
        int rx_progress = progress_receives(rxbuf, rxreq, HIT_T, TAG_HITS, handle_batch, recv_hits);

        // Eagerly try to ship any staged remote work (non-blocking)
        pump_all_nb(out, HIT_T, TAG_HITS, sent_hits);

        // Reap completed sends
        int tx_progress = progress_sends(out);

        // Quick probe: any incoming matchable now?
        int incoming_ready = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_HITS, MPI_COMM_WORLD, &incoming_ready, MPI_STATUS_IGNORE);

        // Local activity (self work is processed immediately; only remote buffers left)
        int local_active = (rx_progress || tx_progress || have_work(out) || incoming_ready) ? 1 : 0;

        if (local_active) {
            idle_spins = 0;
        } else {
            ++idle_spins;
        }

        // Kick a non-blocking global check occasionally or when idle for a while
        if ((iter % REDUCE_INTERVAL_ITERS) == 0 || idle_spins >= IDLE_TRIGGER_ITERS) {
            kick_reduce_nb(local_active);
            fflush(stdout);
        }
        
        // Poll any in-flight non-blocking reduction
        unsigned long long g_sent = 0ULL, g_recv = 0ULL;
        int any_active = 0, red_done = 0;
        poll_reduce_nb(g_sent, g_recv, any_active, red_done);

        // If we have a fresh global view, decide whether to continue
        if (red_done) {
            global_continue = (any_active || (g_sent != g_recv)) ? 1 : 0;
        }

        ++iter;
    } while (global_continue);

    // Clean up any outstanding reduction request (should be null already, but just in case)
    if (reduce_in_flight) {
        MPI_Wait(&red_req, MPI_STATUS_IGNORE);
        reduce_in_flight = false;
    }

    // Drain & close
    for (int d = 0; d < W; ++d) {
        auto& o = out[d];
        for (int k = 0; k < SEND_DEPTH; ++k) {
            if (o.reqs[k] != MPI_REQUEST_NULL) {
                MPI_Wait(&o.reqs[k], MPI_STATUS_IGNORE);
                o.reqs[k] = MPI_REQUEST_NULL;
            }
            o.bufs[k].clear();
        }
        o.staging.clear();
        o.staging.shrink_to_fit();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ----------------------------------------------------------- Discovery -> Solving --------------------------------------------------

    /* Store which tiers were initialized and flush data to disk */
    std::vector<unsigned> initialized_tiers;
    for (int t = game::MAX_TIER; t >= 0; --t) {
        if (ctx_by_tier[t].inited) {
            initialized_tiers.push_back(t);
            pageops::flush_all(ctx_by_tier[t]);
        }
    }

    // Initialize global_db file (in data directory, not owned by any particular rank)
    // Rank 0 creates it first, then all ranks open it
    std::FILE* global_db = nullptr;
    if (R == 0) global_db = pageops::global_db_init("data");
    MPI_Barrier(MPI_COMM_WORLD); // Ensure rank 0 creates it first
    if (R != 0) global_db = pageops::global_db_init("data");

    /* Initializing main to module variables */
    int outcount, indices[BATCH_MAX]; 
    std::uint8_t shipped[BATCH_MAX]; // persisted
    MPI_Request active_reqs[BATCH_MAX]; // persisted
    for (int b = 0; b < BATCH_MAX; b++) active_reqs[b] = MPI_REQUEST_NULL;
    int forward_tier = game::MAX_TIER, backward_tier = game::MAX_TIER;

    /* Ships main to module, handles filled buffer */
    handle_batch_tier = [&](const TierMsg* msgs, int n, int recipient) {
        int completed_sends = 0, next_msg = 0;
        while (completed_sends < n) {
            MPI_Testsome(BATCH_MAX, active_reqs, &outcount, indices, MPI_STATUSES_IGNORE);

            int end = outcount == MPI_UNDEFINED ? n - completed_sends : std::min(outcount, n - completed_sends);
            for (int i = 0; i < end; i++) {
                const TierMsg& m = msgs[next_msg++];
                shipped[indices[i]] = pageops::read_solved_value(m.shape, m.page, m.bit, ctx_by_tier[m.tier]);
                MPI_Isend(shipped + indices[i], 1, MPI_UINT8_T, recipient, TAG_MAIN_TO_MODULE, MPI_COMM_WORLD, active_reqs + indices[i]);
                completed_sends++;
            }
        }
    };

    /* Outstanding Irecv buffer for module to main receives */
    MPI_Datatype TIER_T = pageops::make_TierMsg_type();
    std::vector<TierMsg> trbuf(BATCH_MAX); 
    MPI_Request trreq;

    MPI_Irecv(trbuf.data(), (int)trbuf.size(), TIER_T, MPI_ANY_SOURCE, TAG_MODULE_TO_MAIN,
              MPI_COMM_WORLD, &trreq);

    /* Buffer for skipped positions */
    std::vector<TierMsg> curr_skipped;
    
    /* Necessary variables for token passing scheme (1 rank writing to db at a time) */
    int token = (R == 0) ? 1 : 0, tiers_solved = 0, send_done, tk_flag;
    bool tk_valid = (R == 0) ? true : false;
    MPI_Request tk_rec_req = MPI_REQUEST_NULL, bc_req;
    MPI_Status tk_status;

    /* Outstanding token receive buffer */
    MPI_Irecv(&token, 1, MPI_INT, MPI_ANY_SOURCE, TAG_SEND_TOKEN, MPI_COMM_WORLD, &tk_rec_req);

    /* Main progress loop for solving */
    std::printf("Rank %d: Starting solving stage\n", R);
    do {
        
        /* Complete module to main receives */
        int rx_progress = progress_receives_tier(trbuf, trreq, TIER_T, handle_batch_tier);

        /* Progress main to module sends */
        MPI_Testall(BATCH_MAX, active_reqs, &send_done, MPI_STATUSES_IGNORE);

        /* Receive token, if pending */
        if (tk_rec_req != MPI_REQUEST_NULL) {
            MPI_Test(&tk_rec_req, &tk_flag, &tk_status);
            if (tk_flag) {
                tk_rec_req = MPI_REQUEST_NULL;
                tk_valid = true;
            }
        }

        /* Facilitate the exchange of token (forward/backward tracking) */
        if (token && tk_valid) {
            // Forward pass
            if (forward_tier == backward_tier) {
                if (ctx_by_tier[forward_tier].inited) curr_skipped = pageops::solve_page(ctx_by_tier[forward_tier], forward_tier, MPI_COMM_WORLD);
                else curr_skipped.clear();
                if (R != W - 1) MPI_Send(&token, 1, MPI_INT, R + 1, TAG_SEND_TOKEN, MPI_COMM_WORLD);
                token = R == W - 1 ? 1 : 0;
                tk_valid = false;
                forward_tier--;
            }
            // Backward pass
            if (R == W - 1 || tk_status.MPI_SOURCE == R + 1) {
                int num_skipped = (int)curr_skipped.size();
                std::vector<std::uint8_t> children(num_skipped, 0);
                std::vector<MPI_Request> child_send_requests(num_skipped);
                std::vector<MPI_Request> child_req_requests(num_skipped);
                std::vector<TierMsg> outs_skipped(num_skipped);

                for (int i = 0; i < num_skipped; i++) {
                    TierMsg t = curr_skipped[i];
                    std::uint64_t restored_hash = (t.page << pageops::PAGE_BITS) | t.bit;
                    game::Encoding pos = game::unhash(static_cast<bitops6::Bitboard>(t.shape), restored_hash);
                    const game::Encoding flipped = game::flip(pos);
                    TierMsg neighbor = pageops::make_TierMsg_from_pos(flipped);
                    int dest = owner_of_shape(neighbor.shape, W);
                    
                    if (dest == R) {
                        child_send_requests[i] = child_req_requests[i] = MPI_REQUEST_NULL;
                        uint8_t neighbor_val = read_solved_value(neighbor.shape, neighbor.page, neighbor.bit, ctx_by_tier[neighbor.tier]);
                        write_solved_value(ctx_by_tier[t.tier], t.shape, t.page, t.bit, neighbor_val + 1);
                        continue;
                    } else {
                        outs_skipped[i] = neighbor;
                        MPI_Isend(outs_skipped.data() + i, 1, TIER_T, dest, TAG_MODULE_TO_MAIN, MPI_COMM_WORLD, child_send_requests.data() + i);
                        MPI_Irecv(children.data() + i, 1, MPI_UINT8_T, dest, TAG_MAIN_TO_MODULE, MPI_COMM_WORLD, child_req_requests.data() + i);
                    }
                }
                if (num_skipped) {
                    MPI_Waitall(child_send_requests.size(), child_send_requests.data(), MPI_STATUSES_IGNORE);
                    MPI_Waitall(child_req_requests.size(), child_req_requests.data(), MPI_STATUSES_IGNORE);
                    for (int i = 0; i < num_skipped; i++) {
                        if (children[i]) {
                            TierMsg t = curr_skipped[i];
                            write_solved_value(ctx_by_tier[t.tier], t.shape, t.page, t.bit, children[i] + 1);
                        }
                    }
                }

                if (R != 0) {
                    MPI_Send(&token, 1, MPI_INT, R - 1, TAG_SEND_TOKEN, MPI_COMM_WORLD);
                    token = 0;
                } else {
                    MPI_Send(&token, 1, MPI_INT, 0, TAG_SEND_TOKEN, MPI_COMM_WORLD);
                    MPI_Irecv(&token, 1, MPI_INT, MPI_ANY_SOURCE, TAG_SEND_TOKEN, MPI_COMM_WORLD, &tk_rec_req);
                    if (backward_tier == 0) token = 0;
                }
                tk_valid = false;
                backward_tier--;
            }
        }

        if (!token && tk_rec_req == MPI_REQUEST_NULL) MPI_Irecv(&token, 1, MPI_INT, MPI_ANY_SOURCE, TAG_SEND_TOKEN, MPI_COMM_WORLD, &tk_rec_req);
        int incoming_ready = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &incoming_ready, MPI_STATUS_IGNORE);

        /* Update global termination condition */
        tiers_solved = backward_tier < 0 && !incoming_ready && !send_done && !rx_progress;
    } while (!tiers_solved);

    std::printf("Rank %d: Finished solving stage\n", R);

    /* Close all contexts after solving */
    std::printf("Rank %d: Closing all contexts\n", R);
    for (unsigned t : initialized_tiers) {
        if (ctx_by_tier[t].inited) {
            std::printf("Rank %d: Closing context for tier %u\n", R, t);
            pageops::pages_close(ctx_by_tier[t]);
        }
    }
    std::printf("Rank %d: All contexts closed\n", R);

    // Close global_db
    pageops::global_db_close(global_db);
    std::printf("Rank %d: Closed global_db\n", R);

    // Tear down MPI communication
    MPI_Cancel(&rxreq);
    MPI_Request_free(&rxreq);
    MPI_Cancel(&trreq);
    MPI_Request_free(&trreq);
    MPI_Cancel(&tk_rec_req);
    MPI_Request_free(&tk_rec_req);
    MPI_Type_free(&HIT_T);
    MPI_Type_free(&TIER_T);

    // std::printf("%llu\n", (unsigned long long)how_many_pos);

    MPI_Finalize();
    return 0;
}

// mpic++ -std=c++20 -O3 -march=native -DNDEBUG -Wall -Wextra \
//     -I/opt/homebrew/opt/libomp/include \
//     -L/opt/homebrew/opt/libomp/lib -lomp \
//     -o othello6 main.cpp
// mpirun -np 8 ./othello6