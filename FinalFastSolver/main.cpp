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

            const std::uint64_t h =
                ((std::uint64_t)msg.page << pageops::PAGE_BITS) |
                 (std::uint64_t)msg.bit;
            game::Encoding pos = time_block("unhash", [&] {
                return game::unhash(msg.shape, h);
            });

            std::uint64_t moves = time_block("legal_moves", [&] {
                return game::legal_moves(pos);
            });

            if (moves && !pageops::set_bit(ctx, msg.shape, msg.page, msg.bit)) {
                continue;
            } else if (game::primitive(pos)) {
                printf("discovery side: player: %llu, opponent: %llu \n", pos.player, pos.opponent);
                pageops::set_bit(ctx, msg.shape, msg.page, msg.bit);
                continue;
            }

            ++how_many_pos;
            if ((how_many_pos % 100000) == 0) {
                // printf("Rank %d has found %llu\n", R, how_many_pos);
                fflush(stdout);
                TimerAggregator::instance().report();
            }

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
    // game::Encoding test = {
    //     0b0000000000000000000000000000111111100001100001100001100001111110,  // Player (P)
    //     0b0000000000000000000000000000000000011110011110011110011110000000  // player (O)
    // };

    game::Encoding test = {
        0b0000000000000000000000000000011111100001100001100001100001111110,  // Player (P)
        0b0000000000000000000000000000000000011010001110011010011110000000  // player (O)
    };


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

    // stores offsets of page records per rank per tier
    int global_tier = game::MAX_TIER;
    std::vector<std::vector<uint64_t>> data_offsets(game::MAX_TIER + 1, std::vector<uint64_t>(W, 0));
    std::vector<std::vector<uint64_t>> rec_offsets(game::MAX_TIER + 1, std::vector<uint64_t>(W, 0));
    std::vector<uint64_t> data_flat((game::MAX_TIER + 1) * W, 0);
    std::vector<uint64_t> rec_flat((game::MAX_TIER + 1) * W, 0);

    // Utility function for getting correct offsets
    auto get_next_offset = [&](std::vector<std::uint64_t> vec, std::uint64_t curr) -> std::uint64_t {
        std::vector<std::uint64_t> sorted = vec;
        std::sort(sorted.begin(), sorted.end());
        for (std::uint64_t off : sorted) {if (off > curr) return off;}
        return 0;
    };
    
    // Safe file creation 
    int mode = MPI_MODE_WRONLY | MPI_MODE_CREATE;

    MPI_File tier_db  = pageops::init_db_write_only("data", "dat", global_tier, MPI_COMM_WORLD, mode);
    MPI_File tier_rec = pageops::init_db_write_only("data", "idx", global_tier, MPI_COMM_WORLD, mode);
    MPI_Barrier(MPI_COMM_WORLD);

    /* Main progress loop for solving */
    std::printf("Rank %d: Starting solving stage\n", R);
    do {
        // Solve tier 
        solve_tier(ctx_by_tier[global_tier], global_tier, R, W, tier_db, tier_rec, data_offsets, rec_offsets, MPI_COMM_WORLD);
        pageops::db_close(tier_db); pageops::db_close(tier_rec);

        // Termination condition
        global_tier--;
        if (global_tier < 4) break;

        // Safely initialize a global file for the next tier 
        MPI_Barrier(MPI_COMM_WORLD);
        // Flatten offset array
        for (int i = 0; i < game::MAX_TIER + 1; ++i) {
            std::memcpy(&data_flat[i*W], data_offsets[i].data(), W * sizeof(uint64_t));
            std::memcpy(&rec_flat[i*W], rec_offsets[i].data(), W * sizeof(uint64_t));
        }
        
        MPI_Allreduce(MPI_IN_PLACE, data_flat.data(), (game::MAX_TIER + 1) * W, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, rec_flat.data(), (game::MAX_TIER + 1) * W, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
        
        // Unflatten offset array
        for (int i = 0; i < game::MAX_TIER + 1; ++i) {
            std::memcpy(data_offsets[i].data(), &data_flat[i*W], W * sizeof(uint64_t));
            std::memcpy(rec_offsets[i].data(), &rec_flat[i*W], W * sizeof(uint64_t));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Compression 
        if (global_tier <= game::MAX_TIER - 2) {
            int compressable_tier = global_tier + 2;
            std::uint64_t curr_off = rec_offsets[compressable_tier][R];
            std::uint64_t next_off = get_next_offset(rec_offsets[compressable_tier], curr_off);
            pageops::compress_all_pages(compressable_tier, R, data_offsets[compressable_tier][R], curr_off, next_off);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        tier_db  = pageops::init_db_write_only("data", "dat", global_tier, MPI_COMM_WORLD, mode);
        tier_rec = pageops::init_db_write_only("data", "idx", global_tier, MPI_COMM_WORLD, mode);
    } while (true);

    std::uint64_t curr_off = rec_offsets[5][R];
    std::uint64_t next_off = get_next_offset(rec_offsets[5], curr_off);
    pageops::compress_all_pages(5, R, data_offsets[5][R], curr_off, next_off);

    curr_off = rec_offsets[4][R];
    next_off = get_next_offset(rec_offsets[4], curr_off);
    pageops::compress_all_pages(4, R, data_offsets[4][R], curr_off, next_off);

    std::printf("Rank %d: Finished solving stage\n", R);

    // Close all contexts after solving 
    std::printf("Rank %d: Closing all contexts\n", R);
    for (unsigned t : initialized_tiers) {
        if (ctx_by_tier[t].inited) {
            std::printf("Rank %d: Closing context for tier %u\n", R, t);
            pageops::pages_close(ctx_by_tier[t]);
        }
    }
    std::printf("Rank %d: All contexts closed\n", R);

    // Tear down MPI communication
    MPI_Cancel(&rxreq);
    MPI_Request_free(&rxreq);
    MPI_Type_free(&HIT_T);

    std::printf("number of positions: %llu\n", (unsigned long long)how_many_pos);

    MPI_Finalize();
    return 0;
}

// mpic++ -std=c++20 -O3 -march=native -DNDEBUG -Wall -Wextra \
//     -I/opt/homebrew/opt/libomp/include \
//     -L/opt/homebrew/opt/libomp/lib -lomp \
//     -I/opt/homebrew/include \
//     -L/opt/homebrew/lib -lzstd \
//     -o othello6 main.cpp
// mpirun -np 8 ./othello6