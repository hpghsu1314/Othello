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

#include "utilities.hpp"  // owner_of_shape(...)
#include "page_ops.hpp"   // pageops::Context, pages_init, set_bit, ...
#include "othello6.hpp"   // swap to othello4.hpp when running 4x4

// Throughput knobs (correctness does not depend on specific values)
static constexpr int BATCH_MAX  = 4096; // send/recv chunk size
static constexpr int SEND_DEPTH = 4;    // in-flight slots per destination
static constexpr int TAG_HITS   = 1;

#pragma pack(push,1)
struct HitMsg {
    std::uint64_t shape;
    std::uint32_t page;
    std::uint32_t bit;
};
#pragma pack(pop)

struct Outgoing {
    // In-flight send slots (never modify a buf while its req != MPI_REQUEST_NULL)
    std::array<std::vector<HitMsg>, SEND_DEPTH> bufs;
    std::array<MPI_Request, SEND_DEPTH> reqs;
    // Small staging per destination; non-blocking pumps keep it small
    std::vector<HitMsg> staging;
};

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

// ---- Reap one free slot non-blocking; return slot index or -1 ----
static inline int reap_one_free_slot_nb(Outgoing& o) {
    for (int k = 0; k < SEND_DEPTH; ++k) {
        if (o.reqs[k] == MPI_REQUEST_NULL && o.bufs[k].empty()) return k;
    }
    int index = MPI_UNDEFINED, flag = 0;
    MPI_Testany(SEND_DEPTH, o.reqs.data(), &index, &flag, MPI_STATUS_IGNORE);
    if (flag && index != MPI_UNDEFINED) {
        o.reqs[index] = MPI_REQUEST_NULL;
        o.bufs[index].clear();
        return index;
    }
    return -1;
}

// ---- Try to send one chunk (non-blocking). If no slot, just return. ----
static inline void try_send_dest_nb(Outgoing& o, int dest,
                                    MPI_Datatype HIT_T, int tag,
                                    unsigned long long& sent_hits)
{
    if (o.staging.empty()) return;

    int slot = reap_one_free_slot_nb(o);
    if (slot < 0) return;

    const int n = (int)std::min<std::size_t>(o.staging.size(), (std::size_t)BATCH_MAX);
    o.bufs[slot].resize(n);
    std::memcpy(o.bufs[slot].data(), o.staging.data(), n * sizeof(HitMsg));

    if ((int)o.staging.size() == n) o.staging.clear();
    else o.staging.erase(o.staging.begin(), o.staging.begin() + n);

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
        if (!ctx.inited) pageops::pages_init(ctx, R, tier, "data");
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
    // IMPORTANT CHANGE: we NEVER call MPI here anymore. We only push remote
    // children into per-destination staging; the main loop ships them.
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

            const std::uint64_t h =
                ((std::uint64_t)msg.page << pageops::PAGE_BITS) |
                 (std::uint64_t)msg.bit;
            game::Encoding pos = game::unhash(msg.shape, h);

            std::uint64_t moves = game::legal_moves(pos);
            if (moves) {
                while (moves) {
                    std::uint64_t mv = moves & -moves;
                    moves ^= mv;
                    const game::Encoding nxt = game::flip(game::do_move(pos, mv));
                    const HitMsg m2 = make_HitMsg_from_pos(nxt);
                    const int dest = owner_of_shape(m2.shape, W);
                    if (dest == R) {
                        stack.emplace_back(m2); // keep processing locally
                    } else {
                        out[dest].staging.emplace_back(m2); // NO MPI here
                    }
                }
            } else {
                // pass move
                const HitMsg m2 = make_HitMsg_from_pos(game::flip(pos));
                const int dest = owner_of_shape(m2.shape, W);
                if (dest == R) {
                    stack.emplace_back(m2);
                } else {
                    out[dest].staging.emplace_back(m2); // NO MPI here
                }
            }
        }
    };

    // Batch handler: route each message (self -> process now; remote children are staged)
    handle_batch = [&](const HitMsg* msgs, int n) {
        for (int i = 0; i < n; ++i) {
            const HitMsg& m = msgs[i];
            const int dest = owner_of_shape(m.shape, W);
            if (dest == R) {
                process_self_iter(m);
            } else {
                out[dest].staging.emplace_back(m); // NO MPI here
            }
        }
    };

    // Seed from rank 0
    if (R == 0) {
        const HitMsg m0 = make_HitMsg_from_pos(game::starting_position());
        const int dest0 = owner_of_shape(m0.shape, W);
        if (dest0 == R) {
            process_self_iter(m0);
        } else {
            out[dest0].staging.emplace_back(m0); // NO MPI here
        }
    }

    // Main progress loop
    int global_continue = 1;
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

        // Global in-flight credits
        unsigned long long g_sent = 0ULL, g_recv = 0ULL;
        MPI_Allreduce(&sent_hits, &g_sent, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&recv_hits, &g_recv, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        // Local activity (self work is processed immediately; only remote buffers left)
        int local_active = (rx_progress || tx_progress || have_work(out) || incoming_ready) ? 1 : 0;
        int any_active = 0;
        MPI_Allreduce(&local_active, &any_active, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

        global_continue = (any_active || (g_sent != g_recv)) ? 1 : 0;
    } while (global_continue);

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

    for (unsigned t = 0; t <= game::MAX_TIER; ++t) {
        if (ctx_by_tier[t].inited) {
            pageops::flush_all(ctx_by_tier[t]);
            pageops::pages_close(ctx_by_tier[t]);
        }
    }

    std::printf("%llu\n", (unsigned long long)how_many_pos);

    // Tear down receive
    MPI_Cancel(&rxreq);
    MPI_Request_free(&rxreq);
    MPI_Type_free(&HIT_T);
    MPI_Finalize();
    return 0;
}
