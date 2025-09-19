#include <mpi.h>
#include <vector>
#include <queue>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <array>
#include "utilities.hpp"
#include "page_ops.hpp"
#include "message.hpp"
#include "othello4.hpp"

// Tunables
static constexpr int BATCH_MAX = 4096;
static constexpr int SEND_DEPTH = 4;
static constexpr int TAG_HITS  = 1;
static constexpr int TAG_EXIT  = 9;

struct Outgoing {
    std::array<std::vector<HitMsg>, SEND_DEPTH> bufs;
    std::array<MPI_Request, SEND_DEPTH> reqs;
    int cur = 0;
};

int main(int argc, char** argv) {
    int provided=0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int W=0,R=0;
    MPI_Comm_size(MPI_COMM_WORLD, &W);
    MPI_Comm_rank(MPI_COMM_WORLD, &R);

    MPI_Datatype HIT_T = make_HitMsg_type();

    std::vector<Outgoing> out(W);
    for (int d=0; d<W; ++d) {
        for (int k=0; k<SEND_DEPTH; ++k) {
            out[d].bufs[k].reserve(BATCH_MAX);
            out[d].reqs[k] = MPI_REQUEST_NULL;
        }
        out[d].cur = 0;
    }

    std::vector<HitMsg> rxbuf(BATCH_MAX);
    MPI_Request rxreq = MPI_REQUEST_NULL;
    MPI_Irecv(rxbuf.data(), BATCH_MAX, HIT_T, MPI_ANY_SOURCE, TAG_HITS, MPI_COMM_WORLD, &rxreq);

    HitMsg m{};
    if (R == 0) {
        Game G;
        auto S0 = G.starting_position();
        auto Sc = G.canonical(S0);
        uint64_t sh = G.shape(Sc);
        uint64_t h  = G.hash(Sc);
        m.shape = sh;
        m.page = (uint32_t)(h >> pageops::PAGE_BITS);
        m.bit = (uint32_t)(h & ((1u<<pageops::PAGE_BITS)-1));
    }
    MPI_Bcast(&m, 1, HIT_T, 0, MPI_COMM_WORLD);

    int dest = owner_of_shape(m.shape, W);
    if (R == dest) {
        pageops::set_bit(ctx_t0, m.shape, m.page, m.bit);
    }



    MPI_Cancel(&rxreq);
    MPI_Request_free(&rxreq);

    MPI_Type_free(&HIT_T);
    MPI_Finalize();
    return 0;
}
