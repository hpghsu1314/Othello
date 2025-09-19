// message.hpp
#pragma once
#include <cstdint>
#include <mpi.h>

#pragma pack(push,1)
struct HitMsg {
    std::uint64_t shape;
    std::uint32_t page;
    std::uint32_t bit;
};
#pragma pack(pop)

inline MPI_Datatype make_HitMsg_type() {
    MPI_Datatype T;
    const int      N   = 3;
    int            bl[N] = {1,1,1};
    MPI_Aint       dis[N];
    MPI_Datatype   ty[N] = {MPI_UINT64_T, MPI_UINT32_T, MPI_UINT32_T};

    HitMsg probe{};
    MPI_Aint base = 0, a0=0, a1=0, a2=0;
    MPI_Get_address(&probe, &base);
    MPI_Get_address(&probe.shape, &a0);
    MPI_Get_address(&probe.page,  &a1);
    MPI_Get_address(&probe.bit,   &a2);
    dis[0] = a0 - base; dis[1] = a1 - base; dis[2] = a2 - base;

    MPI_Type_create_struct(N, bl, dis, ty, &T);
    MPI_Type_commit(&T);
    return T;
}
