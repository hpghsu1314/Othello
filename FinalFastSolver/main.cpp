// Minimal MPI + OpenMP bootstrap (no game code)
// build: mpicxx -O3 -fopenmp main.cpp -o solver

#include <mpi.h>
#include <omp.h>
#include <cstdint>
#include <cstdio>
#include "process.hpp"
#include "othello4.hpp"

int main(int argc, char** argv) {
    // Start MPI with threading (OpenMP inside each rank; MPI calls on master)
    MPI_Init_thread(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    
    MPI_Finalize();
    return 0;
}
