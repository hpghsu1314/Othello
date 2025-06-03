FROM gcc:latest

# Installs OpenMPI
RUN apt-get update && apt-get install -y --no-install-recommends \
    openmpi-bin libopenmpi-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . .

RUN make

CMD ["mpirun", "-np", "4", "./MPISolver/solver_mpi"]