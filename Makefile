CC = gcc
MPICC = mpicc

ifeq ($(DEBUG),1)
	CFLAGS = -g -O0 -Wall -Wextra
else
	CFLAGS = -O3 -Wall -Wextra
endif
$(info [INFO] Building with CFLAGS = $(CFLAGS))

PROJECT_ROOT = .
NAIVE_DIR = FinalSolver
MPI_DIR = MPISolver

NAIVE_SRCS = $(wildcard $(NAIVE_DIR)/*.c)
MPI_SRCS = $(wildcard $(MPI_DIR)/*.c)

NAIVE_HDRS = $(wildcard $(NAIVE_DIR)/*.h)
MPI_HDRS = $(wildcard $(MPI_DIR)/*.h)

NAIVE_BIN = final_solver
MPI_BIN = solver_mpi

all: $(if $(strip $(NAIVE_SRCS)),$(NAIVE_BIN)) $(if $(strip $(MPI_SRCS)),$(MPI_BIN))

clean:
	rm -f $(NAIVE_BIN) $(MPI_BIN)

ifneq ($(strip $(NAIVE_SRCS)),)
$(NAIVE_BIN): $(NAIVE_SRCS) $(NAIVE_HDRS)
	$(CC) $(CFLAGS) -o $@ $(NAIVE_SRCS)
endif

ifneq ($(strip $(MPI_SRCS)),)
$(MPI_BIN): $(MPI_SRCS) $(MPI_HDRS)
	$(MPICC) $(CFLAGS) -o $@ $(MPI_SRCS)
endif

compile-naive: $(NAIVE_BIN)
compile-src:   $(MPI_BIN)

DOCKER_CONTAINER_NAME = othello

docker-build:
	docker build -t $(DOCKER_CONTAINER_NAME) $(PROJECT_ROOT)

docker-run:
	docker run --rm -it \
	-v "$(shell pwd):/app" \
	-w /app \
	$(DOCKER_CONTAINER_NAME) bash

.PHONY: all clean compile-naive compile-src docker-build docker-run
