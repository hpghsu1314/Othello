#ifndef SOLVER
#define SOLVER

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "game.h"

#define NUM_WRITERS 9

struct Metadata {
    uint64_t board_shape;
    uint64_t size;
};

#endif