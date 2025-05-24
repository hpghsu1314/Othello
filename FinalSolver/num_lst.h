#ifndef NUM_LST
#define NUM_LST

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct num_entry {
    uint64_t board_shape;
    struct num_entry* next;
} num_entry_t;

extern num_entry_t* num_head;

void num_init();
num_entry_t* num_start();
void num_insert_first(uint64_t board_shape);
num_entry_t* num_remove_first();
num_entry_t* num_first();

#endif