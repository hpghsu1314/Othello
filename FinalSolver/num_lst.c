#include "num_lst.h"

num_entry_t* num_head = NULL;

void num_init() {
    num_head = malloc(sizeof(num_entry_t));
}

num_entry_t* num_start() {
    return num_head;
}

void num_insert_first(uint64_t board_shape) {
    num_entry_t* new_entry = malloc(sizeof(num_entry_t));
    new_entry->board_shape = board_shape;
    new_entry->next = num_head->next;
    num_head->next = new_entry;
}

num_entry_t* num_remove_first() {
    num_entry_t* r_val = num_head->next;
    num_head->next = num_head->next->next;
    r_val->next = NULL;
    return r_val;
}

num_entry_t* num_first() {
    return num_head->next;
}