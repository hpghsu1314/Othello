#ifndef FILE_LST_H
#define FILE_LST_H

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define FOLDER "Data4x4" // Remember to change this for other games
#define COMBINE_FILE "Combined4x4.bin"

typedef struct file_entry {
    uint64_t board_shape;
    FILE* file;
    struct file_entry* next;
    struct file_entry* prev;
} file_entry_t;

extern file_entry_t* file_head;
extern file_entry_t* file_tail;

void file_init();
file_entry_t* file_start();
file_entry_t* file_end();
file_entry_t* file_list_front();
file_entry_t* file_list_end();
bool file_is_end(file_entry_t* file_entry);
FILE* find_file(uint64_t board_shape);
FILE* create_file(uint64_t board_shape);
FILE* file_open_no_trunc(uint64_t board_shape);
void close_first();
void file_remove_elem(file_entry_t* elem);
bool file_is_empty();
long get_file_size(FILE* f);

#endif