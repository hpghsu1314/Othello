#include "file_lst.h"

file_entry_t* file_head = NULL;
file_entry_t* file_tail = NULL;

void file_init() {
    file_head = malloc(sizeof(file_entry_t));
    file_tail = malloc(sizeof(file_entry_t));
}

file_entry_t* file_start() {
    return file_head;
}

file_entry_t* file_end() {
    return file_tail;
}

file_entry_t* file_list_front() {
    return file_head->next;
}

file_entry_t* file_list_end() {
    return file_tail->prev;
}

bool file_is_end(file_entry_t* file_entry) {
    return file_entry->board_shape == 0;
}

FILE* find_file(uint64_t board_shape) {
    for (file_entry_t* elem = file_list_front(); !file_is_end(elem); elem = elem->next) {
        if (elem->board_shape == board_shape) {
            return elem->file;
        }
    }

    return NULL;
}

FILE* create_file(uint64_t board_shape) {
    char filename[100];
    snprintf(filename, sizeof(filename), "%s/%llu.bin", FOLDER, board_shape);
    FILE* file = fopen(filename, "wb+");
    file_entry_t* new_elem = malloc(sizeof(file_entry_t));

    file_tail->prev->next = new_elem;
    new_elem->prev = file_tail->prev;
    new_elem->next = file_tail;
    file_tail->prev = new_elem;
    new_elem->board_shape = board_shape;
    new_elem->file = file;

    return file;
}

FILE* file_open_no_trunc(uint64_t board_shape) {
    char filename[100];
    snprintf(filename, sizeof(filename), "%s/%llu.bin", FOLDER, board_shape);
    FILE* file = fopen(filename, "r+b");
    file_entry_t* new_elem = malloc(sizeof(file_entry_t));

    file_tail->prev->next = new_elem;
    new_elem->prev = file_tail->prev;
    new_elem->next = file_tail;
    file_tail->prev = new_elem;
    new_elem->board_shape = board_shape;
    new_elem->file = file;

    return file;
}

void close_first() {
    file_entry_t* temp = file_head->next;
    file_head->next = temp->next;
    temp->next->prev = file_head;
    fclose(temp->file);
    free(temp);
}

void file_remove_elem(file_entry_t* elem) {
    elem->prev->next = elem->next;
    elem->next->prev = elem->prev;
    fclose(elem->file);
    free(elem);
}

bool file_is_empty() {
    return file_head->next->board_shape == 0;
}

long get_file_size(FILE* f) {
    fseek(f, 0, SEEK_END);
    return ftell(f);
}