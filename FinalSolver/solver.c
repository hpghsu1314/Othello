#include "solver.h"

int main() {
    struct Encoding start;

    file_init();
    num_init();

    file_entry_t* start_of_file = file_start();
    file_entry_t* end_of_file = file_end();
    end_of_file->board_shape = 0;
    start_of_file->board_shape = 0;
    start_of_file->next = end_of_file;
    end_of_file->prev = start_of_file;
    start_of_file->prev = NULL;
    end_of_file->next = NULL;

    num_entry_t* start_of_num = num_start();
    start_of_num->board_shape = 0;
    start_of_num->next = NULL;

    // Player is the bits in the person's perspective
    // i.e. there is no "black" or "white" positions, but only "mine" and "yours"
    // This design is due to possible positions having "no moves" and having to pass onto the other player

    // Print with %llu
    start.player   = 0b0000001001000000;
    start.opponent = 0b0000010000100000;

    start = symmetry(start);

    write_bits_w_hash(start.player | start.opponent, FOUND, hash(start));

    file_entry_t* open_file;
    unsigned char buffer[4096];
    char byte;
    size_t read_bytes;
    uint64_t hash_val;
    struct Encoding temp;
    char filename[100];
    uint64_t temp_board;
    uint64_t c = 0;

    printf("Starting\n");
    // Discovery
    while (!file_is_empty()) {
        open_file = file_list_front();
        fseek(open_file->file, 0, SEEK_SET);
        num_insert_first(open_file->board_shape);
        hash_val = 0;
        while (!feof(open_file->file)) {
            read_bytes = fread(buffer, sizeof(unsigned char), 4096, open_file->file);
            if (read_bytes == 0 && feof(open_file->file)) {
                break;
            } else if (read_bytes != 0) {
                for (uint64_t index = 0; index < read_bytes; index++) {
                    if (buffer[index] == FOUND) {
                        do_moves_and_write(unhash(open_file->board_shape, hash_val));
                    }
                    hash_val++;
                }
            }
            fseek(open_file->file, hash_val, SEEK_SET);
        }
        file_remove_elem(open_file);
    }
    
    num_entry_t* t_num;
    uint8_t prim_val;
    
    // Primitive + Remoteness
    while (start_of_num->next) {
        t_num = num_remove_first();
        FILE* f = file_open_no_trunc(t_num->board_shape);
        fseek(f, 0, SEEK_SET);
        hash_val = 0;
        while (!feof(f)) {
            read_bytes = fread(buffer, sizeof(unsigned char), 4096, f);
            if (read_bytes == 0 && feof(f)) {
                break;
            } else {
                for (uint64_t index = 0; index < read_bytes; index++) {
                    if (buffer[index] == FOUND) {
                        prim_val = primitive(unhash(t_num->board_shape, hash_val));
                        write_bits_w_hash(t_num->board_shape, prim_val, hash_val);
                    }
                    hash_val++;
                }
            }
            fseek(f, hash_val, SEEK_SET);
        }
    }

    struct Metadata header;
    uint8_t prev;
    file_entry_t* file_to_remove;

    prev = file_read_board(start);

    // Condensing
    for (file_entry_t* e = file_list_end(); e->board_shape != 0;) {
        header.board_shape = e->board_shape;
        header.size = get_file_size(e->file);
        fseek(e->file, 0, SEEK_SET);
        
        while (!feof(e->file)) {
            long pos = ftell(e->file);
            read_bytes = fread(buffer, sizeof(unsigned char), 4096, e->file);
            if (read_bytes == 0) {
                if (feof(e->file)) {
                    break;
                }
            } else {
                for (uint64_t index = 0; index < read_bytes; index++) {
                    if (buffer[index] == SKIP) {
                        buffer[index] = prev;
                    } else {
                        prev = buffer[index];
                    }
                }
                fseek(e->file, pos, SEEK_SET);
                fwrite(buffer, sizeof(unsigned char), read_bytes, e->file);
                fseek(e->file, pos + read_bytes, SEEK_SET);
            }
        }
        file_to_remove = e;
        e = e->prev;
        snprintf(filename, sizeof(filename), "%s/%llu.bin", FOLDER, file_to_remove->board_shape);
        file_remove_elem(file_to_remove);
    }
    
    free(start_of_file);
    free(end_of_file);
    free(start_of_num);

    return 0;
}

