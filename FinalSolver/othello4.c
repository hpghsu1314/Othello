#include "othello4.h"

uint64_t extra = 0;

// Change symmetries and stuff

// Remove all of the write operations, put them outside of this file

struct Encoding symmetries[7];

int directions[8] = {-BOARD_SIZE + 1,  1, BOARD_SIZE + 1, -BOARD_SIZE - 1, -1, BOARD_SIZE - 1, BOARD_SIZE, -BOARD_SIZE};

Bitboard valid_move(struct Encoding board) {

    Bitboard player = board.player;
    Bitboard opponent = board.opponent;
    Bitboard empty = FULL_BOARD & (~ (player | opponent)); // Empty squares
    Bitboard valid = 0;

    Bitboard left = player;
    Bitboard right = player;
    Bitboard up = player;
    Bitboard down = player;
    Bitboard l_up = player;
    Bitboard l_down = player;
    Bitboard r_up = player;
    Bitboard r_down = player;

    for (int i = 0; i < BOARD_SIZE - 1; i++) {

        //Checking Left
        left = left & NOT_A_FILE;
        left = (left << 1) & opponent;
        left = left & NOT_A_FILE;
        valid |= (left << 1) & empty;
        
        //Checking Right
        right = right & NOT_F_FILE;
        right = (right >> 1) & opponent;
        right = right & NOT_F_FILE;
        valid |= (right >> 1) & empty;

        //Checking Up
        up = (up << BOARD_SIZE) & opponent;
        valid |= (up << BOARD_SIZE) & empty;

        //Checking Down
        down = (down >> BOARD_SIZE) & opponent;
        valid |= (down >> BOARD_SIZE) & empty;

        //Checking Left-Up
        l_up = l_up & NOT_A_FILE;
        l_up = (l_up << BOARD_SIZE + 1) & opponent;
        l_up = l_up & NOT_A_FILE;
        valid |= (l_up << BOARD_SIZE + 1) & empty;

        //Checking Left-Down
        l_down = l_down & NOT_F_FILE;
        l_down = (l_down >> BOARD_SIZE + 1) & opponent;
        l_down = l_down & NOT_F_FILE;
        valid |= (l_down >> BOARD_SIZE + 1) & empty;

        //Checking Right-Up
        r_up = r_up & NOT_F_FILE;
        r_up = (r_up << BOARD_SIZE - 1) & opponent;
        r_up = r_up & NOT_F_FILE;
        valid |= (r_up << BOARD_SIZE - 1) & empty;

        //Checking Right-Down
        r_down = r_down & NOT_A_FILE;
        r_down = (r_down >> BOARD_SIZE - 1) & opponent;
        r_down = r_down & NOT_A_FILE;
        valid |= (r_down >> BOARD_SIZE - 1) & empty;
    }
    return valid & FULL_BOARD;
}

struct Encoding do_move(struct Encoding board, uint64_t move) {

    Bitboard player = board.player;
    Bitboard opponent = board.opponent;

    struct Encoding r_val;
    uint64_t curr;
    int dir;
    uint64_t to_remove;

    for (int i = 0; i < 8; i++) {
        dir = directions[i];
        curr = move;
        to_remove = 0;
        while (curr > (uint64_t) 0) {
            if (((dir + 2*BOARD_SIZE) % BOARD_SIZE) == 1) {
                curr = curr & NOT_F_FILE;
            } else if (((dir + 2*BOARD_SIZE) % BOARD_SIZE) == (BOARD_SIZE - 1)) {
                curr = curr & NOT_A_FILE;
            }
            curr = (dir < 0) ? (curr << -dir) : (curr >> dir);
            if ((curr & player) > 0) {
                player |= to_remove;
                opponent = opponent & (~to_remove);
                break;
            } else if ((curr & opponent) > 0) {
                to_remove |= curr;
            } else {
                break;
            }
        }
    }

    r_val.player = player | move;
    r_val.opponent = opponent;

    return r_val;
}

// Transpose a 4x4 board
Bitboard transpose4(Bitboard board) {
    Bitboard res = 0;
    uint64_t diag = 0b1000010000100001;
    uint64_t r1 = 0b0100001000010000;
    uint64_t r2 = 0b0010000100000000;
    uint64_t r3 = 0b0001000000000000;
    uint64_t l1 = (r1 >> 3);
    uint64_t l2 = (r2 >> 6);
    uint64_t l3 = (r3 >> 9);
    res |= (r1 & board) >> 3;
    res |= (r2 & board) >> 6;
    res |= (r3 & board) >> 9;
    res |= (l1 & board) << 3;
    res |= (l2 & board) << 6;
    res |= (l3 & board) << 9;
    res |= (diag & board);
    return res;
}

// Vertical flip a 4x4 board
Bitboard vertical4(Bitboard board) {
    Bitboard res = 0;
    uint64_t r1 = 0b0100010001000100;
    uint64_t r2 = 0b1000100010001000;
    uint64_t l1 = (r1 >> 1);
    uint64_t l2 = (r2 >> 3);
    res |= ((r1 & board) >> 1);
    res |= ((r2 & board) >> 3);
    res |= ((l1 & board) << 1);
    res |= ((l2 & board) << 3);
    return res;
}

// Horizontal flip a 4x4 board
Bitboard horizontal4(Bitboard board) {
    Bitboard res = 0;
    uint64_t r2 = 0b1111000000000000;
    uint64_t r1 = 0b0000111100000000;
    uint64_t l2 = (r2 >> 12);
    uint64_t l1 = (r1 >> 4);
    
    res |= ((r2 & board) >> 12);
    res |= ((r1 & board) >> 4);
    res |= ((l2 & board) << 12);
    res |= ((l1 & board) << 4);

    return res;
}

struct Encoding symmetry(struct Encoding board) {
    Bitboard player = board.player;
    Bitboard opponent = board.opponent;
    Bitboard v_player = vertical4(player);
    Bitboard v_opponent = vertical4(opponent);
    Bitboard t_player = transpose4(player);
    Bitboard t_opponent = transpose4(opponent);
    Bitboard hv_player = horizontal4(v_player);
    Bitboard hv_opponent = horizontal4(v_opponent);
    
    symmetries[0].player = v_player;
    symmetries[0].opponent = v_opponent;
    symmetries[1].player = horizontal4(player);
    symmetries[1].opponent = horizontal4(opponent);
    symmetries[2].player = t_player;
    symmetries[2].opponent = t_opponent;
    symmetries[3].player = vertical4(t_player);
    symmetries[3].opponent = vertical4(t_opponent);
    symmetries[4].player = hv_player;
    symmetries[4].opponent = hv_opponent;
    symmetries[5].player = horizontal4(t_player);
    symmetries[5].opponent = horizontal4(t_opponent);
    symmetries[6].player = transpose4(hv_player);
    symmetries[6].opponent = transpose4(hv_opponent);

    struct Encoding r_val;
    r_val.player = player;
    r_val.opponent = opponent;
    struct Encoding t;
    Bitboard curr = player | opponent;
    Bitboard new;

    for (int i = 0; i < 7; i++) {
        t = symmetries[i];
        new = t.player | t.opponent;
        if (curr > new) {
            curr = new;
            r_val.player = t.player;
            r_val.opponent = t.opponent;
        } else if (curr == new) {
            if (r_val.player > t.player) {
                r_val.player = t.player;
                r_val.opponent = t.opponent;
            }
        }
    }

    return r_val;
}

uint64_t hash(struct Encoding b) {
    Bitboard player = b.player;
    Bitboard board_shape = player | b.opponent;
    uint64_t hash = 0;
    uint64_t pos = 0;

    while (board_shape > 0) {
        if (board_shape & 1) {
            hash |= ((player & 1) << pos);
            pos++;
        }
        board_shape >>= 1;
        player >>= 1;
    }
    return hash;
}

struct Encoding unhash(uint64_t board_shape, uint64_t hash_value) {
    uint64_t player = 0;
    uint64_t opponent = 0;

    uint64_t pos = 0;
    while (board_shape > 0) {
        if (board_shape & 1) {
            if (hash_value & 1) {
                player |= (1 << pos);
            } else {
                opponent |= (1 << pos);
            }
            hash_value = hash_value >> 1;
        }
        board_shape = board_shape >> 1;
        pos++;
    }
    
    struct Encoding ret_val;
    ret_val.opponent = opponent;
    ret_val.player = player;

    return ret_val;
}

// uint64_t prim_helper(uint64_t child) {
//     if (child == WIN) {
//         return LOSE;
//     } else if (child == LOSE) {
//         return WIN;
//     }
//     return DRAW;
// }

uint8_t prim_helper(uint8_t child) {
    uint8_t first_bits = 0b11000000;
    uint8_t last_bits = 0b00111111;
    if ((child & first_bits) == WIN) {
        return LOSE | (child & last_bits);
    } else if ((child & first_bits) == LOSE) {
        return WIN | (child & last_bits);
    }
    return DRAW | (child & last_bits);
}

uint8_t file_read_board(struct Encoding board) {
    uint64_t board_shape = board.player | board.opponent;
    FILE* f = find_file(board_shape);
    if (f == NULL) {
        f = file_open_no_trunc(board_shape);
    }
    uint64_t hash_num = hash(board);

    char buffer = 0;
    fseek(f, hash_num, SEEK_SET);
    fread(&buffer, sizeof(char), 1, f);
    return buffer;
}

uint8_t primitive(struct Encoding board) {
    uint64_t valid = valid_move(board);

    if (valid == 0) {
        struct Encoding flip_board = symmetry(flip(board));
        if (valid_move(flip_board) == 0) {
            int curr = count_bits(board.player);
            int oppo = count_bits(board.opponent);
            if (hash(board) > hash(flip_board)) {
                return SKIP;
            } else {
                if (curr > oppo) {
                    return WIN;
                } else if (curr == oppo) {
                    return DRAW;
                } else {
                    return LOSE;
                }
            }
        } else {
            return prim_helper(primitive(flip_board)) + 1;
        }
    }

    uint64_t pos = 0;
    uint8_t curr = LOSE;
    uint8_t temp;
    uint8_t first_bits = 0b11000000;
    uint8_t last_bits = 0b00111111;

    struct Encoding next_board;

    while (valid > 0) {
        if (valid & 1) {
            next_board = symmetry(flip(do_move(board, (1 << pos))));
            uint8_t next_move = file_read_board(next_board);
            if (next_move == SKIP) {
                next_move = file_read_board(symmetry(flip(next_board)));
                if (next_move == SKIP) {
                    printf("FAIL\n");
                }
            }
            temp = prim_helper(next_move) + 1;
            if ((temp & first_bits) != WIN) {
                curr = temp > curr ? temp : curr;
            } else {
                if ((curr & first_bits) != WIN) {
                    curr = temp;
                } else {
                    curr = temp < curr ? temp : curr;
                }
            }
        }
        valid = valid >> 1;
        pos++;
    }

    return curr;
}

void write_bits_w_hash(uint64_t board_shape, uint8_t value, uint64_t hash_value) {

    FILE* file_to_write = find_file(board_shape);
    if (file_to_write == NULL) {
        file_to_write = create_file(board_shape);
    }

    fseek(file_to_write, hash_value, SEEK_SET);
    fwrite(&value, sizeof(uint8_t), 1, file_to_write);
}

void do_moves_and_write(struct Encoding initial_board) {
    uint64_t valid = valid_move(initial_board);
    uint64_t flip_hash;
    uint64_t init_hash;
    struct Encoding flip_board;

    if (valid == 0) {
        flip_board = symmetry(flip(initial_board));
        flip_hash = hash(flip_board);
        init_hash = hash(initial_board);
        
        valid = valid_move(flip_board);
        if (valid == 0) {
            if (flip_hash < init_hash) {
                write_bits_w_hash(flip_board.opponent | flip_board.player, FOUND, flip_hash);
            } else {
                write_bits_w_hash(initial_board.opponent | initial_board.player, FOUND, init_hash);
            }
            return;
        } else {
            write_bits_w_hash(flip_board.opponent | flip_board.player, FOUND, flip_hash);
            initial_board = flip_board;
        }
    }

    uint64_t pos = 0;

    struct Encoding next_board;

    while (valid > 0) {
        if (valid & 1) {
            next_board = symmetry(flip(do_move(initial_board, (1 << pos))));
            init_hash = hash(next_board);
            if (valid_move(next_board) == 0) {
                flip_board = symmetry(flip(next_board));
                flip_hash = hash(flip_board);
                if (valid_move(flip_board)) {
                    write_bits_w_hash(flip_board.opponent | flip_board.player, FOUND, flip_hash);
                } else {
                    if (flip_hash < init_hash) {
                        write_bits_w_hash(flip_board.opponent | flip_board.player, FOUND, flip_hash);
                    } else {
                        write_bits_w_hash(next_board.opponent | next_board.player, FOUND, init_hash);
                    }
                }
            } else {
                write_bits_w_hash(next_board.opponent | next_board.player, FOUND, init_hash);
            }
            
        }
        valid = valid >> 1;
        pos++;
    }
}