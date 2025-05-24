#include "game.h"

struct Encoding flip(struct Encoding board) {
    struct Encoding ret_val;
    ret_val.opponent = board.player;
    ret_val.player = board.opponent;
    return ret_val;
}

int count_bits(Bitboard board) {
    int ret_val = 0;
    while (board > 0) {
        if (board % 2 == 1) {
            ret_val++;
        }
        board = board >> 1;
    }
    return ret_val;
}