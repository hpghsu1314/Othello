#ifndef OTHELLO_4_H
#define OTHELLO_4_H

#include "game.h"

#define BOARD_SIZE 4 // Valid for 4, 6, 8

// Will make a function to make this automatically
#define NOT_A_FILE 0b0111011101110111 // Remember to change this
#define NOT_F_FILE 0b1110111011101110 // Remember to change this
#define FULL_BOARD 0x000000000000FFFF // Remember to change this


Bitboard transpose4(Bitboard board);
Bitboard vertical4(Bitboard board);
Bitboard horizontal4(Bitboard board);

#endif