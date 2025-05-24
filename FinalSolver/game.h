#ifndef GAME_H
#define GAME_H

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "file_lst.h"
#include "num_lst.h"

#ifndef FOLDER
#define FOLDER "DefaultFolder"
#endif

#define WIN 0b11000000
#define DRAW 0b10000000
#define LOSE 0b01000000
#define SKIP 0b00000000

#define FOUND 0b11111111

struct Encoding {
    uint64_t player;
    uint64_t opponent;
};

typedef uint64_t Bitboard;

// gcc -o solver solver.c othello4.c file_lst.c num_lst.c game.c

// Provides valid moves given an Encoding
// Valid moves will be denoted 1 in the return value
Bitboard valid_move(struct Encoding board);

// Given an Encoding and a move, return the result
struct Encoding do_move(struct Encoding board, uint64_t move);

// Returns a deterministic representation of a given Encoding
// i.e. Symmetric Encodings should return the same Encoding
struct Encoding symmetry(struct Encoding board);

// Returns the hash of a given Encoding
uint64_t hash(struct Encoding board);

// Given a board_shape and hash_value, return the corresponding Encoding
struct Encoding unhash(uint64_t board_shape, uint64_t hash_value);

// Given an Encoding, return its primitive value
uint8_t primitive(struct Encoding board);

// Write some information to the corresponding file given hash, the value to write, and board_shape
void write_bits_w_hash(uint64_t board_shape, uint8_t value, uint64_t hash_value);

// Do all possible moves and write corresponding information
// Only used for finding positions
// Not used for finding primitive values
void do_moves_and_write(struct Encoding initial_board);

// Flip player and opponent in Encoding
struct Encoding flip(struct Encoding board);

// Given a Bitboard, return the number of 1s
int count_bits(Bitboard board);

// Read the corresponding file at the position of the hash of the board
uint8_t file_read_board(struct Encoding board);

#endif