// board = initial_board
// board_shape = shape(initial_board)

// board_shapes = dictionary of board shapes
// board_shapes[board_shape] = [some array, with enough size to store all possible positions]
// board_shapes[board_shape][canonical(board)] = Found

// linked_list_of_boards = empty linked list
// linked_list_of_boards.add({board_shape})

// max_positions = dictionary of largest position
// max_positions[board_shape] = canonical(board)

// if linked list not empty:
//     current_board = get the first item in the linked list
//     for each found_position in the current board:
//         for each move for found_position:
//             new_board = canonical(do_move(move, found_position))
//             new_board_shape = shape(new_board)
//             if new_board_shape not in board_shapes.keys():
// 	        board_shapes[new_board_shape] = [some array, with enough size]
//                 linked_list_of_boards.add({board_shape})
// main.cpp
#include <iostream>
#include <cstdint>
#include <iomanip>
#include "othello6.hpp"   // your functions live here

using bitops6::Bitboard;

// Print helpers
static void print_board(Bitboard b) {
    for (int r = 5; r >= 0; --r) {
        for (int c = 0; c < 6; ++c) {
            int idx = r * 6 + c;
            std::cout << (((b >> idx) & 1ULL) ? "1 " : ". ");
        }
        std::cout << "\n";
    }
    std::cout << "----\n";
}

static void print_pos(const games::Encoding& e) {
    std::cout << "Player:\n";   print_board(e.player);
    std::cout << "Opponent:\n"; print_board(e.opponent);
}

int main() {
    auto bit = [](int r, int c){ return 1ULL << (r * 6 + c); };

    // Example: standard 6×6 start (center 4 occupied)
    games::Encoding pos{};
    pos.player   = bit(2,3) | bit(3,2);
    pos.opponent = bit(2,2) | bit(3,3);

    std::cout << "=== Input position ===\n";
    print_pos(pos);

    // 1. legal_moves
    Bitboard moves = othello6::legal_moves(pos);
    std::cout << "Legal moves (hex) = 0x" << std::hex << moves << std::dec << "\n";
    std::cout << "Legal moves grid:\n";
    print_board(moves);

    // 2. Apply each legal move
    Bitboard m = moves;
    while (m) {
        Bitboard mv = m & -m;  // extract lowest set bit
        m ^= mv;

        int idx = __builtin_ctzll(mv);
        int r = idx / 6, c = idx % 6;
        std::cout << "\n-- Move at (r=" << r << ", c=" << c << ") --\n";

        games::Encoding after = othello6::do_move(pos, mv);
        print_pos(after);

        // 3. flip
        std::cout << "After flip:\n";
        print_pos(othello6::flip(after));

        // 4. canonical
        std::cout << "Canonical representative:\n";
        print_pos(othello6::canonical(after));

        // 5. shape
        Bitboard shp = othello6::shape(after);
        std::cout << "Shape (hex) = 0x" << std::hex << shp << std::dec << "\n";
        print_board(shp);

        // 6. hash
        std::uint64_t h = othello6::hash(after);
        std::cout << "Hash = 0x" << std::hex << h << std::dec << "\n";

        // 7. unhash
        games::Encoding round = othello6::unhash(shp, h);
        std::cout << "Unhashed back:\n";
        print_pos(round);

        // Quick consistency check
        if ((after.player == round.player) && (after.opponent == round.opponent)) {
            std::cout << "Round-trip hash/unhash OK\n";
        } else {
            std::cout << "Mismatch in hash/unhash round-trip!\n";
        }
    }

    return 0;
}
