#pragma once
#include <cstdint>
#include <algorithm> 
#include <ranges>     

namespace bitops6 {

    using Bitboard = std::uint64_t;
    inline constexpr int N = 6;
    inline constexpr int CELLS = N * N;
    inline constexpr Bitboard FULL = (CELLS == 64) ? ~0ULL : ((1ULL << CELLS) - 1ULL);

    // --- Row masks (6 bits each) ---
    static constexpr Bitboard R0 = 0x00000003FULL;        // row0 (bits 0..5)
    static constexpr Bitboard R1 = R0 << 6;                  // row1 (6..11)
    static constexpr Bitboard R2 = R0 << 12;                 // row2 (12..17)
    static constexpr Bitboard R3 = R0 << 18;                 // row3 (18..23)
    static constexpr Bitboard R4 = R0 << 24;                 // row4 (24..29)
    static constexpr Bitboard R5 = R0 << 30;                 // row5 (30..35)

    // --- Column masks (6 bits tall) ---
    static constexpr Bitboard C0 = 0x041041041ULL;        // col0
    static constexpr Bitboard C1 = 0x082082082ULL;        // col1
    static constexpr Bitboard C2 = 0x104104104ULL;        // col2
    static constexpr Bitboard C3 = 0x208208208ULL;        // col3
    static constexpr Bitboard C4 = 0x410410410ULL;        // col4
    static constexpr Bitboard C5 = 0x820820820ULL;        // col5

    static constexpr Bitboard A_FILE = C0;
    static constexpr Bitboard F_FILE = C5;
    static constexpr Bitboard NOT_A  = FULL ^ A_FILE;
    static constexpr Bitboard NOT_F  = FULL ^ F_FILE;

    // --- Vertical flip (left <-> right) ---
    inline Bitboard vertical6(Bitboard b) {
        Bitboard res = 0;
        res |= ((b & C0) << 5);  // col0 -> col5
        res |= ((b & C1) << 3);  // col1 -> col4
        res |= ((b & C2) << 1);  // col2 -> col3
        res |= ((b & C3) >> 1);  // col3 -> col2
        res |= ((b & C4) >> 3);  // col4 -> col1
        res |= ((b & C5) >> 5);  // col5 -> col0
        return res;
    }

    // --- Horizontal flip (top <-> bottom) ---
    inline Bitboard horizontal6(Bitboard b) {
        Bitboard res = 0;
        res |= ((b & R0) << 30); // row0 -> row5
        res |= ((b & R1) << 18); // row1 -> row4
        res |= ((b & R2) << 6);  // row2 -> row3
        res |= ((b & R3) >> 6);  // row3 -> row2
        res |= ((b & R4) >> 18); // row4 -> row1
        res |= ((b & R5) >> 30); // row5 -> row0
        return res;
    }

    // --- Transpose across main diagonal ---
    // Swap symmetric pairs across diagonal
    inline Bitboard transpose6(Bitboard b) {
        Bitboard res = 0;

        static constexpr Bitboard DIAG =
            (1ULL<<0) | (1ULL<<7) | (1ULL<<14) | (1ULL<<21) | (1ULL<<28) | (1ULL<<35);
        res |= (b & DIAG);


        static constexpr Bitboard R1 =
            (1ULL<<1) | (1ULL<<8) | (1ULL<<15) | (1ULL<<22) | (1ULL<<29);
        static constexpr Bitboard L1 = (R1 << (5*1));
        res |= ((b & R1) << (5*1));
        res |= ((b & L1) >> (5*1));


        static constexpr Bitboard R2 =
            (1ULL<<2) | (1ULL<<9) | (1ULL<<16) | (1ULL<<23);
        static constexpr Bitboard L2 = (R2 << (5*2));
        res |= ((b & R2) << (5*2));
        res |= ((b & L2) >> (5*2));


        static constexpr Bitboard R3 =
            (1ULL<<3) | (1ULL<<10) | (1ULL<<17);
        static constexpr Bitboard L3 = (R3 << (5*3));
        res |= ((b & R3) << (5*3));
        res |= ((b & L3) >> (5*3));


        static constexpr Bitboard R4 =
            (1ULL<<4) | (1ULL<<11);
        static constexpr Bitboard L4 = (R4 << (5*4));
        res |= ((b & R4) << (5*4));
        res |= ((b & L4) >> (5*4));


        static constexpr Bitboard R5 = (1ULL<<5);
        static constexpr Bitboard L5 = (R5 << (5*5));
        res |= ((b & R5) << (5*5));
        res |= ((b & L5) >> (5*5));

        return res;
    }
}
