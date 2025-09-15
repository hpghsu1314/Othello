#pragma once
#include <cstdint>

namespace bitops4 {

    using Bitboard = std::uint64_t;
    inline constexpr int N = 4;
    inline constexpr int CELLS = N * N;
    inline constexpr Bitboard FULL = (CELLS == 64) ? ~0ULL : ((1ULL << CELLS) - 1ULL); // 0xFFFF

    // --- Row masks (4 bits each) ---
    static constexpr Bitboard R0 = 0x000FULL;   // row0 (bits 0..3)
    static constexpr Bitboard R1 = R0 << 4;     // row1 (4..7)
    static constexpr Bitboard R2 = R0 << 8;     // row2 (8..11)
    static constexpr Bitboard R3 = R0 << 12;    // row3 (12..15)

    // --- Column masks (4 bits tall) ---
    static constexpr Bitboard C0 = 0x00001111ULL; // col0: bits 0,4,8,12
    static constexpr Bitboard C1 = 0x00002222ULL; // col1: bits 1,5,9,13
    static constexpr Bitboard C2 = 0x00004444ULL; // col2: bits 2,6,10,14
    static constexpr Bitboard C3 = 0x00008888ULL; // col3: bits 3,7,11,15

    static constexpr Bitboard A_FILE = C0;
    static constexpr Bitboard D_FILE = C3;
    static constexpr Bitboard NOT_A  = FULL ^ A_FILE;
    static constexpr Bitboard NOT_D  = FULL ^ D_FILE;

    // --- Vertical flip (left <-> right) ---
    inline Bitboard vertical4(Bitboard b) {
        Bitboard res = 0;
        res |= ((b & C0) << 3);  // col0 -> col3
        res |= ((b & C1) << 1);  // col1 -> col2
        res |= ((b & C2) >> 1);  // col2 -> col1
        res |= ((b & C3) >> 3);  // col3 -> col0
        return res;
    }

    // --- Horizontal flip (top <-> bottom) ---
    inline Bitboard horizontal4(Bitboard b) {
        Bitboard res = 0;
        res |= ((b & R0) << 12); // row0 -> row3
        res |= ((b & R1) << 4 ); // row1 -> row2
        res |= ((b & R2) >> 4 ); // row2 -> row1
        res |= ((b & R3) >> 12); // row3 -> row0
        return res;
    }

    // --- Transpose across main diagonal (4x4) ---
    inline Bitboard transpose4(Bitboard b) {
        Bitboard res = 0;

        // Diagonal stays
        static constexpr Bitboard DIAG = (1ULL<<0) | (1ULL<<5) | (1ULL<<10) | (1ULL<<15); // 0x8421
        res |= (b & DIAG);

        // Offset 1: (0,1),(1,2),(2,3) <-> (1,0),(2,1),(3,2)   shift = (N-1)*1 = 3
        static constexpr Bitboard R1 = (1ULL<<1) | (1ULL<<6) | (1ULL<<11);   // 0x842
        static constexpr Bitboard L1 = (R1 << 3);                            // 0x4210
        res |= ((b & R1) << 3);
        res |= ((b & L1) >> 3);

        // Offset 2: (0,2),(1,3) <-> (2,0),(3,1)               shift = 6
        static constexpr Bitboard R2 = (1ULL<<2) | (1ULL<<7);                // 0x84
        static constexpr Bitboard L2 = (R2 << 6);                            // 0x2100
        res |= ((b & R2) << 6);
        res |= ((b & L2) >> 6);

        // Offset 3: (0,3) <-> (3,0)                           shift = 9
        static constexpr Bitboard R3 = (1ULL<<3);                            // 0x8
        static constexpr Bitboard L3 = (R3 << 9);                            // 0x1000
        res |= ((b & R3) << 9);
        res |= ((b & L3) >> 9);

        return res;
    }

} // namespace bitops4
