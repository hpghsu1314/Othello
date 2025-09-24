#pragma once
#include <cstdint>
#include "bitops6.hpp"
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
  #include <immintrin.h>  // for _pext_u64 / _pdep_u64 (BMI2)
#endif
#if defined(_MSC_VER)
  #include <intrin.h>     // __popcnt64 on MSVC
#endif
#if __cpp_lib_bitops >= 201907L
  #include <bit>            // std::popcount
#endif

#if defined(_MSC_VER)
  #define OTH_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
  #define OTH_FORCE_INLINE inline __attribute__((always_inline))
#else
  #define OTH_FORCE_INLINE inline
#endif

namespace game {
    using bitops6::Bitboard;
    using bitops6::FULL;
    using bitops6::NOT_A;
    using bitops6::NOT_F;
    using bitops6::vertical6;
    using bitops6::horizontal6;
    using bitops6::transpose6;
    using bitops6::N;

    // Shared byte encodings
    inline constexpr std::uint8_t WIN  = 0b1100'0000;
    inline constexpr std::uint8_t DRAW = 0b1000'0000;
    inline constexpr std::uint8_t LOSE = 0b0100'0000;
    inline constexpr std::uint8_t SKIP = 0b0000'0000;
    inline constexpr std::uint8_t FOUND = 0b1111'1111;
    inline constexpr std::uint8_t MAX_TIER = 36;

    // Minimal “position” for two-player, bitboard-style games.
    struct Encoding {
        Bitboard player;   // bits for side-to-move
        Bitboard opponent; // bits for opponent
    };

    OTH_FORCE_INLINE Bitboard ks_dir_l(Bitboard me, Bitboard opp, Bitboard empty, int s, Bitboard mask) {
        Bitboard t = opp & mask & (me << s);
        t |= opp & mask & (t << s);
        t |= opp & mask & (t << s);
        t |= opp & mask & (t << s);
        t |= opp & mask & (t << s);   // 5 cascades for 6x6
        return (t << s) & mask & empty;
    }
    OTH_FORCE_INLINE Bitboard ks_dir_r(Bitboard me, Bitboard opp, Bitboard empty, int s, Bitboard mask) {
        Bitboard t = opp & mask & (me >> s);
        t |= opp & mask & (t >> s);
        t |= opp & mask & (t >> s);
        t |= opp & mask & (t >> s);
        t |= opp & mask & (t >> s);
        return (t >> s) & mask & empty;
    }

    OTH_FORCE_INLINE Bitboard legal_moves(const Encoding& pos) {
        const Bitboard me = pos.player & FULL;
        const Bitboard opp = pos.opponent & FULL;
        const Bitboard empty = FULL & ~(me | opp);

        Bitboard m = 0;
        m |= ks_dir_l(me, opp, empty, 1, NOT_A); // east
        m |= ks_dir_r(me, opp, empty, 1, NOT_F); // west
        m |= ks_dir_l(me, opp, empty, 6, FULL); // north
        m |= ks_dir_r(me, opp, empty, 6, FULL); // south
        m |= ks_dir_l(me, opp, empty, 7, NOT_A); // NE
        m |= ks_dir_r(me, opp, empty, 7, NOT_F); // SW
        m |= ks_dir_l(me, opp, empty, 5, NOT_F); // NW
        m |= ks_dir_r(me, opp, empty, 5, NOT_A); // SE
        return m & FULL;
    }

    // Collect flips in the "left-shift" directions (E, N, NE, NW)
    OTH_FORCE_INLINE Bitboard sweep_l(Bitboard mv, Bitboard me, Bitboard opp, int s, Bitboard mask) {
        Bitboard flips = 0;
        Bitboard cur   = (mv << s) & mask;
        // Walk through opponent stones
        while (cur && (cur & opp)) {
            flips |= cur;
            cur = (cur << s) & mask;
        }
        // If we ended on our own stone, the path is valid; otherwise discard
        return (cur & me) ? flips : 0;
    }

    // Collect flips in the "right-shift" directions (W, S, SW, SE)
    OTH_FORCE_INLINE Bitboard sweep_r(Bitboard mv, Bitboard me, Bitboard opp, int s, Bitboard mask) {
        Bitboard flips = 0;
        Bitboard cur   = (mv >> s) & mask;
        while (cur && (cur & opp)) {
            flips |= cur;
            cur = (cur >> s) & mask;
        }
        return (cur & me) ? flips : 0;
    }

    // Apply a single-bit move 'mv' for side-to-move in 'pos'
    inline Encoding do_move(const Encoding& pos, Bitboard mv) {
        Bitboard me  = pos.player & FULL;
        Bitboard opp = pos.opponent & FULL;

        const Bitboard flips =
            sweep_l(mv, me, opp, 1, NOT_A) |  // East
            sweep_r(mv, me, opp, 1, NOT_F) |  // West
            sweep_l(mv, me, opp, N, FULL ) |  // North
            sweep_r(mv, me, opp, N, FULL ) |  // South
            sweep_l(mv, me, opp, N+1, NOT_A) |  // NE
            sweep_r(mv, me, opp, N+1, NOT_F) |  // SW
            sweep_l(mv, me, opp, N-1, NOT_F) |  // NW
            sweep_r(mv, me, opp, N-1, NOT_A);   // SE

        me  ^= flips;
        opp ^= flips;
        me  |= mv;

        return { me & FULL, opp & FULL };
    }

    OTH_FORCE_INLINE Encoding flip(const Encoding& pos) {
        return { pos.opponent, pos.player };
    }

    inline Encoding canonical(const Encoding& pos) noexcept {
        const Bitboard me  = pos.player   & FULL;
        const Bitboard opp = pos.opponent & FULL;
        const Bitboard occ = me | opp;

        // Occupancy transforms (each exactly once)
        const Bitboard occ0 = occ;                   // id
        const Bitboard occ1 = vertical6(occ);        // V
        const Bitboard occ2 = horizontal6(occ);      // H
        const Bitboard occ3 = transpose6(occ);       // T
        const Bitboard occ5 = vertical6(occ3);       // R90 = V∘T
        const Bitboard occ7 = horizontal6(occ3);     // R270 = H∘T
        const Bitboard occ6 = horizontal6(occ1);     // R180 = H∘V
        const Bitboard occ4 = horizontal6(occ5);     // ADIAG= H∘R90

        // Find minimal occupancy; track ties in a bitmask
        Bitboard best_occ = occ0;
        unsigned tie = 1u << 0;
        auto upd = [&](Bitboard x, unsigned b){
            const bool lt = x < best_occ;
            const bool eq = x == best_occ;
            best_occ = lt ? x : best_occ;
            tie = lt ? b : (eq ? (tie | b) : tie);
        };
        upd(occ1,1u<<1); upd(occ2,1u<<2); upd(occ3,1u<<3);
        upd(occ4,1u<<4); upd(occ5,1u<<5); upd(occ6,1u<<6); upd(occ7,1u<<7);

        // Lazily compute 'me' only for tied transforms; avoid repeats where shared
        Bitboard best_me = ~0ULL;

        // id
        if (tie & (1u<<0)) {
            const Bitboard m0 = me;
            if (m0 < best_me) best_me = m0;
        }

        // V
        Bitboard me_v = 0; bool have_v = false;
        if (tie & (1u<<1)) {
            me_v = vertical6(me); have_v = true;
            if (me_v < best_me) best_me = me_v;
        }

        // H
        if (tie & (1u<<2)) {
            const Bitboard m2 = horizontal6(me);
            if (m2 < best_me) best_me = m2;
        }

        // T
        Bitboard me_t = 0; bool have_t = false;
        if (tie & (1u<<3)) {
            me_t = transpose6(me); have_t = true;
            if (me_t < best_me) best_me = me_t;
        }

        // ADIAG = H∘R90 = H∘V∘T
        if (tie & (1u<<4)) {
            if (!have_t) { me_t = transpose6(me); have_t = true; }
            const Bitboard m5tmp = vertical6(me_t);      // R90
            const Bitboard m4 = horizontal6(m5tmp);   // ADIAG
            if (m4 < best_me) best_me = m4;
        }

        // R90 = V∘T
        if (tie & (1u<<5)) {
            if (!have_t) { me_t = transpose6(me); have_t = true; }
            const Bitboard m5 = vertical6(me_t);
            if (m5 < best_me) best_me = m5;
        }

        // R180 = H∘V
        if (tie & (1u<<6)) {
            if (!have_v) { me_v = vertical6(me); have_v = true; }
            const Bitboard m6 = horizontal6(me_v);
            if (m6 < best_me) best_me = m6;
        }

        // R270 = H∘T
        if (tie & (1u<<7)) {
            if (!have_t) { me_t = transpose6(me); have_t = true; }
            const Bitboard m7 = horizontal6(me_t);
            if (m7 < best_me) best_me = m7;
        }

        return { best_me, best_occ ^ best_me };
    }

    // 1) Occupancy "shape"
    OTH_FORCE_INLINE Bitboard shape(const Encoding& s) noexcept {
        return (s.player | s.opponent) & FULL;
    }

    // 2) Pack player's bits over the shape into a compact hash (low bits)
    //    BMI2 _pext_u64 path (fast), with a portable fallback.
    inline std::uint64_t hash(const Encoding& s) noexcept {
        const Bitboard sh = shape(s);
    #if defined(__BMI2__)
        return _pext_u64(s.player, sh);
    #else
        // Portable fallback: iterate set bits of sh (in index order) and pack.
        std::uint64_t h = 0;
        std::uint64_t pos = 0;
        Bitboard m = sh;
        while (m) {
            Bitboard b = m & -m;
            if (s.player & b) h |= (1ULL << pos);
            ++pos;
            m ^= b;                             // clear the bit
        }
        return h;
    #endif
    }

    // 3) Unpack compact hash back into full bitboards (player/opponent) given a shape.
    //    BMI2 _pdep_u64 path (fast), with a portable fallback.
    inline Encoding unhash(Bitboard sh, std::uint64_t h) noexcept {
    #if defined(__BMI2__)
        const Bitboard me = _pdep_u64(h, sh);
        const Bitboard occ = sh;
        return { me & FULL, (occ ^ me) & FULL };
    #else
        Bitboard me = 0;
        std::uint64_t m = sh, src = h;
        while (m) {
            std::uint64_t b = m & -m;
            if (src & 1ull) me |= b;
            src >>= 1;
            m ^= b;
        }
        const Bitboard occ = sh;
        return { me & FULL, (occ ^ me) & FULL };
    #endif
    }

    inline std::uint8_t primitive(const Encoding& pos) noexcept {
        const Bitboard me = pos.player;
        const Bitboard opp = pos.opponent;

        if ((legal_moves(pos) | legal_moves(flip(pos))) != 0) return 0;

        int my_count, opp_count;
    #if defined(__GNUC__) || defined(__clang__)
        my_count = __builtin_popcountll(me);
        opp_count = __builtin_popcountll(opp);
    #elif defined(_MSC_VER)
        my_count = static_cast<int>(__popcnt64(me));
        opp_count = static_cast<int>(__popcnt64(opp));
    #else
        // Portable & fast enough
        auto pc = [](Bitboard x){ int c=0; while (x){ x&=x-1; ++c; } return c; };
        my_count = pc(me);
        opp_count = pc(opp);
    #endif

        if (my_count > opp_count) return WIN;
        if (my_count < opp_count) return LOSE;
        return DRAW;
    }

    OTH_FORCE_INLINE std::uint8_t tier_of(Bitboard shapeMask) noexcept {
    #if __cpp_lib_bitops >= 201907L
        return static_cast<std::uint8_t>(std::popcount(shapeMask));
    #elif defined(_MSC_VER)
        return static_cast<std::uint8_t>(__popcnt64(shapeMask));
    #elif defined(__GNUC__) || defined(__clang__)
        return static_cast<std::uint8_t>(__builtin_popcountll(shapeMask));
    #else
        std::uint64_t x = shapeMask; int c = 0; while (x) { x &= (x - 1); ++c; }
        return static_cast<std::uint8_t>(c);
    #endif
    }

    OTH_FORCE_INLINE Encoding starting_position() {
        return game::Encoding{ 0b000000000000000100001000000000000000, 0b000000000000001000000100000000000000};
    }

    struct Game {
        using Encoding = game::Encoding;

        std::uint64_t legal_moves(const Encoding& s) const noexcept {
            return game::legal_moves(s);
        }
        Encoding do_move(const Encoding& s, std::uint64_t mv) const noexcept {
            return game::do_move(s, mv);
        }
        Encoding flip(const Encoding& s) const noexcept {
            return game::flip(s);
        }
        Encoding canonical(const Encoding& s) const noexcept {
            return game::canonical(s);
        }
        std::uint64_t shape(const Encoding& s) const noexcept {
            return game::shape(s);
        }
        std::uint64_t hash(const Encoding& s) const noexcept {
            return game::hash(s);
        }
        Encoding unhash(std::uint64_t sh, std::uint64_t h) const noexcept {
            return game::unhash(sh, h);
        }
        std::uint8_t primitive(const Encoding& s) const noexcept {
            return game::primitive(s);
        }
        uint8_t tier_of(Bitboard shapemask) const noexcept {
            return game::tier_of(shapemask);
        }
        Encoding starting_position() const noexcept {
            return game::starting_position();
        }
    };
}
#undef OTH_FORCE_INLINE
