#pragma once
#include <cstdint>

namespace games {

    // Shared byte encodings
    inline constexpr std::uint8_t WIN  = 0b1100'0000;
    inline constexpr std::uint8_t DRAW = 0b1000'0000;
    inline constexpr std::uint8_t LOSE = 0b0100'0000;
    inline constexpr std::uint8_t SKIP = 0b0000'0000;
    inline constexpr std::uint8_t FOUND = 0b1111'1111;

    using Bitboard = std::uint64_t;

    // Minimal “position” for two-player, bitboard-style games.
    struct Encoding {
        Bitboard player;   // bits for side-to-move
        Bitboard opponent; // bits for opponent
    };

    // Abstract interface your solver talks to.
    // Keep it small to avoid virtual overhead in hot loops—inline/cache results on the caller side when possible.
    class Game {
    public:
        virtual ~Game() = default;

        // Bitmask of legal moves (1-bit per cell).
        virtual std::uint64_t legal_moves(const Encoding& s) const noexcept = 0;

        // Apply a single-bit move and return new position (still two bitboards).
        virtual Encoding do_move(const Encoding& s, std::uint64_t move) const noexcept = 0;

        // Swap sides (side-to-move flip).
        virtual Encoding flip(const Encoding& s) const noexcept = 0;

        // Canonicalize under the game’s symmetries (for transposition collapsing).
        virtual Encoding canonical(const Encoding& s) const noexcept = 0;

        // Shape = occupancy mask (used for your sharding + page-table layout).
        virtual std::uint64_t shape(const Encoding& s) const noexcept = 0;

        // Compact hash by occupancy order (your “pack only occupied squares” scheme).
        virtual std::uint64_t hash(const Encoding& s) const noexcept = 0;

        // Inverse of hash_compact for a given shape.
        virtual Encoding unhash(std::uint64_t shape, std::uint64_t h) const noexcept = 0;

        // Exact terminal evaluation (+ remoteness in low bits, if you encode it that way).
        // For non-terminals you can either recurse here or have caller orchestrate search.
        virtual std::uint8_t primitive(const Encoding& s) const = 0;
    };

}
