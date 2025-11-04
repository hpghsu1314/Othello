#pragma once          // prevents multiple inclusion
#include <cstdint>     // for std::uint64_t, std::uint32_t, std::uint8_t

// Ensure 1-byte packing
#pragma pack(push, 1)
struct HitMsg {
    std::uint64_t shape;
    std::uint32_t page;
    std::uint32_t bit;
};
struct TierMsg : HitMsg {
    std::uint8_t tier;
};
#pragma pack(pop)