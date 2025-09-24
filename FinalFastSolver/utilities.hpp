#pragma once

#include <cstdint>

static inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

static inline int owner_of_shape(std::uint64_t shape, int W) {
    return (int)(splitmix64(shape) % (std::uint64_t)W);
}
