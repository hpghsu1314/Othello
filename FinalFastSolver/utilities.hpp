#pragma once
#include <cstdint>

#if defined(_MSC_VER) && defined(_M_X64)
  #include <intrin.h>
#endif

// MurmurHash3 64-bit finalizer: good avalanche/mixing for 64->64.
static inline std::uint64_t mix64(std::uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

// Lemire's fast unbiased reduction: maps x uniformly into [0, n)
static inline std::uint64_t fast_reduce(std::uint64_t x, std::uint64_t n) {
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned __int64 hi;
    (void)_umul128(x, n, &hi);      // high 64 bits of x*n
    return static_cast<std::uint64_t>(hi);
#else
    return static_cast<std::uint64_t>(((__uint128_t)x * (__uint128_t)n) >> 64);
#endif
}

// Deterministic owner mapping with a baked-in salt to decorrelate low-entropy shapes.
static inline int owner_of_shape(std::uint64_t shape, int W) {
    if (W <= 0) return 0;
    constexpr std::uint64_t kSalt = 0x9E3779B97F4A7C15ULL; // golden ratio-based
    const std::uint64_t h = mix64(shape ^ kSalt);
    return static_cast<int>(fast_reduce(h, static_cast<std::uint64_t>(W)));
}
