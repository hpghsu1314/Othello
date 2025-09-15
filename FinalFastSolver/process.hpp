#pragma once
#include <cstdint>

namespace process {
    static inline uint64_t hash_shape(uint64_t shapeID) {
        uint64_t z = shapeID + 0x9e3779b97f4a7c15ull;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        return z ^ (z >> 31);
    }
    struct Owner {
        int world, rank;
        int operator()(uint64_t shapeID) const { return int(hash_shape(shapeID) % uint64_t(world)); }
        bool is_owner(uint64_t shapeID) const { return (*this)(shapeID) == rank; }
    };

}