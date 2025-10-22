#pragma once
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <array>
#include <vector>
#include <filesystem>
#include <functional>
#include "othello6.hpp"
#include "bitops6.hpp"

namespace pageops {

namespace fs = std::filesystem;

static constexpr std::size_t PAGE_BITS   = 12;
static constexpr std::size_t PAGE_BYTES  = (1u << (PAGE_BITS - 3));              // 2^(PAGE_BITS)/8
static constexpr std::size_t PAGE_QWORDS = PAGE_BYTES / sizeof(std::uint64_t);
static constexpr std::uint8_t REMOTENESS_BITS = 5;

struct Key {
    std::uint64_t shape;
    std::uint32_t page;
    inline bool operator==(const Key& o) const noexcept { return shape==o.shape && page==o.page; }
};
struct KeyHash {
    inline std::size_t operator()(const Key& k) const noexcept {
        std::uint64_t x = k.shape ^ (std::uint64_t(k.page) * 0x9e3779b97f4a7c15ULL);
        x ^= (x>>33); x *= 0xff51afd7ed558ccdULL;
        x ^= (x>>33); x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= (x>>33);
        return (std::size_t)x;
    }
};

#pragma pack(push,1)
struct PageIdxRec { std::uint64_t shape; std::uint32_t page; std::uint64_t off; };
#pragma pack(pop)
static_assert(sizeof(PageIdxRec)==20, "PageIdxRec must be 20 bytes");

struct PageBuf {
    alignas(64) std::array<std::uint64_t, PAGE_QWORDS> w{};
    bool dirty{false};
    PageBuf() = default;
    PageBuf(const PageBuf&)            = delete;
    PageBuf& operator=(const PageBuf&) = delete;
};

using MetaMap  = std::unordered_map<Key, std::uint64_t, KeyHash>;
using CacheMap = std::unordered_map<Key, PageBuf,   KeyHash>;

struct Context {
    // Hot file handles (kept open)
    std::FILE* data = nullptr;  // pages.dat
    std::FILE* idx  = nullptr;  // pages.idx
    std::FILE* solved = nullptr; // pages.solved
    std::size_t file_size = 0;  // current end of pages.dat (bytes)

    // Paths
    std::string data_path;
    std::string idx_path;
    std::string solved_path;

    // In-RAM maps
    MetaMap  meta_map;          // (shape,page) -> offset
    CacheMap cache;             // (shape,page) -> page buffer (in-place, no moves/copies)

    // Cache size knob (pages)
    std::size_t cache_cap_pages = 4096u * 16u;

    bool inited{false};
};

// One shared zero page for append
inline const std::uint8_t ZERO_PAGE[PAGE_BYTES] = {0};

// ---------- init / shutdown ----------

inline void pages_init(Context& ctx, int owner_num, unsigned tier, const std::string& root_dir) {
    if (ctx.inited) return;

    char dir[512], dpath[512], ipath[512], spath[512];
    std::snprintf(dir,   sizeof(dir),   "%s/owner_%02d/tier_%02u", root_dir.c_str(), owner_num, tier);
    std::snprintf(dpath, sizeof(dpath), "%s/pages.dat", dir);
    std::snprintf(ipath, sizeof(ipath), "%s/pages.idx", dir);
    std::snprintf(spath, sizeof(spath), "%s/pages.solved", dir);
    fs::create_directories(dir);

    ctx.data_path = dpath;
    ctx.idx_path  = ipath;
    ctx.solved_path = spath;

    ctx.data = std::fopen(ctx.data_path.c_str(), "rb+");
    if (!ctx.data) ctx.data = std::fopen(ctx.data_path.c_str(), "wb+");
    if (!ctx.data) { std::perror("fopen pages.dat"); std::exit(1); }
    std::setvbuf(ctx.data, nullptr, _IOFBF, 1<<20); // 1 MiB

    ctx.idx = std::fopen(ctx.idx_path.c_str(), "ab+");
    if (!ctx.idx) { std::perror("fopen pages.idx"); std::exit(1); }
    std::setvbuf(ctx.idx, nullptr, _IOFBF, 1<<20);

    ctx.solved = std::fopen(ctx.solved_path.c_str(), "rb+");
    if (!ctx.solved) ctx.solved = std::fopen(ctx.solved_path.c_str(), "wb+");
    if (!ctx.solved) { std::perror("fopen pages.solved"); std::exit(1); }
    std::setvbuf(ctx.solved, nullptr, _IOFBF, 1<<20);

    std::fseek(ctx.data, 0, SEEK_END);
    ctx.file_size = (std::size_t)std::ftell(ctx.data);

    ctx.inited = true;
}

inline void pages_flush(Context& ctx) {
    if (ctx.data) std::fflush(ctx.data);
    if (ctx.idx)  std::fflush(ctx.idx);
    if (ctx.solved) std::fflush(ctx.solved);
}

inline void pages_close(Context& ctx) {
    pages_flush(ctx);
    if (ctx.data) { std::fclose(ctx.data); ctx.data=nullptr; }
    if (ctx.idx)  { std::fclose(ctx.idx);  ctx.idx=nullptr;  }
    if (ctx.solved) { std::fclose(ctx.solved); ctx.solved=nullptr; }
    ctx.inited = false;
}

// ---------- flush (no threading) ----------

inline void flush_page(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    const Key k{shape,page_id};
    auto it = ctx.cache.find(k);
    if (it == ctx.cache.end()) return;
    PageBuf& p = it->second;

    if (p.dirty) {
        const auto mit = ctx.meta_map.find(k);
        if (mit == ctx.meta_map.end()) { std::fprintf(stderr, "flush_page: missing meta\n"); std::exit(1); }
        const std::uint64_t off = mit->second;

        std::fseek(ctx.data, (long)off, SEEK_SET);
        if (std::fwrite(p.w.data(), 1, PAGE_BYTES, ctx.data) != PAGE_BYTES) {
            std::perror("fwrite pages.dat"); std::exit(1);
        }
        p.dirty = false;
    }
}

inline void flush_all(Context& ctx) {
    for (auto &kv : ctx.cache) {
        const Key& k = kv.first;
        flush_page(ctx, k.shape, k.page);
    }
    if (ctx.data) std::fflush(ctx.data);
    if (ctx.idx)  std::fflush(ctx.idx);
}

inline std::uint64_t offset_of(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    auto it = ctx.meta_map.find(Key{shape,page_id});
    return (it==ctx.meta_map.end()) ? ~0ULL : it->second;
}

// ---------- ensure page exists (single-threaded), return offset ----------

inline std::uint64_t ensure_page_exists_fast(Context& ctx,
                                             std::uint64_t shape,
                                             std::uint32_t page_id)
{
    const Key k{shape,page_id};

    if (auto it = ctx.meta_map.find(k); it != ctx.meta_map.end())
        return it->second;

    // Reserve unique offset and publish mapping
    const std::uint64_t off = ctx.file_size;
    ctx.file_size += PAGE_BYTES;
    ctx.meta_map.try_emplace(k, off);

    // Write zero page at the reserved offset
    std::fseek(ctx.data, (long)off, SEEK_SET);
    if (std::fwrite(ZERO_PAGE, 1, PAGE_BYTES, ctx.data) != PAGE_BYTES) {
        std::perror("fwrite pages.dat"); std::exit(1);
    }

    // Append index record
    const PageIdxRec rec{shape, page_id, off};
    if (std::fwrite(&rec, 1, sizeof(rec), ctx.idx) != sizeof(rec)) {
        std::perror("fwrite pages.idx"); std::exit(1);
    }

    if (ctx.solved) {
        std::fseek(ctx.solved, (long)off, SEEK_SET);
        std::array<std::uint8_t, PAGE_BYTES> zeros{};
        if (std::fwrite(zeros.data(), 1, PAGE_BYTES, ctx.solved) != PAGE_BYTES) {
            std::perror("fwrite pages.solved"); std::exit(1);
        }
    }

    return off;
}

// ---------- cache management (single-threaded) ----------

inline PageBuf* find_cached(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    auto it = ctx.cache.find(Key{shape,page_id});
    return (it==ctx.cache.end()) ? nullptr : &it->second;
}

inline PageBuf& ensure_cached(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    const Key k{shape, page_id};

    if (PageBuf* p = find_cached(ctx, shape, page_id)) return *p;
    const std::uint64_t off = ensure_page_exists_fast(ctx, shape, page_id);

    // Insert and initialize from disk
    auto [it, inserted] = ctx.cache.try_emplace(k);
    PageBuf* ref = &it->second;
    if (inserted) {
        std::array<std::uint64_t, PAGE_QWORDS> tmp;
        std::fseek(ctx.data, (long)off, SEEK_SET);
        if (std::fread(tmp.data(), 1, PAGE_BYTES, ctx.data) != PAGE_BYTES) {
            std::perror("fread pages.dat"); std::exit(1);
        }
        std::memcpy(ref->w.data(), tmp.data(), sizeof(ref->w));
        ref->dirty = false;

        // Evict if over capacity (~25% of cache), prefer clean pages first
        if (ctx.cache.size() > ctx.cache_cap_pages) {
            const std::size_t target = ctx.cache.size() / 4;
            std::size_t evicted = 0;

            // First pass: erase clean pages
            for (auto e = ctx.cache.begin(); e != ctx.cache.end() && evicted < target; ) {
                if (e->first.shape == k.shape && e->first.page == k.page) { ++e; continue; }
                if (!e->second.dirty) { e = ctx.cache.erase(e); ++evicted; }
                else ++e;
            }

            // Second pass: flush and erase some dirty pages if needed
            for (auto e = ctx.cache.begin(); e != ctx.cache.end() && evicted < target; ) {
                if (e->first.shape == k.shape && e->first.page == k.page) { ++e; continue; }
                if (e->second.dirty) {
                    flush_page(ctx, e->first.shape, e->first.page);
                    e = ctx.cache.erase(e);
                    ++evicted;
                } else ++e;
            }
        }
    }
    return *ref;
}

inline bool set_bit(Context& ctx,
                    std::uint64_t shape,
                    std::uint32_t page_id,
                    std::uint32_t bit_in_page)
{
    PageBuf& pb = ensure_cached(ctx, shape, page_id);

    const std::uint32_t wi = bit_in_page >> 6;
    const std::uint32_t bi = bit_in_page & 63u;
    const std::uint64_t m  = 1ULL << bi;

    const std::uint64_t old = pb.w[wi];
    const bool was_new = (old & m) == 0;
    if (was_new) { pb.w[wi] = old | m; pb.dirty = true; }
    return was_new;
}

// ---------- solved values storage ----------

inline void write_solved_value(Context& ctx, std::uint64_t shape, std::uint32_t page_id, 
                               std::uint32_t bit_in_page, std::uint8_t value) {
    if (!ctx.solved) {
        std::fprintf(stderr, "Error: solved file not open\n");
        return;
    }
    
    const std::uint64_t offset = offset_of(ctx, shape, page_id);
    if (offset == ~0ULL) {
        std::fprintf(stderr, "Error: page not found for shape %llu, page %u\n", shape, page_id);
        return; // page doesn't exist
    }
    
    const std::uint64_t solved_offset = offset + bit_in_page;
    if (std::fseek(ctx.solved, (long)solved_offset, SEEK_SET) != 0) {
        std::perror("fseek solved file");
        return;
    }
    if (std::fwrite(&value, 1, 1, ctx.solved) != 1) {
        std::perror("fwrite solved value");
        return;
    }
}

inline std::uint8_t read_solved_value(Context& ctx, std::uint64_t shape, std::uint32_t page_id, 
                                     std::uint32_t bit_in_page) {
    if (!ctx.solved) {
        std::fprintf(stderr, "Error: solved file not open\n");
        return 0;
    }
    
    const std::uint64_t offset = offset_of(ctx, shape, page_id);
    if (offset == ~0ULL) {
        std::fprintf(stderr, "Error: page not found for shape %llu, page %u\n", shape, page_id);
        return 0; // page doesn't exist
    }
    
    const std::uint64_t solved_offset = offset + bit_in_page;
    if (std::fseek(ctx.solved, (long)solved_offset, SEEK_SET) != 0) {
        std::perror("fseek solved file");
        return 0;
    }
    std::uint8_t value;
    if (std::fread(&value, 1, 1, ctx.solved) != 1) {
        std::perror("fread solved value");
        return 0;
    }
    return value;
}

inline void solve_pos(Context& ctx,
                    std::uint64_t shape,
                    std::uint32_t page_id,
                    std::uint32_t bit_in_page,
                    const std::function<Context&(std::uint8_t)>& get_ctx_func) 
{
    // Safety checks
    if (bit_in_page >= (1u << PAGE_BITS)) {
        std::fprintf(stderr, "Error: bit_in_page %u exceeds PAGE_BITS %zu\n", bit_in_page, PAGE_BITS);
        return;
    }
    
    std::uint64_t restored_hash = (page_id << PAGE_BITS) | bit_in_page;
    game::Encoding pos = game::unhash(static_cast<bitops6::Bitboard>(shape), restored_hash);
    if (std::uint8_t prim = game::primitive(pos)) {
        write_solved_value(ctx, shape, page_id, bit_in_page, prim);
        return;
    }

    std::uint64_t moves = game::legal_moves(pos);
    std::vector<uint8_t> children; 

    if (moves) {
        while (moves) {
            std::uint64_t mv = moves & -moves;
            moves ^= mv;
            const game::Encoding nxt = game::flip(game::do_move(pos, mv));
            const game::Encoding c  = game::canonical(nxt);
            const std::uint64_t sh = game::shape(c);
            const std::uint64_t h  = game::hash(c);
            const std::uint8_t tier = game::tier_of(sh);
            std::printf("Debug: Accessing tier %d for child position\n", tier);
            Context& child_ctx = get_ctx_func(tier);

            std::uint32_t curr_page = static_cast<std::uint32_t>(h >> PAGE_BITS);
            ensure_page_exists_fast(child_ctx, sh, curr_page);
            std::uint32_t curr_bit  = static_cast<std::uint32_t>(h & ((1u << PAGE_BITS) - 1u));
            std::uint8_t child = read_solved_value(child_ctx, sh, curr_page, curr_bit);
            children.push_back(child);
        }
    } else {
        const game::Encoding flipped = game::flip(pos);
        const game::Encoding c  = game::canonical(flipped);
        const std::uint64_t sh = game::shape(c);
        const std::uint64_t h  = game::hash(c);
        const std::uint8_t tier = game::tier_of(sh);
        std::printf("Debug: Accessing tier %d for pass move\n", tier);
        Context& child_ctx = get_ctx_func(tier);

        std::uint32_t curr_page = static_cast<std::uint32_t>(h >> PAGE_BITS);
        ensure_page_exists_fast(child_ctx, sh, curr_page);
        std::uint32_t curr_bit  = static_cast<std::uint32_t>(h & ((1u << PAGE_BITS) - 1u));
        std::uint8_t child = read_solved_value(child_ctx, sh, curr_page, curr_bit);
        children.push_back(child);
    }

    // Value and Remoteness logic
    std::uint8_t v_and_r;
    if (bitops6::all_losing(children)) {
        v_and_r = game::WIN | (bitops6::find_closest_prim(children) + 1);
    } else if (bitops6::tie_exists(children)) {
        v_and_r = game::DRAW | (bitops6::find_farthest_prim(children) + 1);
    } else {
        v_and_r = game::LOSE | (bitops6::find_farthest_prim(children) + 1);
    }
    
    write_solved_value(ctx, shape, page_id, bit_in_page, v_and_r);
    std::printf("Restored hash: %llu\n", (unsigned long long)restored_hash);
    std::printf("Value and remoteness: %d\n", v_and_r);
    return;
}

inline void solve_page(Context& ctx, const std::function<Context&(std::uint8_t)>& get_ctx_func) {
    std::fseek(ctx.data, 0, SEEK_SET);
    std::uint64_t shape, word_buf = 0, wi = 0; // shape, word buffer, word index

    while (true) {
        if (std::fseek(ctx.idx, (wi >> (PAGE_BITS - 6)) * sizeof(PageIdxRec), SEEK_SET) != 0) break;
        if (std::fread(&shape, sizeof(shape), 1, ctx.idx) != 1) {
            if (std::feof(ctx.idx)) break; 
            std::perror("fread shape");
            std::exit(1);
        }
        if (std::fread(&word_buf, sizeof(word_buf), 1, ctx.data) != 1) {
            if (std::feof(ctx.data)) break; 
            std::perror("fread data");
            std::exit(1);
        }

        std::uint64_t bi = 0;
        while (word_buf) {
            while (!(word_buf & 1)) {
                word_buf >>= 1;
                bi++; 
            }
            std::uint32_t pid = wi >> (PAGE_BITS - 6); // page id
            std::uint64_t bip = (wi << 6) + bi; // bit in page
            solve_pos(ctx, shape, pid, bip & ((1 << PAGE_BITS) - 1), get_ctx_func);
            word_buf >>= 1;
            bi++;
        }
        wi++;
    }
}

} // namespace pageops
