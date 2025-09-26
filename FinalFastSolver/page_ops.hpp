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

namespace pageops {

namespace fs = std::filesystem;

static constexpr std::size_t PAGE_BITS   = 12;
static constexpr std::size_t PAGE_BYTES  = (1u << (PAGE_BITS - 3));              // 2^(PAGE_BITS)/8
static constexpr std::size_t PAGE_QWORDS = PAGE_BYTES / sizeof(std::uint64_t);

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
    std::size_t file_size = 0;  // current end of pages.dat (bytes)

    // Paths
    std::string data_path;
    std::string idx_path;

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

    char dir[512], dpath[512], ipath[512];
    std::snprintf(dir,   sizeof(dir),   "%s/owner_%02d/tier_%02u", root_dir.c_str(), owner_num, tier);
    std::snprintf(dpath, sizeof(dpath), "%s/pages.dat", dir);
    std::snprintf(ipath, sizeof(ipath), "%s/pages.idx", dir);
    fs::create_directories(dir);

    ctx.data_path = dpath;
    ctx.idx_path  = ipath;

    ctx.data = std::fopen(ctx.data_path.c_str(), "rb+");
    if (!ctx.data) ctx.data = std::fopen(ctx.data_path.c_str(), "wb+");
    if (!ctx.data) { std::perror("fopen pages.dat"); std::exit(1); }
    std::setvbuf(ctx.data, nullptr, _IOFBF, 1<<20); // 1 MiB

    ctx.idx = std::fopen(ctx.idx_path.c_str(), "ab+");
    if (!ctx.idx) { std::perror("fopen pages.idx"); std::exit(1); }
    std::setvbuf(ctx.idx, nullptr, _IOFBF, 1<<20);

    std::fseek(ctx.data, 0, SEEK_END);
    ctx.file_size = (std::size_t)std::ftell(ctx.data);

    ctx.inited = true;
}

inline void pages_flush(Context& ctx) {
    if (ctx.data) std::fflush(ctx.data);
    if (ctx.idx)  std::fflush(ctx.idx);
}

inline void pages_close(Context& ctx) {
    pages_flush(ctx);
    if (ctx.data) { std::fclose(ctx.data); ctx.data=nullptr; }
    if (ctx.idx)  { std::fclose(ctx.idx);  ctx.idx=nullptr;  }
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

} // namespace pageops
