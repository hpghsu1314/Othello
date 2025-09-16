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
#include <atomic>
#include <omp.h>

namespace pageops {

namespace fs = std::filesystem;

static constexpr std::size_t PAGE_BITS      = 12;                           // 2^16 bits/page
static constexpr std::size_t PAGE_BYTES     = (1u << (PAGE_BITS - 3));      // 8 KiB
static constexpr std::size_t PAGE_QWORDS    = PAGE_BYTES / sizeof(std::uint64_t);

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
    alignas(64) std::array<std::uint64_t, PAGE_QWORDS> w;
    std::atomic_flag lock;
    bool dirty;
    PageBuf() : w{}, lock(ATOMIC_FLAG_INIT), dirty(false) {}    // default-constructible; no copy/move
    PageBuf(const PageBuf&)            = delete;
    PageBuf& operator=(const PageBuf&) = delete;
};

using MetaMap  = std::unordered_map<Key, std::uint64_t, KeyHash>;
using CacheMap = std::unordered_map<Key, PageBuf, KeyHash>;

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
    std::size_t cache_cap_pages = 64*1024;

    // OpenMP locks
    omp_lock_t meta_lock;
    omp_lock_t cache_lock;
    omp_lock_t data_lock;
    omp_lock_t idx_lock;

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

    omp_init_lock(&ctx.meta_lock);
    omp_init_lock(&ctx.cache_lock);
    omp_init_lock(&ctx.data_lock);
    omp_init_lock(&ctx.idx_lock);

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
    if (ctx.inited) {
        omp_destroy_lock(&ctx.meta_lock);
        omp_destroy_lock(&ctx.cache_lock);
        omp_destroy_lock(&ctx.data_lock);
        omp_destroy_lock(&ctx.idx_lock);
        ctx.inited = false;
    }
}

// ---------- flush (call outside parallel regions) ----------

inline void flush_page(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    // Snapshot pointer
    omp_set_lock(&ctx.cache_lock);
    auto it = ctx.cache.find(Key{shape,page_id});
    PageBuf* p = (it==ctx.cache.end()) ? nullptr : &it->second;
    omp_unset_lock(&ctx.cache_lock);
    if (!p) return;

    // Serialize with writers
    while (p->lock.test_and_set(std::memory_order_acquire)) { /* spin */ }
    if (p->dirty) {
        std::uint64_t off;
        omp_set_lock(&ctx.meta_lock);
        off = ctx.meta_map[Key{shape,page_id}];
        omp_unset_lock(&ctx.meta_lock);

        omp_set_lock(&ctx.data_lock);
        std::fseek(ctx.data, (long)off, SEEK_SET);
        if (std::fwrite(p->w.data(), 1, PAGE_BYTES, ctx.data) != PAGE_BYTES) {
            std::perror("fwrite pages.dat"); std::exit(1);
        }
        omp_unset_lock(&ctx.data_lock);
        p->dirty = false;
    }
    p->lock.clear(std::memory_order_release);
}

inline void flush_all(Context& ctx) {
    // snapshot keys
    std::vector<Key> keys;
    omp_set_lock(&ctx.cache_lock);
    keys.reserve(ctx.cache.size());
    for (auto &kv : ctx.cache) keys.push_back(kv.first);
    omp_unset_lock(&ctx.cache_lock);

    for (const Key& k : keys) flush_page(ctx, k.shape, k.page);
    if (ctx.data) std::fflush(ctx.data);
    if (ctx.idx)  std::fflush(ctx.idx);
}

inline std::uint64_t offset_of(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    omp_set_lock(&ctx.meta_lock);
    auto it = ctx.meta_map.find(Key{shape,page_id});
    const std::uint64_t off = (it==ctx.meta_map.end()) ? ~0ULL : it->second;
    omp_unset_lock(&ctx.meta_lock);
    return off;
}

// ---------- ensure page exists (thread-safe), return offset ----------

inline std::uint64_t ensure_page_exists_fast(Context& ctx,
                                             std::uint64_t shape,
                                             std::uint32_t page_id)
{
    const Key k{shape,page_id};

    // Fast path under meta lock
    omp_set_lock(&ctx.meta_lock);
    if (auto it = ctx.meta_map.find(k); it != ctx.meta_map.end()) {
        const std::uint64_t off = it->second;
        omp_unset_lock(&ctx.meta_lock);
        return off;
    }

    // Reserve a unique offset and PUBLISH the mapping immediately,
    // so other threads see it and don't try to create again.
    const std::uint64_t off = ctx.file_size;
    ctx.file_size += PAGE_BYTES;
    ctx.meta_map.try_emplace(k, off);
    omp_unset_lock(&ctx.meta_lock);

    // Write the page at the RESERVED offset (not append), so offset is stable.
    omp_set_lock(&ctx.data_lock);
    std::fseek(ctx.data, (long)off, SEEK_SET);
    if (std::fwrite(ZERO_PAGE, 1, PAGE_BYTES, ctx.data) != PAGE_BYTES) {
        std::perror("fwrite pages.dat"); std::exit(1);
    }
    omp_unset_lock(&ctx.data_lock);

    // Append one index record
    omp_set_lock(&ctx.idx_lock);
    const PageIdxRec rec{shape, page_id, off};
    if (std::fwrite(&rec, 1, sizeof(rec), ctx.idx) != sizeof(rec)) {
        std::perror("fwrite pages.idx"); std::exit(1);
    }
    omp_unset_lock(&ctx.idx_lock);

    return off;
}

// ---------- cache management (thread-safe) ----------

inline PageBuf* find_cached(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    omp_set_lock(&ctx.cache_lock);
    auto it = ctx.cache.find(Key{shape,page_id});
    PageBuf* ptr = (it==ctx.cache.end()) ? nullptr : &it->second;
    omp_unset_lock(&ctx.cache_lock);
    return ptr;
}

inline PageBuf& ensure_cached(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    const Key k{shape, page_id};

    // Fast probe
    if (PageBuf* p = find_cached(ctx, shape, page_id)) return *p;

    // Ensure exists on disk and get its offset (unique, published)
    const std::uint64_t off = ensure_page_exists_fast(ctx, shape, page_id);

    // Insert (or find existing) in cache
    omp_set_lock(&ctx.cache_lock);
    auto [it, inserted] = ctx.cache.try_emplace(k);   // default-constructs PageBuf if missing
    PageBuf* ref = &it->second;
    if (inserted) {
        // Pin page so concurrent set_bit will wait until load is done
        while (ref->lock.test_and_set(std::memory_order_acquire)) { /* spin */ }
    }
    omp_unset_lock(&ctx.cache_lock);

    if (inserted) {
        // Load from disk without holding cache_lock
        std::array<std::uint64_t, PAGE_QWORDS> tmp;
        omp_set_lock(&ctx.data_lock);
        std::fseek(ctx.data, (long)off, SEEK_SET);
        if (std::fread(tmp.data(), 1, PAGE_BYTES, ctx.data) != PAGE_BYTES) {
            std::perror("fread pages.dat"); std::exit(1);
        }
        omp_unset_lock(&ctx.data_lock);

        // Copy into the pinned page; use sizeof(ref->w) to keep analyzer happy
        std::memcpy(ref->w.data(), tmp.data(), sizeof(ref->w));
        ref->dirty = false;
        ref->lock.clear(std::memory_order_release);

        // Evict (possibly) after the page is fully initialized
        // Evict ~25% but never the page we just loaded
        std::size_t evicted = 0;
        std::size_t target  = 0;

        // First pass: clean victims
        omp_set_lock(&ctx.cache_lock);
        if (ctx.cache.size() > ctx.cache_cap_pages) {
            target = ctx.cache.size() / 4;
            for (auto e = ctx.cache.begin(); e != ctx.cache.end() && evicted < target; ) {
                if (e->first.shape == k.shape && e->first.page == k.page) { ++e; continue; }
                if (!e->second.dirty) { e = ctx.cache.erase(e); ++evicted; }
                else ++e;
            }
        }
        omp_unset_lock(&ctx.cache_lock);

        if (evicted < target) {
            // Second pass: choose dirty victims, flush, then evict
            const std::size_t need = target - evicted;
            std::vector<Key> victims; victims.reserve(need);

            omp_set_lock(&ctx.cache_lock);
            for (auto &kv : ctx.cache) {
                if (kv.first.shape == k.shape && kv.first.page == k.page) continue;
                if (kv.second.dirty) {
                    victims.push_back(kv.first);
                    if (victims.size() == need) break;
                }
            }
            omp_unset_lock(&ctx.cache_lock);

            for (const Key& vk : victims) flush_page(ctx, vk.shape, vk.page);

            omp_set_lock(&ctx.cache_lock);
            for (const Key& vk : victims) {
                auto e = ctx.cache.find(vk);
                if (e != ctx.cache.end() && !e->second.dirty) ctx.cache.erase(e);
            }
            omp_unset_lock(&ctx.cache_lock);
        }
    }

    // Return current entry (re-find under lock to avoid stale iterators)
    omp_set_lock(&ctx.cache_lock);
    PageBuf& out = ctx.cache.find(k)->second;
    omp_unset_lock(&ctx.cache_lock);
    return out;
}



// ---------- bit set (thread-safe per page). Returns true if 0->1 ----------

inline bool set_bit(Context& ctx,
                    std::uint64_t shape,
                    std::uint32_t page_id,
                    std::uint32_t bit_in_page)
{
    PageBuf& pb = ensure_cached(ctx, shape, page_id);

    // Per-page spin
    while (pb.lock.test_and_set(std::memory_order_acquire)) { /* spin */ }

    const std::uint32_t wi = bit_in_page >> 6;
    const std::uint32_t bi = bit_in_page & 63u;
    const std::uint64_t m  = 1ULL << bi;

    const std::uint64_t old = pb.w[wi];
    const bool was_new = (old & m) == 0;
    if (was_new) { pb.w[wi] = old | m; pb.dirty = true; }

    pb.lock.clear(std::memory_order_release);
    return was_new;
}

} // namespace pageops
