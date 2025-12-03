#pragma once
#include <mpi.h>
#include <cstdint>
#include <cstdio>
#include <list>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <array>
#include <vector>
#include <deque>
#include <filesystem>
#include <functional>
#include "othello6.hpp" // must be changed for 4x4
#include "bitops6.hpp"
#include "message.hpp"
#include "utilities.hpp"
#include <zstd.h> 

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

using MetaMap  = std::unordered_map<std::uint64_t, std::unordered_map<Key, std::uint64_t, KeyHash>>; 
using CacheMap = std::unordered_map<Key, PageBuf,   KeyHash>;

class LRUCache {
public:
    struct Entry {
        Key key;
        std::array<std::uint8_t, (PAGE_BYTES << 3)> buf;
    };

    unsigned int capacity;
    std::list<Entry> order;
    std::unordered_map<Key, std::list<Entry>::iterator, KeyHash> lookup;

    LRUCache(unsigned int capacity)
        : capacity(capacity) {}

    std::array<std::uint8_t, (PAGE_BYTES << 3)>* find(Key key) {
        auto it = lookup.find(key);
        if (it == lookup.end()) return nullptr;

        order.splice(order.begin(), order, it->second);
        return &it->second->buf;
    }

    void put(Key key, const std::array<std::uint8_t, (PAGE_BYTES << 3)>& buf) {
        auto it = lookup.find(key);
        if (it != lookup.end()) {
            it->second->buf = buf;
            order.splice(order.begin(), order, it->second);
            return;
        }

        if (order.size() == capacity) {
            auto lru = std::prev(order.end());
            lookup.erase(lru->key);
            order.erase(lru);
        }

        order.push_front({key, buf});
        lookup[key] = order.begin();
    }
};


template <typename Key>
std::uint64_t* nested_find(MetaMap& map, std::uint64_t outer, const Key& inner_key)
{
    auto it_outer = map.find(outer);
    if (it_outer == map.end()) return nullptr;

    auto& inner_map = it_outer->second;
    auto it_inner = inner_map.find(inner_key); 
    if (it_inner == inner_map.end()) return nullptr;

    return &it_inner->second;
}

struct Context {
    // Hot file handles (kept open)
    std::FILE* data = nullptr;  // pages.dat
    std::FILE* idx  = nullptr;  // pages.idx
    std::FILE* solved = nullptr; // pages.solved, per tier
    std::size_t file_size = 0;  // current end of pages.dat (bytes)

    // Paths
    std::string data_path;
    std::string idx_path;
    std::string solved_path;

    // In-RAM maps
    MetaMap  meta_map;          // tier -> (shape, page) -> offset
    CacheMap cache;             // (shape,page) -> page buffer (in-place, no moves/copies)

    // Cache size knob (pages)
    std::size_t cache_cap_pages = 4096u * 16u;

    // Read-only LRU cache for solving
    LRUCache rec_cache;             // (shape, page) -> offset (no flushing)

    bool inited{false};
    Context() : rec_cache(4096u * 16u) {}
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

inline MPI_File init_db_write_only(const std::string& root_dir, const std::string& suffix, int tier, MPI_Comm comm, int mpi_mode = MPI_MODE_WRONLY) {
    // Create data directory
    char gpath[512];
    std::snprintf(gpath, sizeof(gpath), "%s/tier_%02u_db/tier.%s", root_dir.c_str(), tier, suffix.c_str());

    // Ensure path exists on all ranks
    fs::path dir = fs::path(gpath).parent_path();
    fs::create_directories(dir);
    
    // Use MPI File utilities for file operations 
    MPI_File fh;
    int rc = MPI_File_open(comm, gpath, mpi_mode, MPI_INFO_NULL, &fh);
    if (rc != MPI_SUCCESS) {
        char errstr[MPI_MAX_ERROR_STRING];
        int resultlen;
        MPI_Error_string(rc, errstr, &resultlen);
        fprintf(stderr, "Rank: MPI_File_open failed: %s\n", errstr);
        MPI_Abort(comm, rc);
    }
    return fh;
}

inline std::FILE* rank_init_db_write_only(const std::string& root_dir, const std::string& suffix, int tier, int rank) {
    // Create data directory
    char gpath[512];
    std::snprintf(gpath, sizeof(gpath), "%s/owner_%02u/tier_%02u_db/tier.%s", root_dir.c_str(), rank, tier, suffix.c_str());

    // Ensure path exists on all ranks
    fs::path dir = fs::path(gpath).parent_path();
    fs::create_directories(dir);
    
    std::FILE* fp = std::fopen(gpath, "wb");
    if (!fp) {
        std::perror("fopen");
    }
    return fp;
}

inline std::FILE* init_db_read_only(const std::string& root_dir, const std::string& suffix, int tier) {
    // Ensure data directory exists
    fs::create_directories(root_dir);
    
    char gpath[512];
    std::snprintf(gpath, sizeof(gpath), "%s/tier_%02u_db/tier.%s", root_dir.c_str(), tier, suffix.c_str());
    
    // File already exists from below tiers
    std::FILE* fp = std::fopen(gpath, "rb");
    if (!fp) {
        std::perror("fopen");
    }
    return fp;
}

inline void db_close(MPI_File& fh) {
    MPI_File_close(&fh);
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
        const auto* mit = nested_find(ctx.meta_map, game::tier_of(shape), k);
        if (!mit) { std::fprintf(stderr, "flush_page: missing meta\n"); std::exit(1); }
        const std::uint64_t off = *mit;

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
    if (ctx.solved) std::fflush(ctx.solved);
}

inline std::uint64_t offset_of(Context& ctx, std::uint64_t shape, std::uint32_t page_id) {
    const Key k{shape,page_id};
    const auto* it = nested_find(ctx.meta_map, game::tier_of(shape), k);
    return (!it) ? ~0ULL : *it;
}

// ---------- ensure page exists (single-threaded), return offset ----------

inline std::uint64_t ensure_page_exists_fast(Context& ctx,
                                             std::uint64_t shape,
                                             std::uint32_t page_id)
{
    const Key k{shape,page_id};

    const auto* it = nested_find(ctx.meta_map, game::tier_of(shape), k);
    if (it) return *it;

    // Reserve unique offset and publish mapping
    const std::uint64_t off = ctx.file_size;
    ctx.file_size += PAGE_BYTES;

    auto it_outer = ctx.meta_map.find(game::tier_of(shape));   
    if (it_outer != ctx.meta_map.end()) {       
        auto& inner_map = it_outer->second;     
        inner_map.try_emplace(k, off);          
    } else ctx.meta_map.emplace(game::tier_of(shape), MetaMap::mapped_type{{k, off}});

    // Write zero page at the reserved offset
    std::fseek(ctx.data, (long)off, SEEK_SET);
    if (std::fwrite(ZERO_PAGE, 1, PAGE_BYTES, ctx.data) != PAGE_BYTES) {
        std::perror("fwrite pages.dat"); std::exit(1);
    }

    std::fseek(ctx.solved, (long)(off << 3), SEEK_SET);
    std::array<std::uint8_t, PAGE_BYTES> extended_zero_page{};
    if (std::fwrite(extended_zero_page.data(), 1, (PAGE_BYTES << 3), ctx.solved) != (PAGE_BYTES << 3)) {
        std::perror("fwrite pages.solved"); std::exit(1);
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

inline TierMsg make_TierMsg_from_pos(const game::Encoding& position) noexcept {
    TierMsg m{};
    const auto c  = game::canonical(position);
    const auto sh = game::shape(c);
    const auto h  = game::hash(c);
    const auto tier = game::tier_of(sh);
    m.tier = tier;
    m.shape = sh;   
    m.page  = static_cast<std::uint32_t>(h >> pageops::PAGE_BITS);
    m.bit   = static_cast<std::uint32_t>(h & ((1u << pageops::PAGE_BITS) - 1u));
    m.hash = h;
    return m;
}

inline std::uint8_t load_into_cache(Context& ctx, 
                            TierMsg& child, 
                            int num_processes,
                            std::FILE* data_fp, 
                            std::FILE* rec_fp, 
                            std::vector<std::vector<std::uint64_t>>& data_lookup,
                            std::vector<std::vector<std::uint64_t>>& rec_lookup) 
{
    const Key k{child.shape, child.page};
    std::array<std::uint8_t, (PAGE_BYTES << 3)>* buffer = ctx.rec_cache.find(k);
    if (buffer) return (*buffer)[child.bit];

    // Find page respective offset
    int dest = owner_of_shape(child.shape, num_processes);
    std::fseek(rec_fp, rec_lookup[child.tier][dest], SEEK_SET);
    PageIdxRec rec;
    while (std::fread(&rec, sizeof(rec), 1, rec_fp) == 1) {
        if (rec.shape == child.shape && rec.page == child.page) break;
    }

    // load page into cache (lookup + page offset)
    std::array<std::uint8_t, (PAGE_BYTES << 3)> newBuf;

    std::fseek(data_fp, data_lookup[child.tier][dest] + (rec.off << 3), SEEK_SET);
    if (std::fread(newBuf.data(), 1, (PAGE_BYTES << 3), data_fp) != (PAGE_BYTES << 3)) {
        std::perror("fread child value");
        return 0;
    }

    ctx.rec_cache.put(k, newBuf);
    return newBuf[child.bit];
}

inline std::uint8_t get_prim(std::uint8_t val) {
    std::uint8_t prim = val >> 6, rem = val & ((1 << 6) - 1), res = val + 1;
    if (prim == (game::WIN >> 6)) res = game::LOSE | (rem + 1);
    if (prim == (game::LOSE >> 6)) res = game::WIN | (rem + 1);
    return res;
}

inline void write_solved_value(Context& ctx, std::uint64_t shape, std::uint32_t page_id, 
                               std::uint32_t bit_in_page, std::uint8_t value) {
    if (!ctx.solved) {
        std::fprintf(stderr, "Error: solved file not open on write\n");
        return;
    }
    
    const std::uint64_t page_offset = offset_of(ctx, shape, page_id);
    const std::uint64_t solved_offset = (page_offset << 3) + bit_in_page;

    if (std::fseek(ctx.solved, (long)solved_offset, SEEK_SET) != 0) {
        std::perror("fseek solved file");
        return;
    }

    if (std::fwrite(&value, 1, 1, ctx.solved) != 1) {
        std::perror("fwrite solved value");
        return;
    }

    printf("successful write \n");
}

inline void solve_pos(Context& ctx,
                      std::uint64_t shape,
                      std::uint32_t page_id,
                      std::uint32_t bit_in_page,
                      std::uint8_t tier,
                      int num_processes,
                      std::FILE* data_fp,
                      std::FILE* rec_fp,
                      std::vector<std::vector<std::uint64_t>>& data_lookup,
                      std::vector<std::vector<std::uint64_t>>& rec_lookup) 
{
    std::uint64_t restored_hash = (page_id << PAGE_BITS) | bit_in_page;
    game::Encoding pos = game::unhash(static_cast<bitops6::Bitboard>(shape), restored_hash);
    if (std::uint8_t prim = game::primitive(pos)) {
        if (tier == 35) {
            printf("player: %llu, opponent: %llu \n", pos.player, pos.opponent);
        }
        write_solved_value(ctx, shape, page_id, bit_in_page, prim); 
        return;
    }

    std::uint64_t moves = game::legal_moves(pos);
    int num_children = moves ? std::popcount(moves): 1;
    std::vector<std::uint8_t> children(num_children);

    if (moves) {
        while (moves) {
            std::uint64_t mv = moves & -moves;
            moves ^= mv;
            game::Encoding nxt = game::flip(game::do_move(pos, mv));
            if (!game::primitive(nxt) && !game::legal_moves(nxt)) nxt = game::flip(nxt);
            TierMsg child = make_TierMsg_from_pos(nxt);

            std::uint8_t child_val = load_into_cache(ctx, child, num_processes, data_fp, rec_fp, data_lookup, rec_lookup);
            children.push_back(child_val); 
        }
    } else return;

    // Value and Remoteness logic
    std::uint8_t v_and_r;
    std::vector<std::uint8_t> res = game::points_of_interest(children);
    if (res[1] == 0xFF && res[2] == 0xFF) v_and_r = game::LOSE | (res[0] + 1);
    else if (res[2] != 0xFF) v_and_r = game::WIN | (res[2] + 1);
    else v_and_r = game::DRAW | (res[1] + 1);
    
    write_solved_value(ctx, shape, page_id, bit_in_page, v_and_r);
    std::printf("Restored hash: %llu\n", (unsigned long long)restored_hash);
    std::printf("The value and remoteness: %lu for tier %d \n", (long)v_and_r, tier);
}

inline void solve_tier(Context& ctx, 
                       std::uint8_t tier, 
                       int rank, 
                       int num_processes,
                       MPI_File& tier_db, 
                       MPI_File& tier_rec, 
                       std::vector<std::vector<std::uint64_t>>& data_lookup,
                       std::vector<std::vector<std::uint64_t>>& rec_lookup,
                       MPI_Comm comm) 
{
    auto it_outer = ctx.meta_map.find(tier);

    std::FILE* next_tier = nullptr;
    std::FILE* next_rec = nullptr;
    bool accessed = false;
    if (tier != game::MAX_TIER) {
        next_tier = init_db_read_only("data", "dat", tier+1);
        next_rec = init_db_read_only("data", "idx", tier+1);
        accessed = true;
    }

    if (it_outer != ctx.meta_map.end()) {
        auto& inner_map = it_outer->second; 
        for (const auto& [inner_key, value] : inner_map) {
            std::uint64_t word_buf;

            std::fseek(ctx.data, value, SEEK_SET);
            for (std::uint8_t wi = 0; wi < (1 << (PAGE_BITS - 6)); wi++) {
                if (std::fread(&word_buf, sizeof(word_buf), 1, ctx.data) != 1) break;

                if (word_buf) {
                    std::uint8_t bi = 0;
                    while (word_buf) {
                        while (!(word_buf & 1)) {
                            word_buf >>= 1;
                            bi++; 
                        }
                        std::uint32_t bip = (wi << 6) + bi; 
                        solve_pos(ctx, inner_key.shape, inner_key.page, bip, tier, num_processes, next_tier, next_rec, data_lookup, rec_lookup);
                        word_buf >>= 1;
                        bi++;
                    }
                } 
            }
        }
    }

    // page.solved -> tier_db | page.idx -> tier_rec
    size_t data_size = 0, rec_size = 0; 
    std::vector<std::uint8_t> data_buf, rec_buf;
    if (ctx.solved) {
        fseek(ctx.solved, 0, SEEK_END); fseek(ctx.idx, 0, SEEK_END);
        data_size = ftell(ctx.solved); rec_size = ftell(ctx.idx);
        data_buf.resize(data_size); rec_buf.resize(rec_size);
        rewind(ctx.solved); rewind(ctx.idx);

        if (fread(data_buf.data(), 1, data_size, ctx.solved) != data_size) perror("Failed to read solved.dat");
        if (fread(rec_buf.data(), 1, rec_size, ctx.idx) != rec_size) perror("Failed to read idx.dat");
    }

    MPI_Offset data_size_mpi = static_cast<MPI_Offset>(data_size);
    MPI_Offset rec_size_mpi = static_cast<MPI_Offset>(rec_size);

    MPI_Offset off_data = 0, off_rec = 0;
    MPI_Exscan(&data_size_mpi, &off_data, 1, MPI_OFFSET, MPI_SUM, comm);
    MPI_Exscan(&rec_size_mpi, &off_rec, 1, MPI_OFFSET, MPI_SUM, comm);

    if (rank == 0) {
        off_data = 0;
        off_rec = 0;
    }

    MPI_File_write_ordered(tier_db, data_buf.data(), data_size, MPI_BYTE, MPI_STATUS_IGNORE);
    MPI_File_write_ordered(tier_rec, rec_buf.data(), rec_size, MPI_BYTE, MPI_STATUS_IGNORE);

    data_lookup[tier][rank] = off_data;
    rec_lookup[tier][rank] = off_rec;

    if (accessed) {
        fclose(next_tier); fclose(next_rec);
    }
}

// Method 1:
inline std::vector<std::uint8_t> serialize_key_map(std::uint8_t* page_bytes) {
    std::vector<std::uint8_t> buffer;

    for (std::uint16_t i = 0; i < (1 << PAGE_BITS); i++) {
        if (page_bytes[i]) {
            // serialize key (i) as two bytes, big-endian
            buffer.push_back(static_cast<std::uint8_t>(i >> 8));
            buffer.push_back(static_cast<std::uint8_t>(i & 0xFF));
            buffer.push_back(page_bytes[i]);
        }
    }

    return buffer;
}

// Method 2:
inline std::array<std::uint8_t, 4096> impute_invalid_values(std::uint8_t* page_bytes) {
    std::array<std::uint8_t, 4096> result;
    result[0] = page_bytes[0];
    bool is_valid = page_bytes[0] ? true : false;

    for (std::uint16_t i = 1; i < (1 << PAGE_BITS); i++) {
        if (!page_bytes[i] && is_valid) result[i] = result[i-1];
        else {
            if (page_bytes[i]) is_valid = true; 
            result[i] = page_bytes[i];
        }
    }
    return result;
}

struct Node {
    std::uint8_t value;
    std::uint16_t left;
    std::uint16_t right;
};

// Method 3:
inline std::vector<Node> merge_page_bytes(std::uint8_t* page_bytes) {
    std::vector<Node> curr_depth;
    int end_index = 1 << PAGE_BITS;
    for (std::uint16_t i = 0; i < end_index; i++) {
        Node byte_node = {page_bytes[i], i, i};
        curr_depth.push_back(byte_node);
    }

    std::uint16_t front_index = 0;
    int num_merged = 0;

    while (true) {
        Node curr_node = curr_depth[front_index++];

        if (curr_node.right == end_index - 1) {
            curr_depth.push_back(curr_node);
            if (!num_merged) break;
            else {num_merged = 0; continue;}
        }

        Node adj_node = curr_depth[front_index++];

        if ((curr_node.right - curr_node.left == adj_node.right - adj_node.left) && (curr_node.value == adj_node.value)) {
            Node new_node = {curr_node.value, curr_node.left, adj_node.right};
            curr_depth.push_back(new_node);
            num_merged++;
        } else {
            curr_depth.push_back(curr_node);
            curr_depth.push_back(adj_node);
            if (adj_node.right == end_index - 1 && !num_merged) break;
        }

        if (adj_node.right == end_index - 1) num_merged = 0;
    }
    
    curr_depth.erase(curr_depth.begin(), curr_depth.begin() + front_index);
    return curr_depth;
}

std::vector<uint8_t> serialize_nodes(const std::vector<Node>& nodes) {
    std::vector<uint8_t> buffer;
    buffer.reserve(nodes.size() * 5); 

    for (const auto& n : nodes) {
        buffer.push_back(n.value);     
        buffer.push_back(n.left >> 8);      // left high byte     
        buffer.push_back(n.left & 0xFF);    // left low byte
        buffer.push_back(n.right >> 8);     // right high byte
        buffer.push_back(n.right & 0xFF);   // right low byte
    }
    return buffer;
}

inline void compress_all_pages(int compressable_tier, int rank, std::uint64_t file_offset, std::uint64_t start_rec_offset, std::uint64_t end_rec_offset) {
    std::FILE* unused_tier_db = init_db_read_only("data", "dat", compressable_tier); 
    std::FILE* unused_rec_db = init_db_read_only("data", "idx", compressable_tier);
    std::FILE* new_tier_db = rank_init_db_write_only("data", "dat", compressable_tier, rank);
    std::FILE* new_rec_db = rank_init_db_write_only("data", "idx", compressable_tier, rank);
    size_t running_offset = 0;

    if (!end_rec_offset) {
        std::fseek(unused_rec_db, start_rec_offset, SEEK_SET);

        PageIdxRec rec;
        while (true) {
            size_t nread = std::fread(&rec, sizeof(rec), 1, unused_rec_db);
            if (nread != 1) {
                if (feof(unused_rec_db)) break;       
                if (ferror(unused_rec_db)) {     
                    perror("Error reading rec_db");
                    break;
                }
            }

            std::fseek(unused_tier_db, file_offset + rec.off, SEEK_SET);

            std::uint8_t buf[4096];
            std::fread(buf, 1, 4096, unused_tier_db);

            std::vector<std::uint8_t> m1 = serialize_key_map(buf);
            std::array<std::uint8_t, 4096> m2 = impute_invalid_values(buf);
            std::vector<std::uint8_t> m3 = serialize_nodes(merge_page_bytes(buf));

            size_t m1_size_max = ZSTD_compressBound(m1.size());
            size_t m2_size_max = ZSTD_compressBound(m2.size());
            size_t m3_size_max = ZSTD_compressBound(m3.size());

            std::vector<uint8_t> m1_c(m1_size_max), m2_c(m2_size_max), m3_c(m3_size_max);

            size_t m1_csize = ZSTD_compress(m1_c.data(), m1_size_max, m1.data(), m1.size(), 3);
            size_t m2_csize = ZSTD_compress(m2_c.data(), m2_size_max, m2.data(), m2.size(), 3);
            size_t m3_csize = ZSTD_compress(m3_c.data(), m3_size_max, m3.data(), m3.size(), 3);

            size_t minimum = std::min({m1_csize, m2_csize, m3_csize});
            if (minimum == m1_csize) {
                m1_c.resize(m1_csize);
                rec.off = (1ULL << 62) | running_offset;
                running_offset += m1_csize;

                std::fwrite(m1_c.data(), 1, m1_c.size(), new_tier_db);
                std::fwrite(&rec, sizeof(rec), 1, new_rec_db);
            } else if (minimum == m2_csize) {
                m2_c.resize(m2_csize);
                rec.off = (2ULL << 62) | running_offset;
                running_offset += m2_csize;

                std::fwrite(m2_c.data(), 1, m2_c.size(), new_tier_db);
                std::fwrite(&rec, sizeof(rec), 1, new_rec_db);
            } else {
                m3_c.resize(m3_csize);
                rec.off = (3ULL << 62) | running_offset;
                running_offset += m3_csize;

                std::fwrite(m3_c.data(), 1, m3_c.size(), new_tier_db);
                std::fwrite(&rec, sizeof(rec), 1, new_rec_db);
            }
        }
    } else {
        std::uint64_t num_iterations = (end_rec_offset - start_rec_offset) >> 12;
        std::fseek(unused_rec_db, start_rec_offset, SEEK_SET);
        for (std::uint64_t i = 0; i < num_iterations; i ++) {
            PageIdxRec rec;
            std::fread(&rec, sizeof(rec), 1, unused_rec_db);
            std::fseek(unused_tier_db, file_offset + rec.off, SEEK_SET);

            std::uint8_t buf[4096];
            std::fread(buf, 1, 4096, unused_tier_db);

            // Try all compression schemes
            std::vector<std::uint8_t> m1 = serialize_key_map(buf);
            std::array<std::uint8_t, 4096> m2 = impute_invalid_values(buf);
            std::vector<std::uint8_t> m3 = serialize_nodes(merge_page_bytes(buf)); 

            size_t m1_size_max = ZSTD_compressBound(m1.size()), m2_size_max = ZSTD_compressBound(m2.size()), m3_size_max = ZSTD_compressBound(m3.size());
            std::vector<uint8_t> m1_c(m1_size_max), m2_c(m2_size_max), m3_c(m3_size_max);
            
            // Compression level 3, can be tuned
            size_t m1_csize = ZSTD_compress(m1_c.data(), m1_size_max, m1.data(), m1.size(), 3);
            size_t m2_csize = ZSTD_compress(m2_c.data(), m2_size_max, m2.data(), m2.size(), 3);
            size_t m3_csize = ZSTD_compress(m3_c.data(), m3_size_max, m3.data(), m3.size(), 3);

            size_t minimum = std::min({m1_csize, m2_csize, m3_csize});
            if (minimum == m1_csize) {
                m1_c.resize(m1_csize);
                rec.off = (1 << 62) | running_offset;
                running_offset += m1_csize;
                
                std::fwrite(m1_c.data(), 1, m1_c.size(), new_tier_db);
                std::fwrite(&rec, sizeof(rec), 1, new_rec_db);
            } else if (minimum == m2_csize) {
                m2_c.resize(m2_csize);
                rec.off = (2 << 62) | running_offset;
                running_offset += m2_csize;
                
                std::fwrite(m2_c.data(), 1, m2_c.size(), new_tier_db);
                std::fwrite(&rec, sizeof(rec), 1, new_rec_db);
            } else {
                m3_c.resize(m3_csize);
                rec.off = (3 << 62) | running_offset;
                running_offset += m3_csize;
                
                std::fwrite(m3_c.data(), 1, m3_c.size(), new_tier_db);
                std::fwrite(&rec, sizeof(rec), 1, new_rec_db);
            }
        }
    }

    std::fflush(new_tier_db); std::fflush(new_rec_db);
    std::fclose(unused_tier_db); std::fclose(unused_rec_db);
    std::fclose(new_tier_db); std::fclose(new_rec_db);
}

} // namespace pageops

/*

Test 1:
┌───┬───┬───┬───┬───┬───┐
│ . │ P │ P │ P │ P │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ P │ P │ P │ P │ P │
└───┴───┴───┴───┴───┴───┘

Test 2: 
┌───┬───┬───┬───┬───┬───┐
│ . │ P │ P │ P │ P │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ . │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ O │ . │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ . │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ P │ P │ P │ P │ . │
└───┴───┴───┴───┴───┴───┘   

*/