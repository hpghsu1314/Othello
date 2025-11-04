#pragma once
#include <mpi.h>
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
#include "message.hpp"
#include "utilities.hpp"

namespace pageops {

namespace fs = std::filesystem;

static constexpr std::size_t PAGE_BITS   = 12;
static constexpr std::size_t PAGE_BYTES  = (1u << (PAGE_BITS - 3));              // 2^(PAGE_BITS)/8
static constexpr std::size_t PAGE_QWORDS = PAGE_BYTES / sizeof(std::uint64_t);
static constexpr std::uint8_t REMOTENESS_BITS = 5;
static constexpr int TAG_MODULE_TO_MAIN = 2;
static constexpr int TAG_MAIN_TO_MODULE = 3;

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

template <typename Key>
uint64_t* nested_find(MetaMap& map, std::uint64_t outer, const Key& inner_key)
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
    std::FILE* solved = nullptr; // pages.solved
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

    auto it_outer = ctx.meta_map.find(game::tier_of(shape));   // look for existing outer key
    if (it_outer != ctx.meta_map.end()) {       // only if outer exists
        auto& inner_map = it_outer->second;     // reference to inner map
        inner_map.try_emplace(k, off);          // insert inner key if missing
    } else {
        ctx.meta_map.emplace(game::tier_of(shape), MetaMap::mapped_type{{k, off}});
    }

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
        std::fseek(ctx.solved, (long)(off << 3), SEEK_SET);
        std::array<std::uint8_t, PAGE_BYTES> zeros{};
        if (std::fwrite(zeros.data(), 1, (PAGE_BYTES << 3), ctx.solved) != (PAGE_BYTES << 3)) {
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
        std::fprintf(stderr, "Error: solved file not open on write\n");
        return;
    }
    
    const std::uint64_t page_offset = offset_of(ctx, shape, page_id);

    printf("write page offset: %llu with bit in page %u\n", page_offset, bit_in_page);
    const std::uint64_t solved_offset = (page_offset << 3) + bit_in_page;
    if (std::fseek(ctx.solved, (long)solved_offset, SEEK_SET) != 0) {
        std::perror("fseek solved file");
        return;
    }

    printf("Writing file: %s \n", ctx.solved_path.c_str());
    printf("Writing file (ptr): %p \n", ctx.solved);
    printf("Writing this value: %d \n", value);
    printf("Writing at this offset: %llu \n", solved_offset);
    printf("write side -- shape: %llu, page: %u, bit: %u \n", shape, page_id, bit_in_page);

    if (std::fwrite(&value, 1, 1, ctx.solved) != 1) {
        std::perror("fwrite solved value");
        return;
    }

    printf("sucessful write!! \n");
}

inline std::uint8_t read_solved_value(std::uint64_t shape, std::uint32_t page_id, 
                                     std::uint32_t bit_in_page, Context& ctx) {
    if (!ctx.solved) {
        std::fprintf(stderr, "Error: solved file not open on read\n");
        return 0;
    }
    
    const std::uint64_t page_offset = offset_of(ctx, shape, page_id);

    const std::uint64_t solved_offset = (page_offset << 3) + bit_in_page;

    printf("Reading file: %s \n", ctx.solved_path.c_str());
    printf("Reading file size: %lu \n", ctx.file_size);
    printf("offset: %llu \n", solved_offset);
    
    if (std::fseek(ctx.solved, (long)solved_offset, SEEK_SET) != 0) {
        std::perror("fseek solved file");
        return 0;
    }

    std::uint8_t value;
    if (std::fread(&value, 1, 1, ctx.solved) != 1) {
        std::perror("fread solved value");
        return 0;
    }

    printf("fread succesfull! The value is: %d \n", value);

    return value;
}

inline MPI_Datatype make_TierMsg_type() {
    MPI_Datatype T;
    const int N = 4;
    int          bl[N]  = {1,1,1,1};
    MPI_Aint     dis[N] = {};
    MPI_Datatype ty[N]  = {MPI_UINT8_T, MPI_UINT64_T, MPI_UINT32_T, MPI_UINT32_T};

    TierMsg probe{};
    MPI_Aint base=0, a0=0, a1=0, a2=0, a3=0;
    MPI_Get_address(&probe,        &base);
    MPI_Get_address(&probe.tier,  &a0);
    MPI_Get_address(&probe.shape,  &a1);
    MPI_Get_address(&probe.page,   &a2);
    MPI_Get_address(&probe.bit,    &a3);
    dis[0] = a0 - base; dis[1] = a1 - base; dis[2] = a2 - base, dis[3] = a3 - base;

    MPI_Type_create_struct(N, bl, dis, ty, &T);
    MPI_Type_commit(&T);
    return T;
}

inline TierMsg make_TierMsg_from_pos(const game::Encoding& position) noexcept {
    TierMsg m{};
    const auto c  = game::canonical(position);
    const auto sh = game::shape(c);
    const auto h  = game::hash(c);
    const auto tier = game::tier_of(sh);
    m.shape = sh;
    m.page  = static_cast<std::uint32_t>(h >> pageops::PAGE_BITS);
    m.bit   = static_cast<std::uint32_t>(h & ((1u << pageops::PAGE_BITS) - 1u));
    m.tier = tier;
    return m;
}

inline void solve_pos(Context (&ctx_by_tier)[game::MAX_TIER + 1],
                    std::uint64_t shape,
                    std::uint32_t page_id,
                    std::uint32_t bit_in_page,
                    std::uint8_t tier, 
                    MPI_Comm comm) 
{
    Context& ctx = ctx_by_tier[tier];
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Safety checks
    if (bit_in_page >= (1u << PAGE_BITS)) {
        std::fprintf(stderr, "Error: bit_in_page %u exceeds PAGE_BITS %zu\n", bit_in_page, PAGE_BITS);
        return;
    }
    
    std::uint64_t restored_hash = (page_id << PAGE_BITS) | bit_in_page;
    game::Encoding pos = game::unhash(static_cast<bitops6::Bitboard>(shape), restored_hash);
    printf("Game position (player): %llu \n", pos.player);
    if (std::uint8_t prim = game::primitive(pos)) {
        printf("stats before solving -- shape: %llu, page: %u, bit: %u\n", shape, page_id, bit_in_page);
        write_solved_value(ctx, shape, page_id, bit_in_page, prim);
        return;
    }

    std::uint64_t moves = game::legal_moves(pos);
    int num_children = (moves ? std::popcount(moves) : 1);
    std::vector<std::uint8_t> children(num_children);
    std::vector<MPI_Request> child_send_requests(num_children);
    std::vector<MPI_Request> child_req_requests(num_children);
    std::vector<TierMsg> outgoing_msgs(num_children);
    std::vector<int> skipped_msgs(num_children);
    int head = 0;
    MPI_Datatype TIER_T = make_TierMsg_type(); 

    if (moves) {
        printf("Number of moves is %d", std::popcount(moves));
        int dest;
        while (moves) {
            std::uint64_t mv = moves & -moves;
            moves ^= mv;
            const game::Encoding nxt = game::flip(game::do_move(pos, mv));
            TierMsg child = make_TierMsg_from_pos(nxt);
            printf("send side -- shape: %llu, page: %u, bit: %u, tier: %d \n", child.shape, child.page, child.bit, child.tier);
            dest = owner_of_shape(child.shape, size);

            outgoing_msgs[head] = child;
            if (dest == rank) {
                if (child.tier == tier) {
                    skipped_msgs.push_back(head);
                } else {
                    children[head] = read_solved_value(child.shape, child.page, child.bit, ctx_by_tier[child.tier]);
                }
                child_send_requests[head] = child_req_requests[head] = MPI_REQUEST_NULL;
                head++;
                continue;
            }

            MPI_Send(outgoing_msgs.data() + head, 1, TIER_T, dest, TAG_MODULE_TO_MAIN, comm);
            printf("reached past the send request blocking \n"); 
            MPI_Irecv(children.data() + head, 1, MPI_UINT8_T, dest, TAG_MAIN_TO_MODULE, comm, child_req_requests.data() + head);
            head++;
        }
    } else {
        const game::Encoding flipped = game::flip(pos);
        TierMsg child = make_TierMsg_from_pos(flipped);
        int dest = owner_of_shape(child.shape, size);

        outgoing_msgs[head] = child;
        if (dest == rank) {
            if (child.tier == tier) {
                skipped_msgs.push_back(head);
            } else {
                children[head] = read_solved_value(child.shape, child.page, child.bit, ctx_by_tier[child.tier]);
            }
            child_send_requests[head] = child_req_requests[head] = MPI_REQUEST_NULL;
        } else {
            MPI_Send(outgoing_msgs.data() + head, 1, TIER_T, dest, TAG_MODULE_TO_MAIN, comm);
            MPI_Irecv(children.data() + head, 1, MPI_UINT8_T, dest, TAG_MAIN_TO_MODULE, comm, child_req_requests.data() + head);
        }
    }
    MPI_Waitall(child_req_requests.size(), child_req_requests.data(), MPI_STATUSES_IGNORE);
    printf("reached past waitall \n");
    for (int idx : skipped_msgs) {
        TierMsg neighbor = outgoing_msgs[idx];
        children[idx] = read_solved_value(neighbor.shape, neighbor.page, neighbor.bit, ctx_by_tier[tier]);
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

inline void solve_page(Context (&ctx_by_tier)[game::MAX_TIER + 1], std::uint8_t tier, MPI_Comm comm) {
    Context& ctx = ctx_by_tier[tier];
    auto it_outer = ctx.meta_map.find(tier);

    if (it_outer != ctx.meta_map.end()) {
        auto& inner_map = it_outer->second; 
        for (const auto& [inner_key, value] : inner_map) {
            std::uint64_t word_buf;

            std::fseek(ctx.data, value, SEEK_SET);
            for (uint8_t wi = 0; wi < (1 << (PAGE_BITS - 6)); wi++) {
                if (std::fread(&word_buf, sizeof(word_buf), 1, ctx.data) != 1) break;

                if (word_buf) {
                    std::uint8_t bi = 0;
                    while (word_buf) {
                        while (!(word_buf & 1)) {
                            word_buf >>= 1;
                            bi++; 
                        }
                        std::uint32_t bip = (wi << 6) + bi; 
                        solve_pos(ctx_by_tier, inner_key.shape, inner_key.page, bip, tier, comm);
                        word_buf >>= 1;
                        bi++;
                    }
                } 
            }
        }
    }
}

} // namespace pageops

/*

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

┌───┬───┬───┬───┬───┬───┐
│ P │ P │ P │ P │ P │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ P │ O │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ P │ O │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ P │ O │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ O │ O │ O │ P │ P │
├───┼───┼───┼───┼───┼───┤
│ P │ P │ P │ P │ P │ P │
└───┴───┴───┴───┴───┴───┘

*/