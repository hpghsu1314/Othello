#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <filesystem>
#include <random>
#include <cassert>
#include <omp.h>
#include "page_ops.hpp"

namespace fs = std::filesystem;

#define REQ(cond, msg) do { if(!(cond)){ std::fprintf(stderr,"FAIL @%s:%d: %s\n",__FILE__,__LINE__, msg); std::exit(2);} } while(0)

static std::uint64_t fsize(const std::string& p) {
    return fs::exists(p) ? (std::uint64_t)fs::file_size(p) : 0ULL;
}

static void read_bytes(const std::string& path, std::uint64_t off, void* dst, std::size_t n) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::perror("open"); std::exit(1); }
    if (std::fseek(f, (long)off, SEEK_SET) != 0) { std::perror("fseek"); std::exit(1); }
    if (std::fread(dst, 1, n, f) != n) { std::perror("fread"); std::exit(1); }
    std::fclose(f);
}

#pragma pack(push,1)
struct PageIdxRec { std::uint64_t shape; std::uint32_t page; std::uint64_t off; };
#pragma pack(pop)

static std::vector<PageIdxRec> read_all_idx(const std::string& idx_path) {
    std::vector<PageIdxRec> v;
    if (!fs::exists(idx_path)) return v;
    const std::uint64_t sz = fsize(idx_path);
    REQ(sz % sizeof(PageIdxRec) == 0, "pages.idx size not multiple of record");
    v.resize(sz / sizeof(PageIdxRec));
    read_bytes(idx_path, 0, v.data(), (std::size_t)sz);
    return v;
}

int main() {
    // Unique output dir per run (so idx size checks are clean)
    const int owner = 0;
    const unsigned tier = 77; // pick a throwaway tier
    const std::string root = "out_stresstest";
    const std::string dir  = root + "/owner_00/tier_77";
    fs::create_directories(dir);
    const std::string dpth = dir + "/pages.dat";
    const std::string ipth = dir + "/pages.idx";

    // Clean any leftovers from a previous failed run
    if (fs::exists(dpth)) fs::remove(dpth);
    if (fs::exists(ipth)) fs::remove(ipth);

    std::printf("PAGE_BYTES = %zu (set PAGE_BITS in page_ops.hpp to shrink for tests)\n", pageops::PAGE_BYTES);

    pageops::Context ctx;
    ctx.cache_cap_pages = 128;  // keep small to force eviction; scale if you want
    pageops::pages_init(ctx, owner, tier, root);

    // ---------- Test 0: baseline single-thread ----------
    {
        const std::uint64_t A = 0xAAAALLu;
        auto off0 = pageops::ensure_page_exists_fast(ctx, A, 0);
        auto off1 = pageops::ensure_page_exists_fast(ctx, A, 1);
        REQ(off1 % pageops::PAGE_BYTES == 0, "off1 align");
        REQ(off0 % pageops::PAGE_BYTES == 0, "off0 align");

        // bit 0 and last bit
        bool n1 = pageops::set_bit(ctx, A, 0, 0);
        bool n2 = pageops::set_bit(ctx, A, 0, (1u<<pageops::PAGE_BITS) - 1);
        bool n1dup = pageops::set_bit(ctx, A, 0, 0);
        REQ(n1 && n2 && !n1dup, "basic set_bit");

        pageops::flush_all(ctx);

        const unsigned char first = [&]{
            unsigned char b=0; read_bytes(dpth, off0 + 0, &b, 1); return b;
        }();
        const unsigned char last  = [&]{
            unsigned char b=0; read_bytes(dpth, off0 + pageops::PAGE_BYTES-1, &b, 1); return b;
        }();
        std::printf("T0: first=0x%02X last=0x%02X (expect 01,80)\n", first, last);
        REQ(first == 0x01 && last == 0x80, "disk bytes mismatch");
    }

    // ---------- Test 1: creation race on same (shape,page) ----------
    {
        const std::uint64_t SHAPE = 0xC0FFEEu;
        const std::uint32_t PAGE  = 42;
        const int T = std::max(8, omp_get_max_threads()*2);
        std::vector<std::uint64_t> offs(T);

        #pragma omp parallel for schedule(static)
        for (int i=0;i<T;++i) {
            offs[i] = pageops::ensure_page_exists_fast(ctx, SHAPE, PAGE);
        }

        // All offsets must be identical
        for (int i=1;i<T;++i) REQ(offs[i]==offs[0], "creation race produced different offsets");

        pageops::flush_all(ctx);

        // pages.idx should have exactly one record for (SHAPE,PAGE)
        auto recs = read_all_idx(ipth);
        int cnt = 0;
        for (auto &r: recs) if (r.shape==SHAPE && r.page==PAGE) ++cnt;
        REQ(cnt == 1, "idx recorded duplicate creations for same key");
        std::printf("T1: creation race OK (one idx record, stable offset)\n");
    }

    // ---------- Test 2: concurrent bit sets on SAME page (data race) ----------
    {
        const std::uint64_t SHAPE = 0xDEADBEAFu;
        const std::uint32_t PAGE  = 7;
        auto off = pageops::ensure_page_exists_fast(ctx, SHAPE, PAGE);

        // Expected bitmap (word-wise OR), using atomics to build reference under concurrency
        std::vector<std::atomic<std::uint64_t>> expected(pageops::PAGE_BYTES/8);
        for (auto &a: expected) a.store(0, std::memory_order_relaxed);

        const int T = std::max(8, omp_get_max_threads()*2);
        const int OPS_PER_THREAD = 5000;  // scale up as you like

        #pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            std::mt19937_64 rng(0x12345678ULL + tid*0x9e3779b97f4a7c15ULL);
            std::uniform_int_distribution<std::uint32_t> bitdist(0, (1u<<pageops::PAGE_BITS)-1);

            #pragma omp for schedule(static)
            for (int t=0; t<T; ++t) {
                for (int k=0;k<OPS_PER_THREAD;++k) {
                    const std::uint32_t bit = bitdist(rng);
                    // write to file-backed page
                    (void)pageops::set_bit(ctx, SHAPE, PAGE, bit);
                    // build reference
                    const std::uint32_t wi = bit >> 6, bi = bit & 63u;
                    const std::uint64_t m = 1ULL << bi;
                    expected[wi].fetch_or(m, std::memory_order_relaxed);
                }
            }
        }

        pageops::flush_all(ctx);

        // Read back and compare whole page
        std::vector<std::uint64_t> back(expected.size());
        read_bytes(dpth, off, back.data(), pageops::PAGE_BYTES);

        for (std::size_t i=0;i<back.size();++i) {
            const auto exp = expected[i].load(std::memory_order_relaxed);
            REQ(back[i] == exp, "concurrent set_bit content mismatch");
        }
        std::printf("T2: concurrent set_bit on same page OK\n");
    }

    // ---------- Test 3: cache pressure & eviction (including dirty flush) ----------
    {
        // Make cache tiny so we force eviction/flush
        ctx.cache_cap_pages = 32;

        const std::uint64_t SHAPE = 0xFACEB00Cu;
        const int PAGES = 32 * 6; // 6x over cap
        std::vector<std::uint64_t> offs(PAGES);

        // Create many pages in parallel
        #pragma omp parallel for schedule(static)
        for (int p=0;p<PAGES;++p) {
            offs[p] = pageops::ensure_page_exists_fast(ctx, SHAPE, (std::uint32_t)p);
            // dirty them immediately: set 3 bits
            pageops::set_bit(ctx, SHAPE, (std::uint32_t)p, 0);
            pageops::set_bit(ctx, SHAPE, (std::uint32_t)p, 11);
            pageops::set_bit(ctx, SHAPE, (std::uint32_t)p, (1u<<pageops::PAGE_BITS)-1);
        }

        // Force further churn: touch more distinct pages so eviction triggers
        #pragma omp parallel for schedule(static)
        for (int p=PAGES; p<PAGES+128; ++p) {
            pageops::ensure_page_exists_fast(ctx, SHAPE, (std::uint32_t)p);
        }

        // Flush and verify a few random victims were persisted
        pageops::flush_all(ctx);

        // Sample and check bytes on disk for several pages spread across range
        auto check_page = [&](int p) {
            auto off = pageops::offset_of(ctx, SHAPE, (std::uint32_t)p);
            REQ(off != ~0ULL, "missing offset after churn");
            unsigned char b0=0, b1=0, bL=0;
            read_bytes(dpth, off + 0, &b0, 1);
            read_bytes(dpth, off + (11>>3), &b1, 1);
            read_bytes(dpth, off + pageops::PAGE_BYTES - 1, &bL, 1);
            REQ((b0 & 0x01) != 0, "bit0 not set after eviction+flush");
            REQ((b1 & (1u << (11&7))) != 0, "bit11 not set after eviction+flush");
            REQ((bL & 0x80) != 0, "last bit not set after eviction+flush");
        };

        for (int p : {0, 7, 31, 63, 95, PAGES-1}) check_page(p);
        std::printf("T3: cache pressure + eviction + dirty flush OK\n");
    }

    // ---------- Test 4: many shapes/pages in parallel; consistency + idx sanity ----------
    {
        const int SHAPES = 200;
        const int PAGES_PER_SHAPE = 20;

        std::vector<std::uint64_t> shapes(SHAPES);
        for (int i=0;i<SHAPES;++i) shapes[i] = 0xBEEF0000ULL + (std::uint64_t)i;

        #pragma omp parallel for schedule(static)
        for (int s=0;s<SHAPES;++s) {
            for (int p=0;p<PAGES_PER_SHAPE;++p) {
                (void)pageops::ensure_page_exists_fast(ctx, shapes[s], (std::uint32_t)p);
                // sprinkle a few bits
                pageops::set_bit(ctx, shapes[s], (std::uint32_t)p, (p*3) & ((1u<<pageops::PAGE_BITS)-1));
            }
        }

        pageops::flush_all(ctx);

        auto recs = read_all_idx(ipth);
        // Count how many records for our shape block (they must exist; exactly one per created page)
        std::size_t found = 0;
        for (auto &r : recs) {
            if (r.shape >= shapes.front() && r.shape <= shapes.back()) ++found;
        }
        REQ(found >= (std::size_t)SHAPES*PAGES_PER_SHAPE, "idx missing records (parallel create)");
        std::printf("T4: parallel create across many shapes/pages OK\n");
    }

    // ---------- Test 5: repeated ensure/set do NOT grow files ----------
    {
        const auto dat_before = fsize(dpth);
        const auto idx_before = fsize(ipth);

        const std::uint64_t S = 0x12345678ULL;
        for (int i=0;i<5000;++i) {
            (void)pageops::ensure_page_exists_fast(ctx, S, 0);
            (void)pageops::ensure_page_exists_fast(ctx, S, 1);
            (void)pageops::set_bit(ctx, S, 0, 0);
        }
        pageops::flush_all(ctx);

        const auto dat_after = fsize(dpth);
        const auto idx_after = fsize(ipth);
        REQ(dat_after == dat_before, "pages.dat grew on repeats");
        REQ(idx_after == idx_before, "pages.idx grew on repeats");
        std::printf("T5: repeat ensure/set doesn't bloat files OK\n");
    }

    pageops::pages_close(ctx);
    puts("\nALL PAGE OPS TESTS PASSED.");
    return 0;
}
