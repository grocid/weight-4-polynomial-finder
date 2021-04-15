#include <iostream>
#include <iomanip>
#include <sstream>
#include <thread>

#include <fstream>

#include <execution>
#include <algorithm>

#include <boost/functional/hash.hpp>
#include "parallel_hashmap/phmap.h"


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "config.h"

#include "stringops.cpp"
#include "masks.hpp"

#define LAYER_SIZE (THREADS*THREADS)

// Chunk size is the total number of elements contained in memory
// at the same time, which equates to how much RAM is available.
#define CHUNK_SIZE (1 << 25)

using namespace std;

static constexpr uint64_t total_map_size = (uint64_t)(((uint64_t)1 << (polynomial_degree/3 + BETA)) * ALPHA);
static constexpr uint64_t exp_len = (get_degree(total_map_size) + 8) / 8;

PUREFUN inline constexpr uint32_t phi(uint128_t val)
{
    return val % THREADS;
}

#ifdef __GNUC__
#define PACKSTRUCT __attribute__((__packed__))
#else
#define PACKSTRUCT
#endif

typedef struct PACKSTRUCT exp_t {
    uint8_t exponent[exp_len];
} exp_t;

typedef struct PACKSTRUCT mask_t {
    uint8_t mask[masklenb];
} mask_t;

typedef struct PACKSTRUCT imask_t {
    uint8_t imask[imasklenb];
} imask_t;

typedef struct PACKSTRUCT cmap_poly {
    imask_t imask;
    exp_t exponent;
} cmap_poly;

typedef struct PACKSTRUCT clay_poly {
    exp_t e1;
    exp_t e2;
} clay_poly;

PUREFUN constexpr inline exp_t pack_exp(uint64_t exp) {
    exp_t rv = {};
    for (uint i = 0; i < exp_len; i++)
        rv.exponent[i] = exp >> (8*i);
    return rv;
}


PUREFUN constexpr inline uint64_t unpack_exp(const exp_t &exp) {
    uint64_t rv = 0;
    for (uint i = 0; i < exp_len; i++)
        rv |= ((uint64_t)exp.exponent[i]) <<  (8*i);
    return rv;
}

PUREFUN constexpr inline mask_t pack_mask(uint128_t mask) {
    mask_t rv = {};
    for (uint i = 0; i < masklenb; i++)
        rv.mask[i] = mask >> (8*i);
    return rv;
}


PUREFUN constexpr inline uint128_t unpack_mask(const mask_t &mask) {
    uint128_t rv = 0;
    for (uint i = 0; i < masklenb; i++)
        rv |= ((uint128_t)mask.mask[i]) <<  (8*i);
    return rv;
}

PUREFUN constexpr inline imask_t pack_imask(uint128_t imask) {
   imask_t rv = {};
    for (uint i = 0; i < imasklenb; i++)
        rv.imask[i] = imask >> (8*i);
    return rv;
}


PUREFUN constexpr inline uint128_t unpack_imask(const imask_t &imask) {
    uint128_t rv = 0;
    for (uint i = 0; i < imasklenb; i++)
        rv |= ((uint128_t)imask.imask[i]) <<  (8*i);
    return rv;
}

PUREFUN constexpr inline imask_t xor_imask(const imask_t &imask1, const imask_t &imask2) {
    imask_t rv = {};
    for (uint i = 0; i < imasklenb; i++)
        rv.imask[i] = imask1.imask[i] ^ imask2.imask[i];
    return rv;
}


namespace std {
    template<>
    struct hash<mask_t>
    {
        typedef mask_t argument_type;
        typedef size_t result_type;
        PUREFUN constexpr inline result_type operator () (const argument_type& x) const
        {
            result_type rv = 0;
            for (uint i = 0; i < masklenb; i++) {
                rv = rv << 8 | rv >> (sizeof(result_type)-1*8);
                rv ^= x.mask[i];
            }
            return rv;
        }
    };
    template<>
    struct hash<imask_t>
    {
        typedef imask_t argument_type;
        typedef size_t result_type;
        PUREFUN constexpr inline result_type operator () (const argument_type& x) const
        {
            result_type rv = 0;
            for (uint i = 0; i < imasklenb; i++) {
                rv = rv << 8 | rv >> (sizeof(result_type)-1*8);
                rv ^= x.imask[i];
            }
            return rv;
        }
    };
}

PUREFUN constexpr inline bool operator==(const mask_t& lhs, const mask_t& rhs)
{
    return !memcmp(lhs.mask, rhs.mask, masklenb);
}

PUREFUN constexpr inline bool operator==(const imask_t& lhs, const imask_t& rhs)
{
    return !memcmp(lhs.imask, rhs.imask, imasklenb);
}

#define map_t phmap::flat_hash_map
#define EXTRAARGS phmap::container_internal::hash_default_hash<K>, \
                  phmap::container_internal::hash_default_eq<K>, \
                  std::allocator<std::pair<const K, V>>, 32, phmap::NullMutex

template <class K, class V> using layer = map_t<K, V>;

map_t<imask_t, clay_poly> collision_layer[LAYER_SIZE];


void in_memory_generate(uint128_t mask, uint32_t thread)
{
    uint64_t exponent = 0;
    uint32_t idx;
#ifdef DEBUG_MESSAGES
    uint32_t added = 0;
#endif
    uint128_t px = 1;
    map_t<mask_t, cmap_poly> collision_map;
    
    for(uint64_t i = 0; i < total_map_size; ++i)
    {
        // every thread considers each monomial in the range
        // but it will only work with monomials that match
        // the phi condition.
        uint128_t mpx = get_mask_bits(px);
        idx = phi(mpx);

        if (unlikely(idx == thread))
        {
            exp_t pexp = pack_exp(exponent);
            mask_t mbits = pack_mask(mpx);
            imask_t imbits = pack_imask(get_imask_bits(px));
            auto [it, result] = collision_map.try_emplace(
                mbits, cmap_poly{imbits, pexp}
            );

            if(unlikely(!result))
            {
                imask_t imaskxor = xor_imask(it->second.imask, imbits);
                uint128_t py = unpack_imask(imaskxor);
                exp_t exponent2 = it->second.exponent;
                auto [it2, result] = collision_layer[THREADS * thread + phi(py)].try_emplace(
                    imaskxor, clay_poly{pexp, exponent2}
                );
#ifdef DEBUG_MESSAGES
                added++;
#endif
                // This is rather unlikely to happen, but if we don't check
                // then we might miss a multiple.
                if(unlikely(!result))
                {
                    vector<uint64_t> exponents;
                    exponents.push_back(unpack_exp(pexp));
                    exponents.push_back(unpack_exp(exponent2));
                    exponents.push_back(unpack_exp(it2->second.e1));
                    exponents.push_back(unpack_exp(it2->second.e2));

                    cout << polynomial_representation(exponents).str() << endl;
                }
            }
        }

        px = next_monomial(px, POLY, polynomial_degree);
        exponent++;
    }
#ifdef DEBUG_MESSAGES
    cout << "\tThread " << thread  << " produced " << added << " in first step." << endl;
#endif
}

void in_memory_merge(uint128_t mask, int thread)
{
    auto base = collision_layer[THREADS * 0 + thread];
 
    for(int i = 1; i < THREADS; ++i)
    {
        for (auto& it: collision_layer[THREADS * i + thread])
        {
            auto [it2, result] = base.try_emplace(
                it.first, it.second
            );

            if(!result)
            {
                vector<uint64_t> exponents;
                exponents.push_back(unpack_exp(it.second.e1));
                exponents.push_back(unpack_exp(it.second.e2));
                exponents.push_back(unpack_exp(it2->second.e1));
                exponents.push_back(unpack_exp(it2->second.e2));

                cout << polynomial_representation(exponents).str() << endl;
            }
        }
        collision_layer[THREADS * i + thread].clear();
    }
}


void in_memory_search (uint128_t mask)
{
    cout << endl << "Running in-memory (square algorithm) search..." << endl;
#if THREADS > 1
    cout << "[1/2]\t";
#endif

    cout << "Generating 2^" << log2(total_map_size) << " monomials..." << endl;
    cout << "Using " << (exp_len*8) << " bits for exponents..." << endl;
    thread t[THREADS];    
    for (uint32_t i = 0; i < THREADS; ++i) {
        t[i] = thread(in_memory_generate, mask, i);
    }

    for (uint32_t i = 0; i < THREADS; ++i) {
        t[i].join();
    }

#if THREADS > 1
    cout << "[2/2]\tMerging binomials..." << endl << endl;

    for (uint32_t i = 0; i < THREADS; ++i) {
        t[i] = thread(in_memory_merge, mask, i);
    }

    for (uint32_t i = 0; i < THREADS; ++i) {
        t[i].join();
    }
#endif
}

typedef tuple<uint128_t, uint64_t, uint64_t> element_t;

bool comparePolynomialMask(element_t i1, element_t i2) 
{ 
    return ((get<0>(i1) & mask) < (get<0>(i2) & mask)); 
} 

/*

Idea: 

1. use the phi condition to determine subsets
2. each thread generates N entries
3. entries are split into different files based on subsets, hence every file generates M files.
4. files belonging to a thread are concatenated
5. the concatenated files are sorted individually by threads
6. using matching, new entries are generated and thus subsets based on phi condition
7. the files are picked up by corresponding thread and concatenated
8. final merge can be done individually

*/


void write(string fname, vector<element_t> data)
{
    ofstream ofile;

    ofile.open(fname, ios::binary);
    if (!data.empty())
    {
        ofile.write(
            reinterpret_cast<char*>(&data[0]), 
            data.size() * sizeof(element_t)
        );
    }
}

vector<element_t> merge_files(string fname[])
{

    /*
        Open all files as streams
        Use iterator and add all chunks to minheap with elements <element_t, chunk>
        Pop element to output stream
        Read additional element to minheap from corresponding filestream
    */
}


void sequential_generate(uint128_t mask, uint64_t thread, uint64_t starting_exponent, uint32_t chunk)
{

    uint64_t idx;
    uint64_t exponent = starting_exponent;
    uint128_t px = gf2x_exp(2, exponent, POLY);
    uint32_t f = polynomial_degree;

    vector<element_t> sequential;

    for(uint64_t i = 0; i < CHUNK_SIZE; ++i)
    {
        idx = phi(px & mask);

        if (idx == thread)
        {
            sequential.push_back({px, exponent, 0});
        }

        px = next_monomial(px, POLY, f);
        exponent++;
    }

    sort(sequential.begin(), sequential.end(), comparePolynomialMask);
    write(generate_fname(thread, chunk).str(), sequential);
}

void thread_generate(uint128_t mask)
{
    cout << endl << "Running sequential search..." << endl;

#if THREADS > 1
    cout << "[1/2]\t";
#endif

    cout << "Generating 2^" << log2(total_map_size) << " monomials..." << endl;
    
    thread t[THREADS];

    uint32_t num_chunks = total_map_size/CHUNK_SIZE;

    for (uint32_t chunk = 0; chunk < 1; ++chunk)
    {
        for (uint32_t i = 0; i < THREADS; ++i) {
            t[i] = thread(sequential_generate, mask, i, CHUNK_SIZE * chunk, chunk);
        }

        for (uint32_t i = 0; i < THREADS; ++i) {
            t[i].join();
        }
        cout << "\tChunk (" << chunk + 1 << "/" << num_chunks << ") : " 
             << CHUNK_SIZE * sizeof(element_t) / 1024 / 1024 
             << " MBs written to disk." << endl;
    }
}

void thread_merge()
{

}


int main() 
{ 
    cout << "Polynomial:  " << polynomial_representation(POLY).str() << endl;
    cout << "Degree:      " << polynomial_degree << endl;
    cout << "Alpha:       " << ALPHA << endl;
    cout << "Beta:        " << BETA << endl;
#if THREADS > 1
    cout << "Threads:     " << THREADS << endl;
#ifdef IN_MEM_GENERATION
    cout << "Buckets:     " << LAYER_SIZE << endl;
#else
    cout << "Chunksize:   " << CHUNK_SIZE << endl;
    cout << "Chunks:      " << total_map_size/CHUNK_SIZE << endl;
#endif
#endif
    cout << "Seed:        " << SEED << endl;
    cout << "Mask:        " << hexmask_representation(mask).str() << endl;

#ifdef IN_MEM_GENERATION
    in_memory_search(mask);
#else
    thread_generate(mask);
#endif

    return 0;
} 
