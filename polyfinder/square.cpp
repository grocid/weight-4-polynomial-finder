#include <iostream>
#include <atomic>
#include <thread>
#include <new>
#include <fstream>

#include <execution>
#include <algorithm>

#include "parallel_hashmap/phmap.h"

#include "config.h"

#include "gf2_monomial.cpp"
#include "stringops.cpp"
#include "masks.hpp"

#define LAYER_SIZE (THREADS*COLL_BUCKETS)

// Chunk size is the total number of elements contained in memory
// at the same time, which equates to how much RAM is available.
#define CHUNK_SIZE (1 << 25)

using namespace std;

#ifdef __GNUC__
#define PACKSTRUCT __attribute__((__packed__))
#else
#define PACKSTRUCT
#endif

typedef struct PACKSTRUCT exp_t {
    uint8_t exponent[exp_len];
} exp_t;

typedef struct PACKSTRUCT imask_t {
    uint8_t imask[imasklenb];
} imask_t;

typedef struct PACKSTRUCT cmap_poly {
#if !CMAP_DROP_IMASK
    imask_t imask;
#endif
    exp_t exponent;
} cmap_poly;

//TODO: use a bitpacked representation to save a bit more of space
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

PUREFUN constexpr inline imask_t pack_imask(imask_int_t imask) {
   imask_t rv = {};
    for (uint i = 0; i < imasklenb; i++)
        rv.imask[i] = imask >> (8*i);
    return rv;
}


PUREFUN constexpr inline imask_int_t unpack_imask(const imask_t &imask) {
    imask_int_t rv = 0;
    for (uint i = 0; i < imasklenb; i++)
        rv |= ((imask_int_t)imask.imask[i]) <<  (8*i);
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
    template<>
    struct hash<clay_poly>
    {
        typedef clay_poly argument_type;
        typedef size_t result_type;
        PUREFUN inline result_type operator () (const argument_type& x) const
        {
            imask_int_t t = get_imask_bits(gf_exp2(unpack_exp(x.e1))^gf_exp2(unpack_exp(x.e2)));
            if constexpr (sizeof(result_type) >= imasklen) {
                return t;
            } else {
                result_type rv = t;
                for (uint i = 1; i < (masklen + sizeof(result_type) -1)  / sizeof(result_type); i++) {
                    rv ^= t >> (i*(sizeof(result_type)*8));
                }
                return rv;
            }
        }
    };
    template<>
    struct hash<cmap_poly>
    {
        typedef cmap_poly argument_type;
        typedef size_t result_type;
        PUREFUN inline result_type operator () (const argument_type& x) const
        {
            mask_int_t t = get_mask_bits(gf_exp2(unpack_exp(x.exponent)));
            if constexpr (sizeof(result_type) >= masklen) {
                return t;
            } else {
                result_type rv = t;
                for (uint i = 1; i < (masklen + sizeof(result_type) -1)  / sizeof(result_type); i++) {
                    rv ^= t >> (i*(sizeof(result_type)*8));
                }
                return rv;
            }
        }
    };
}

PUREFUN constexpr inline bool operator==(const imask_t& lhs, const imask_t& rhs)
{
    return !memcmp(lhs.imask, rhs.imask, imasklenb);
}

PUREFUN inline bool operator==(const clay_poly& x, const clay_poly& y)
{
    return (gf_exp2(unpack_exp(x.e1))^gf_exp2(unpack_exp(x.e2))) == (gf_exp2(unpack_exp(y.e1))^gf_exp2(unpack_exp(y.e2)));
}

PUREFUN inline bool operator==(const cmap_poly& x, const cmap_poly& y)
{
    return (gf_exp2(unpack_exp(x.exponent))& mask)==((gf_exp2(unpack_exp(y.exponent)))& mask);
}

#define map_t phmap::flat_hash_map
#define set_t phmap::flat_hash_set
#define EXTRAARGS phmap::container_internal::hash_default_hash<K>, \
                  phmap::container_internal::hash_default_eq<K>, \
                  std::allocator<std::pair<const K, V>>, 32, phmap::NullMutex

template <class K, class V> using layer = map_t<K, V>;

//TODO: This is heavily limited by things like bss size, there are probably better ways to go around this
//Maybe being smart in using mmap can help with this issue
vector <clay_poly> collision_layer[THREADS][COLL_BUCKETS];

static inline void in_memory_generate(int threadid, uint64_t task)
{
    uint64_t exponent = 0;
    uint32_t idx;
#ifdef DEBUG_MESSAGES
    uint32_t added = 0;
    uint32_t candidated = 0;
#endif
    uint128_t px = 1;
    //When total_map_size is significantly larger than 2^masksize using an array is significantly faster and memory efficient
    //Here we use the trick that exponent will be 0 when the element has not been placed
    //TODO: we can save even more bits if we use bit sizes instead of byte sizes
    //TODO: I suspect the major slowdown here comes from TLB misses. We could use hugepages to see an improvement at the expense of physical memory.
    cmap_poly *collision_map = new cmap_poly[stage1_table_len];
    
    for(uint64_t i = 0; i < total_map_size; ++i)
    {
        // every thread considers each monomial in the range
        // but it will only work with monomials that match
        // the phi condition.
        mask_int_t mpx = get_mask_bits(px);
        idx = mpx % STAGE1_TASKS;
        //TODO: we may be able to achieve greater polynomial exponentiation speeds by using arrays of polynomials instead of single ones
        if (unlikely(idx == task))
        {
            mpx /= STAGE1_TASKS; //Reduce mpx to the real number of elements
#ifdef DEBUG_MESSAGES
            //Ensure exponentiation code is working as it should
            if (gf_exp2(exponent) != px)
                cerr << "Exponentiation error" << exponent << " " << hexmask_representation(gf_exp2(exponent)).str() << " " << hexmask_representation(px).str() << endl;
            candidated++;
#endif
#if !CMAP_DROP_IMASK
            imask_t imbits = pack_imask(get_imask_bits(px));
#endif
            uint64_t exponent2 = unpack_exp(collision_map[mpx].exponent);
            //If element with that mask not present already
            if (unlikely(exponent2 == 0)) {
                collision_map[mpx].exponent = pack_exp(exponent+1); //So 0 is kept correctly
#if !CMAP_DROP_IMASK
                collision_map[mpx].imask = imbits;
#endif
            } else {
                exponent2 -= 1; // Normalize to the right value
#if !CMAP_DROP_IMASK
                imask_t imaskxor = xor_imask(collision_map[mpx].imask, imbits);
                imask_int_t py = unpack_imask(imaskxor);
#else
                imask_int_t py = get_imask_bits(px^gf_exp2(exponent2));
#endif
#ifdef DEBUG_MESSAGES
                added++;
#endif
                //Since we know the size these will have, maybe use an mmap?
                collision_layer[threadid][py % COLL_BUCKETS].emplace_back(
                    clay_poly{pack_exp(exponent2), pack_exp(exponent)}
                );
            }
        }
        px = next_monomial(px);
        exponent++;
    }
    delete [] collision_map;
#ifdef DEBUG_MESSAGES
    cout << "\tThread " << thread  << " had " << candidated << " and produced " << added << " in first step." << endl;
#endif
}

static atomic_bool found_poly (false);
static inline void in_memory_merge(int threadid, uint64_t bucket)
{
    UNUSED(threadid);
    map_t <imask_t,clay_poly> base;
    size_t nelems = 0;
    for(int i = 1; i < THREADS; ++i) {
        nelems += collision_layer[i][bucket].size();
    }
    //By preallocating we ensure we avoid hashes
    base.reserve(nelems);
    for(int i = 1; i < THREADS; ++i) {
        for (auto& it: collision_layer[i][bucket])
        {
            imask_t imask = pack_imask(get_imask_bits(gf_exp2(unpack_exp(it.e1))^gf_exp2(unpack_exp(it.e2))));
            auto [it2, result] = base.try_emplace(
                imask, clay_poly {it}
            );

            if (unlikely(!result))
            {
                vector<uint64_t> exponents;
                exponents.push_back(unpack_exp(it.e1));
                exponents.push_back(unpack_exp(it.e2));
                exponents.push_back(unpack_exp(it2->second.e1));
                exponents.push_back(unpack_exp(it2->second.e2));
                cout << polynomial_representation(exponents).str() << endl;
                found_poly.store(true,memory_order_relaxed);
            }
        }
        collision_layer[i][bucket].clear();
    }
}


static atomic_uint_least64_t last_task (0);
void in_memory_generate_thread(int threadid) {
    uint64_t task = 0;
    //HACK: This will only work as long as STAGE1_TASKS + threads < 2^64 -1
    while ((task=last_task.fetch_add(1,memory_order_relaxed)) < STAGE1_TASKS) {
            in_memory_generate(threadid,task);
    }
}


static atomic_uint_least64_t last_bucket (0);
void in_memory_merge_thread(int threadid) {
    uint64_t bucket = 0;
    //HACK: This will only work as long as COLL_BUCKETS + threads < 2^64 -1
    while ((bucket=last_bucket.fetch_add(1,memory_order_relaxed)) < COLL_BUCKETS) {
            in_memory_merge(threadid,bucket);
    }
}


void in_memory_search (void)
{
    cout << endl << "Running in-memory (square algorithm) search..." << endl;
    cout << "[1/2]\t Generating 2^" << log2(total_map_size) << " monomials..." << endl;

    for (int i = 0; i < THREADS; ++i) {
        for (uint32_t k = 0; k < COLL_BUCKETS; ++k) {
            collision_layer[i][k].reserve(coll_set_size);
        }
    }
    last_task.store(0,memory_order_relaxed);
    thread t[THREADS];
    for (int i = 0; i < THREADS; ++i) {
        t[i] = thread(in_memory_generate_thread, i);
    }

    for (int i = 0; i < THREADS; ++i) {
        t[i].join();
    }

    //Free any unallocated space
    for (int i = 0; i < THREADS; ++i) {
        for (uint32_t k = 0; k < COLL_BUCKETS; ++k) {
            collision_layer[i][k].reserve(collision_layer[i][k].size());
        }
    }

    cout << "[2/2]\tMerging binomials..." << endl << endl;

    last_bucket.store(0,memory_order_relaxed);
    for (uint32_t i = 0; i < THREADS; ++i) {
        t[i] = thread(in_memory_merge_thread, i);
    }

    for (uint32_t i = 0; i < THREADS; ++i) {
        t[i].join();
    }
}

#ifndef IN_MEM_GENERATION
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
#endif

int main() 
{ 
    cout << "Polynomial:     " << polynomial_representation(POLY).str() << endl;
    cout << "Degree:         " << polynomial_degree << endl;
    cout << "Alpha:          " << ALPHA << endl;
    cout << "Beta:           " << BETA << endl;
    cout << "Stage1 tasks:   " << STAGE1_TASKS << endl;
    cout << "Stage2 buckets: " << COLL_BUCKETS << endl;
#if THREADS > 1
    cout << "Threads:        " << THREADS << endl;
#ifdef IN_MEM_GENERATION
    cout << "Buckets:        " << LAYER_SIZE << endl;
#else
    cout << "Chunksize:      " << CHUNK_SIZE << endl;
    cout << "Chunks:         " << total_map_size/CHUNK_SIZE << endl;
#endif
#endif
    cout << "Seed:           " << SEED << endl;
    cout << "Mask:           " << hexmask_representation(mask).str() << endl;
    cout << "Mask bits:      " << masklen << endl;
    cout << "~Mask bits:     " << imasklen << endl;
    cout << "Exp bit-size:   " << exp_bit_len << endl;
    cout << "Exp size:       " << sizeof(exp_t) << endl;
    cout << "Cmap size:      " << sizeof(cmap_poly) << endl;
    cout << "Binomial size:  " << sizeof(clay_poly) << endl;
    cout << "Clayer size:    " << sizeof(pair<imask_t, clay_poly>) << endl;
    cout << "Generating 2^" << log2(total_map_size) << " monomials..." << endl;
    found_poly.store(false,memory_order_relaxed);
#ifdef IN_MEM_GENERATION
    in_memory_search();
#else
    thread_generate(mask);
#endif
    if(found_poly.load(memory_order_relaxed))
        return 0;
    else
        return 1;
}
