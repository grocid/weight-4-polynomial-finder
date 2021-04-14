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

#include "gf2_monomial.cpp"
#include "common.cpp"
#include "stringops.cpp"

#define THREADS 23
#define BETA    1
#define ALPHA   1
//#define SEED    (time(NULL)+1)
#define SEED    (1618440951)

//#define DEBUG_MESSAGES 1
#define IN_MEM_GENERATION 1

// LILI-128
//#define POLY ((((uint128_t) 257) << 64) + (uint128_t)(299733420677929))
// Bluetooth E0 x^1243366916+x^210123252+x^139501803+1
//#define POLY ((((uint128_t) 0x2090000) << 64) + 0xa0048000000003)
// Bluetooth E0
//#define POLY ((((uint128_t) 0x212011231) << 64) + 0x2fea1a9dc028a229)
#define POLY ((((uint128_t)1) << 80)|(0xa195150d15*2+1))

#define LAYER_SIZE (THREADS*THREADS)

// Chunk size is the total number of elements contained in memory
// at the same time, which equates to how much RAM is available.
#define CHUNK_SIZE (1 << 25)

using namespace std;
typedef unsigned __int128 uint128_t;

namespace std {
    template<>
    struct hash<uint128_t>
    {
       typedef uint128_t argument_type;
       typedef size_t result_type;

       result_type operator () (const argument_type& x) const
       {
          return boost::hash<uint128_t>()(static_cast<uint128_t>(x));
       }
    };
}

#define map_t phmap::flat_hash_map
#define EXTRAARGS phmap::container_internal::hash_default_hash<K>, \
                  phmap::container_internal::hash_default_eq<K>, \
                  std::allocator<std::pair<const K, V>>, 32, phmap::NullMutex

template <class K, class V> using layer = map_t<K, V>;

map_t<uint128_t, tuple<uint64_t,  uint64_t>> collision_layer[LAYER_SIZE];

static uint64_t total_map_size;
static uint32_t polynomial_degree;
static uint128_t mask;

inline uint32_t phi(uint128_t val)
{
    return val % THREADS;
}

void in_memory_generate(uint128_t mask, uint32_t thread)
{
    uint64_t exponent = 0;
    uint32_t idx, added = 0;
    uint128_t py, px = 1;
    uint128_t f = (((uint128_t) 1) << get_degree(POLY));

    map_t<uint128_t, tuple<uint128_t, uint64_t>> collision_map;
    map_t<uint128_t, tuple<uint128_t, uint64_t>>::iterator it;

    for(uint64_t i = 0; i < total_map_size; ++i)
    {
        // every thread considers each monomial in the range
        // but it will only work with monomials that match
        // the phi condition.
        idx = phi(px & mask);

        if (idx == thread)
        {
            auto [it, result] = collision_map.try_emplace(
                px & mask, tuple{px, exponent}
            );

            if(!result)
            {
                py = get<0>(it->second);                
                auto [it2, result] = collision_layer[THREADS * thread + phi(px ^ py)].try_emplace(
                    py ^ px, tuple{exponent, get<1>(it->second)}
                );
#ifdef DEBUG_MESSAGES
                added++;
#endif
                // This is rather unlikely to happen, but if we don't check
                // then we might miss a multiple.
                if(!result)
                {
                    vector<uint64_t> exponents;
                    exponents.push_back(exponent);
                    exponents.push_back(get<1>(it->second));
                    exponents.push_back(get<0>(it2->second));
                    exponents.push_back(get<1>(it2->second));

                    cout << polynomial_representation(exponents).str() << endl;
                    return;
                }
            }
        }

        px = next_monomial(px, f, POLY);
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
                exponents.push_back(get<0>(it.second));
                exponents.push_back(get<1>(it.second));
                exponents.push_back(get<0>(it2->second));
                exponents.push_back(get<1>(it2->second));

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
    uint128_t f = (((uint128_t) 1) << polynomial_degree);

    vector<element_t> sequential;

    for(uint64_t i = 0; i < CHUNK_SIZE; ++i)
    {
        idx = phi(px & mask);

        if (idx == thread)
        {
            sequential.push_back({px, exponent, 0});
        }

        px = next_monomial(px, f, POLY);
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
    polynomial_degree = get_degree(POLY);
    mask = random_mask(0, polynomial_degree/3-2, polynomial_degree, SEED);
    total_map_size = (uint64_t)(((uint64_t)1 << (polynomial_degree/3 + BETA)) * ALPHA);

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
