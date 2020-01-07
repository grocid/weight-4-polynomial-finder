#include <iostream>
#include <iomanip>
#include <sstream>
#include <thread>

#include <execution>
#include <algorithm>

#include <boost/functional/hash.hpp>
#include "parallel_hashmap/phmap.h"

#include "gf2_monomial.cpp"
#include "common.cpp"
#include "stringops.cpp"

#define THREADS 23
#define BETA    1
#define ALPHA   1.3
#define SEED    0*(time(NULL)+1)

// LILI-128
//#define POLY ((((uint128_t) 257) << 64) + (uint128_t)(299733420677929))
// Bluetooth E0
#define POLY ((((uint128_t) 0x2090000) << 64) + 0xa0048000000003)
// Bluetooth E0
//#define POLY ((((uint128_t) 0x212011231) << 64) + 0x2fea1a9dc028a229)

#define LAYER_SIZE (THREADS*THREADS)

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

inline uint32_t phi(uint128_t val)
{
    return val % THREADS;
}

void in_memory_generate(uint128_t mask, int thread)
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

                added++;
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

        px <<= 1;
        if ((px & f) != 0){
            px ^= POLY;
        }
        exponent++;
    }

    cout << "Thread " << thread  << " produced " << added << " in first step." << endl;
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

        //delete collision_layer[THREADS * i + thread];
    }
}




/* 
void sequential_generate()
{

    vector<tuple<uint128_t, uint64_t, uint64_t>> elements;

    uint128_t f = (((uint128_t) 1) << DEGREE);
    for(uint64_t exponent = 0; exponent < SIZE; ++exponent)
    {
        px <<= 1;
        if ((px & f) != 0){
            px ^= POLY;
        }
        elements.push_back ({px, exponent, 0});
    }

    sort(execution::par_unseq, input.begin(), input.end());
}*/

int main() 
{ 
    polynomial_degree = get_degree(POLY);
    uint128_t mask = random_mask(0, polynomial_degree/3, polynomial_degree, SEED);
    total_map_size = (uint64_t)(((uint64_t)1 << (polynomial_degree/3 + BETA)) * ALPHA);

    cout << "Polynomial:  " << polynomial_representation(POLY).str() << endl;
    cout << "Degree:      " << polynomial_degree << endl;
    cout << "Alpha:       " << ALPHA << endl;
    cout << "Beta:        " << BETA << endl;
    cout << "Threads:     " << THREADS << endl;
    cout << "Buckets:     " << LAYER_SIZE << endl;
    cout << "Seed:        " << SEED << endl;
    cout << "Mask:        " << hexmask_representation(mask).str() << endl;

    thread t[THREADS];

    cout << endl << "Generating 2^" << log2(total_map_size) << " monomials..." << endl;
    
    for (int i = 0; i < THREADS; ++i) {
        t[i] = thread(in_memory_generate, mask, i);
    }

    for (int i = 0; i < THREADS; ++i) {
        t[i].join();
    }

    cout << "Merging binomials..." << endl << endl;

    for (int i = 0; i < THREADS; ++i) {
        t[i] = thread(in_memory_merge, mask, i);
    }

    for (int i = 0; i < THREADS; ++i) {
        t[i].join();
    }

    return 0;
} 
