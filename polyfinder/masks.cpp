#include <iostream>

#include "config.h"
using namespace std;

uint128_t random_mask(uint128_t mask, int weight, int degree, unsigned int seed) 
{
    uint64_t index = 0;
    uint128_t new_mask = 0;
    srand(seed);

    for(int i = 0; i < weight; ++i)
    {
        int t = 0;
        do
        {
            index = rand() % degree;
            t = ((new_mask >> index) & (uint128_t)0x1);
            t |= ((mask >> index) & (uint128_t)0x1);
        }
        while (t);

        new_mask = new_mask | ((uint128_t)1 << index);
    }

    return (uint128_t)new_mask;
}


//Find the relative positions of the mask bits
inline void gen_mask_bits(uint128_t mask, uint maskbits[128], uint32_t f, bool match) {
    uint masklen = 0;
    for (uint i = 0; i < f+1; i++) {
        if((mask & 1) == match) {
            maskbits[masklen++]=i;
        }
        mask >>= 1;
    }
}

inline void gen_mask_fun(const char *name, const int masklen, const uint maskbits[128]) {
    cout << "PUREFUN constexpr inline uint128_t " << name << "(uint128_t px) { return (uint128_t) 0 ";
    for (int i = 0; i < masklen; i++) {
        cout << "| (((px >> " << maskbits[i]-i << ") & (((uint128_t)1) << "<< i <<"))) ";
    }
    cout << "; }" << endl;
}

int main() 
{ 
    uint64_t maskseed = SEED;
    uint128_t mask = random_mask(0, masklen, polynomial_degree, maskseed);
    uint maskbits[masklen];
    uint imaskbits[imasklen];
    gen_mask_bits(mask, maskbits,polynomial_degree,1);
    gen_mask_bits(mask, imaskbits,polynomial_degree,0);
    cout << "#pragma once" << endl;
    cout << "static constexpr uint128_t mask = (((uint128_t)"<< ((uint64_t)(mask >> 64)) << ") << 64)|((uint128_t)" << ((uint64_t)mask) << ");"<< endl;
    cout << "static constexpr uint64_t maskseed = "<< maskseed << ";"<< endl;
    gen_mask_fun("get_mask_bits",masklen,maskbits);
    gen_mask_fun("get_imask_bits",imasklen,imaskbits);
    return 0;
} 
