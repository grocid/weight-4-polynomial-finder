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
    int maskbound = 0; // Shift for the high level bits mask
    unsigned long long masklo=0, maskhi=0;
    const char * maskinttype = masklen > 64 ? "uint128_t" : "uint64_t";
    while (maskbound < masklen && maskbits[maskbound] < 64) maskbound ++;
    for (int i = 0; i < maskbound; i++) {
        masklo |= ((unsigned long long)1) <<  maskbits[i];
    }
    for (int i = maskbound; i < masklen; i++) {
        maskhi |= ((unsigned long long)1) <<  maskbits[i];
    }
    cout << "typedef "<< maskinttype <<" " << name << "_int_t;" << endl;
    cout << "#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__)) && defined(__BMI2__)" <<endl;
    cout << "PUREFUN constexpr inline " << maskinttype << " get_" << name << "_bits(uint128_t px) { return ";
    if (maskbound < masklen)
        cout << "(((" << maskinttype << ")__builtin_ia32_pext_di((unsigned long long)(px>>64),(unsigned long long)" << maskhi << "U)) << " << maskbound << ") |";
    cout << "((" << maskinttype << ")__builtin_ia32_pext_di((unsigned long long)px,(unsigned long long)" << masklo << "U))";
    cout << "; }" << endl;
    cout << "#else" <<endl;
    cout << "PUREFUN constexpr inline " << maskinttype << " get_" << name << "_bits(uint128_t px) { return (" << maskinttype << ") 0 ";
    for (int i = 0; i < masklen; i++) {
        cout << "| ((((" << maskinttype << ")(px >> " << maskbits[i]-i << ")) & (((" << maskinttype << ")1) << "<< i <<"))) ";
    }
    cout << "; }" << endl;
    cout << "#endif" <<endl;
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
    cout << "static constexpr uint128_t mask = (((uint128_t)"<< ((uint64_t)(mask >> 64)) << "U) << 64)|((uint128_t)" << ((uint64_t)mask) << "U);"<< endl;
    cout << "static constexpr uint64_t maskseed = "<< maskseed << ";"<< endl;
    gen_mask_fun("mask",masklen,maskbits);
    gen_mask_fun("imask",imasklen,imaskbits);
    return 0;
} 
