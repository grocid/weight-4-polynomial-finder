#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;
#ifdef __GNUC__
__extension__ typedef unsigned __int128 uint128_t;
#else
typedef unsigned __int128 uint128_t;
#endif

stringstream polynomial_representation(uint128_t poly) {
    stringstream polynomial;
    for(size_t i = 127; i >= 1; i--)
    {
        if ((poly >> i) & 0x1)
        {
            polynomial << "x^" << i << "+";
        }
    }
    polynomial << "1";

    return polynomial;
}

stringstream polynomial_representation(std::vector<uint64_t> exponents)
{
    sort(exponents.begin(), exponents.end());
    reverse(exponents.begin(), exponents.end());

    stringstream polynomial;
    for(size_t i = 0; i < exponents.size() - 1; ++i)
    {
        polynomial << "x^" << exponents[i] - exponents[exponents.size()-1] << "+";
    }
    polynomial << "1";

    return polynomial;
}

stringstream hexmask_representation(uint128_t mask)
{
    stringstream hexmask;
    hexmask << internal << setfill('0');
    hexmask << hex << setw(16) << (uint64_t)(mask >> 64);
    hexmask << hex << setw(16) << (uint64_t)(mask);
    return hexmask;
}

stringstream generate_fname(uint32_t thread, uint32_t step)
{
    stringstream fname;
    fname << thread << "_" << step << ".cached";
    return fname;
}
