#include <stdint.h>
#include <inttypes.h>

typedef unsigned __int128 uint128_t;

/*
    GF(2) polynomial multiplication mod P(x)
*/
uint128_t multiply(uint128_t a, uint128_t b, uint128_t f, uint128_t poly, int degree)
{
    uint128_t r = 0;
    uint128_t q = b;
    for(int i = 0; i < degree; i++)
    {
        if((a & ((uint128_t)1 << i)) != 0)
        {
            r ^= q;
        }
        q = (q << 1);
        if ((q & f) != 0)
        {
            q ^= poly;
        }
    }
    return r;
}

/*
    GF(2) polynomial exponentiation mod P(x)
*/
uint128_t gf2_exp(uint128_t x, uint128_t n, uint128_t f, uint128_t poly, int degree)
{
    if(n == 0)
    {
        return 1;
    }
    uint128_t y = 1;
    while(n > 1)
    {
        if(n % 2 == 0)
        {
            x = multiply(x, x, f, poly, degree);
            n = n / 2;
        }
        else
        {
            y = multiply(x, y, f, poly, degree);
            x = multiply(x, x, f, poly, degree);
            n = (n - 1) / 2;
        }
    }
    return multiply(x, y, f, poly, degree);
}
