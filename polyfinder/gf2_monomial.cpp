typedef unsigned __int128 uint128_t;


#ifdef __GNUC__
#define PUREFUN [[gnu::pure]]
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#else
#define PUREFUN
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

#if __GNUC__
PUREFUN inline uint32_t get_degree(uint128_t px)
{
    unsigned long long hi = px >> 64, lo = px;
    //NOTE: DO NOT REMOVE the != if you do this becomes a long on 32-bits and results get messed up!
    if (likely(hi != 0))
        return 127-__builtin_clzll(hi);
    else
        return 63-__builtin_clzll(lo|1);
}
#else
PUREFUN inline uint32_t get_degree(uint128_t px)
{
    uint32_t degree = 0;

    for (uint32_t i = 0; i < 128; ++i)
    {
        if (px & 0x1)
        {
            degree = i;
        }
        px >>= 1;
    }

    return degree;
}
#endif

PUREFUN inline uint128_t next_monomial(uint128_t px, uint128_t f, uint128_t poly)
{
    px <<= 1;
    //NOTE: DO NOT REMOVE the != if you do this becomes a long and results get messed up!
    if (likely((px & f) != 0))
        px ^= poly;
    return px;
}


/*
    GF(2) polynomial multiplication mod P(x)
*/
PUREFUN inline uint128_t gf2x_multiply(uint128_t a, uint128_t b, uint128_t f, uint128_t poly)
{
    uint128_t q = b, r = 0;

    for (int i = 0; i < 128; ++i)
    {
        if ((a & ((uint128_t)1 << i)))
            r ^= q;
        q = next_monomial(q,f,poly);
    }
    return r;
}

/*
    GF(2) polynomial exponentiation mod P(x)
*/
PUREFUN inline uint128_t gf2x_exp(uint128_t x, uint128_t n, uint128_t poly)
{
    uint128_t y = 1, f = (((uint128_t) 1) << get_degree(poly));
    for (;n != 0;n >>= 1)
    {
        if (n & 1) // n %2  == 1
            y = gf2x_multiply(x, y, f, poly);
        x = gf2x_multiply(x, x, f, poly);
    }
    return y;
}
