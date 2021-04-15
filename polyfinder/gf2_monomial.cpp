#ifdef __GNUC__
#define PUREFUN [[gnu::pure]]
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#else
#define PUREFUN
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

#ifdef __GNUC__
PUREFUN constexpr inline uint32_t get_degree(uint128_t px)
{
    unsigned long long hi = px >> 64, lo = px;
    if (likely(hi != 0))
        return 127-__builtin_clzll(hi);
    else
        return 63-__builtin_clzll(lo|1);
}

PUREFUN inline uint32_t popcnt128(uint128_t px)
{
    unsigned long long hi = px >> 64, lo = px;
    return __builtin_popcountll(hi) + __builtin_popcountll(lo);
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

PUREFUN inline uint32_t popcnt128(uint128_t px)
{
    uint32_t rv = 0;

    for (uint32_t i = 0; i < 128; ++i)
    {
        rv += (px & 0x1);
        px >>= 1;
    }
    return rv;
}
#endif

// px*2
PUREFUN inline uint128_t next_monomial(uint128_t px, uint128_t poly, uint32_t f)
{
    px <<= 1;
    if ((px & (((uint128_t)1) << f)) != 0)
        px ^= poly;
    return px;
}


/*
    GF(2) polynomial multiplication mod P(x)
*/
PUREFUN inline uint128_t gf2x_multiply(uint128_t a, uint128_t b, uint128_t poly, uint32_t f)
{
    uint128_t r = 0;

    for (int i = 0; i < 128; ++i)
    {
        if (a&1)
            r ^= b;
        a>>=1;
        b = next_monomial(b,poly,f);
    }
    return r;
}

/*
    GF(2) polynomial exponentiation mod P(x)
*/
PUREFUN inline uint128_t gf2x_exp(uint128_t x, uint128_t n, uint128_t poly, uint32_t f)
{
    uint128_t y = 1;
    for (;n != 0;n >>= 1)
    {
        if (n & 1) // n %2  == 1
            y = gf2x_multiply(x, y, f, poly);
        x = gf2x_multiply(x, x, f, poly);
    }
    return y;
}
PUREFUN inline uint128_t gf2x_exp(uint128_t x, uint128_t n, uint128_t poly)
{
    return gf2x_exp(x, n, poly, get_degree(poly));
}
