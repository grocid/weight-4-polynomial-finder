typedef unsigned __int128 uint128_t;

inline uint32_t get_degree(uint128_t px)
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

/*
    GF(2) polynomial multiplication mod P(x)
*/
inline uint128_t gf2x_multiply(uint128_t a, uint128_t b, uint128_t f, uint128_t poly)
{
    uint128_t q = b, r = 0;

    for (int i = 0; i < 128; ++i)
    {
        if ((a & ((uint128_t)1 << i)) != 0)
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
inline uint128_t gf2x_exp(uint128_t x, uint128_t n, uint128_t poly)
{
    if (n == 0)
    {
        return 1;
    }

    uint128_t y = 1, f = (((uint128_t) 1) << get_degree(poly));
    while (n > 1)
    {
        if (n % 2 == 0)
        {
            x = gf2x_multiply(x, x, f, poly);
            n = n / 2;
        }
        else
        {
            y = gf2x_multiply(x, y, f, poly);
            x = gf2x_multiply(x, x, f, poly);
            n = (n - 1) / 2;
        }
    }

    return gf2x_multiply(x, y, f, poly);
}