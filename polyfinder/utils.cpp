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
