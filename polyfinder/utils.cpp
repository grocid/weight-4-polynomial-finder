#ifdef __GNUC__
#define PUREFUN [[gnu::pure]]
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
__extension__ typedef unsigned __int128 uint128_t;
#else
#define PUREFUN
#define likely(x)       (x)
#define unlikely(x)     (x)
typedef unsigned __int128 uint128_t;
#endif

typedef struct uint256_t {
    uint128_t hi;
    uint128_t lo;
} uint256_t;

#ifdef __GNUC__
PUREFUN constexpr inline uint32_t get_degree(uint128_t px)
{
    unsigned long long hi = px >> 64, lo = px;
    if (likely(hi != 0))
        return 127-__builtin_clzll(hi);
    else
        return 63-__builtin_clzll(lo|1);
}

PUREFUN constexpr inline int popcnt128(uint128_t px)
{
    unsigned long long hi = px >> 64, lo = px;
    return __builtin_popcountll(hi) + __builtin_popcountll(lo);
}
PUREFUN constexpr inline int popcnt64(uint64_t px)
{
    return __builtin_popcountll((unsigned long long)px);
}
#else
PUREFUN constexpr inline uint32_t get_degree(uint128_t px)
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

PUREFUN constexpr inline int popcnt128(uint128_t px)
{
    int rv = 0;

    for (uint32_t i = 0; i < 128; ++i)
    {
        rv += (px & 0x1);
        px >>= 1;
    }
    return rv;
}

PUREFUN constexpr inline int popcnt64(uint64_t px)
{
    int rv = 0;

    for (uint32_t i = 0; i < 64; ++i)
    {
        rv += (px & 0x1);
        px >>= 1;
    }
    return rv;
}
#endif

//For now this is a very trivial optimization, if the map_size is going to
//be a power of two drop one element. This is because we use element 0
//to indicate that the element is not present on the allocated memory
PUREFUN constexpr inline uint64_t optimize_map_size (uint64_t map_size)
{
    if (map_size > (1<<20) && popcnt64(map_size) == 1)
        return map_size-1;
    else
        return map_size;
}


PUREFUN constexpr inline uint32_t get_degree(uint256_t px)
{
    if (likely(px.hi != 0))
        return 128+get_degree(px.hi);
    else
        return get_degree(px.lo);
}

PUREFUN constexpr inline uint256_t u256c (uint128_t lo)
{
    uint256_t rv = {}; rv.hi = 0; rv.lo = lo; return rv;
}

PUREFUN constexpr inline uint256_t u256c (uint128_t hi, uint128_t lo)
{
    uint256_t rv = {}; rv.hi = hi; rv.lo = lo; return rv;
}

PUREFUN constexpr inline uint256_t pxor (uint256_t px, uint256_t py)
{
    uint256_t rv = {}; rv.hi = px.hi ^ py.hi; rv.lo = px.lo ^ py.lo; return rv;
}

PUREFUN constexpr inline uint256_t shr (uint256_t px, uint64_t pos)
{
    uint256_t rv = {};
    if (pos == 0) return px;
    else if (pos >= 256) { rv.hi = 0; rv.lo = 0; return rv; }
    else if (pos == 128) { rv.hi = 0; rv.lo = px.hi; return rv; }
    else if (pos > 128) { rv.hi = 0; rv.lo =  px.hi >> (pos-128); return rv; }
    else { rv.hi = px.hi >> pos; rv.lo = (px.hi << (128 - pos)) | (px.lo >> pos); return rv; }
}

PUREFUN constexpr inline uint256_t shl (uint256_t px, uint64_t pos)
{
    uint256_t rv = {};
    if (pos == 0) return px;
    else if (pos >= 256) { rv.hi = 0; rv.lo = 0; return rv; }
    else if (pos == 128) { rv.hi = px.lo; rv.lo = 0; return rv; }
    else if (pos > 128) { rv.hi =  px.lo << (pos-128); rv.lo = 0; return rv; }
    else { rv.hi = px.hi << pos | (px.lo >> (128 - pos)); rv.lo = px.lo << pos; return rv; }
}
