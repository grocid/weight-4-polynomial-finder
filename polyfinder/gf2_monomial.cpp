// This code contains a mix of optimizations along with code blatantly "inspired" on QMechanic's
// amazing work. Thanks man! You gave me the bits I was missing to make this go faster than GME
// https://github.com/QMechanic/GF2_128_lib/blob/master/gf2x_arithmetic.cpp

// px*2
__attribute__((always_inline)) PUREFUN static constexpr inline uint128_t next_monomial(uint128_t px)
{
    px <<= 1;
    if ((px & (((uint128_t)1) << polynomial_degree)) != 0)
        px ^= POLY;
    return px;
}


/*
    GF(2) POLYnomial multiplication mod P(x)
*/
PUREFUN static constexpr inline uint128_t gf2x_multiply(uint128_t a, uint128_t b)
{
    uint128_t r = 0;
    for (uint32_t i = 0; i < polynomial_degree+1; ++i)
    {
        if (a&(((uint128_t)1) << i))
            r ^= b;
        b = next_monomial(b);
    }
    return r & pmask;
}

/*
    GF(2) POLYnomial exponentiation mod P(x)
*/
PUREFUN static constexpr inline uint128_t gf2x_exp(uint128_t x, uint64_t n)
{
    uint128_t y = 1;
#pragma GCC unroll 64
    for (uint32_t i = 0; i < 64 && n; ++i,n/=2)
    {
        if (n & 1) // n %2  == 1
            y = gf2x_multiply(x, y);
        x = gf2x_multiply(x, x);
    }
    return y;
}

PUREFUN static constexpr inline uint128_t gf2x_divide(uint256_t a, uint128_t b)
{
    const uint32_t a_deg = get_degree(a);
    const uint32_t b_deg = get_degree(b);

    uint256_t t = u256c(b);
    uint128_t out = 0;
    t = shl(t, a_deg - b_deg);
    for (int32_t i = a_deg - b_deg; i >= 0; --i)
    {
        if (get_degree(a) == get_degree(t))
        {
            a = pxor(a,t);
            out |= ((uint128_t)1) << i;
        }
        t = shr(t, 1);
    }
    return out;
}


#if defined(__x86_64__) || defined(__i386__)
#include <emmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#include <wmmintrin.h>

typedef struct gfmm256 {
    __m128i lo;
    __m128i hi;
} gfmm256;

PUREFUN static constexpr inline __m128i tovector(uint128_t a)
{
    uint64_t hi = a >> 64, lo = a;
    return __extension__ (__m128i)(__v2du){ lo, hi };
}

PUREFUN static constexpr inline __m128i tovector64(uint64_t a)
{
    //Use little endian so shifts are efficient
    return __extension__ (__m128i)(__v2du){ a, 0 };
}


PUREFUN  static constexpr inline uint128_t fromvector(__m128i a)
{
    uint64_t hi = __extension__ ((__v2du) a)[1], lo = __extension__ ((__v2du) a)[0];
    return (((uint128_t)hi)<<64)|((uint128_t)lo);
}

// karatsuba algorithm
// a*b, 128bit input, 256 bit output
PUREFUN  static inline gfmm256 gf2x_mul(__m128i a, __m128i b)
{
    gfmm256 out = {};
    __m128i z0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i z2 = _mm_clmulepi64_si128(a, b, 0x11);
    __m128i t1 = _mm_xor_si128(_mm_unpackhi_epi64(a,b),_mm_unpacklo_epi64(a,b));
    __m128i z1 = _mm_xor_si128(_mm_clmulepi64_si128(t1, t1, 0x10),_mm_xor_si128(z0, z2));
    //Order of unpack is first low bits then high bits
    // Return: z0lo, z0hi ^ z1lo, z1hi ^ z2lo, z2hi
    __m128i t2 = _mm_xor_si128(_mm_unpackhi_epi64(z0,z1),_mm_unpacklo_epi64(z1,z2));
    out.lo = _mm_unpacklo_epi64(z0,t2);
    out.hi = _mm_unpackhi_epi64(t2,z2);
    return out;
}

// Same but dropping the higher 64-bits use schoolbook instead
// TODO: maybe there is a smart trick to get this with only two multiplications
PUREFUN  static inline __m128i gf2x_mullo(__m128i a, __m128i b)
{
    return _mm_xor_si128(_mm_bslli_si128(_mm_xor_si128(_mm_clmulepi64_si128(a, b, 0x01),_mm_clmulepi64_si128(a, b, 0x10)),8),_mm_clmulepi64_si128(a, b, 0x00));
}

// a², 128 bit input, 256 bit output
PUREFUN static inline gfmm256 gf2x_square(__m128i a)
{
    gfmm256 out = {};
    out.lo = _mm_clmulepi64_si128(a, a, 0x00);
    out.hi = _mm_clmulepi64_si128(a, a, 0x11);
    return out;
}

//Shift right by polynomial degree
PUREFUN static inline __m128i pdgshr(gfmm256 in) {
    constexpr uint32_t div = (polynomial_degree / 64) * 8;
    constexpr uint32_t div2 = polynomial_degree / 8;
    constexpr uint32_t rem = polynomial_degree % 64;
    constexpr uint32_t rem2 = 64-rem;
    if constexpr (polynomial_degree >= 256) {
        return _mm_setzero_si128();
    } else if constexpr (polynomial_degree == 128) {
        return in.hi;
    } else if constexpr (polynomial_degree % 8 == 0) {
        //Combine hi and lo bytewise
        if constexpr (polynomial_degree > 128)
            return _mm_bsrli_si128(in.hi,div2);
        else
            return _mm_alignr_epi8(in.hi,in.lo,div2);
    } else {
        //Combine hi and lo bitwise slow!
        if constexpr (polynomial_degree > 192)
            return _mm_srli_epi64(_mm_bsrli_si128(in.hi, 8),rem);
        else if constexpr (polynomial_degree > 128)
            return _mm_or_si128(_mm_srli_epi64(in.hi,rem),_mm_slli_epi64(_mm_bsrli_si128(in.hi,8),rem2));
        else if constexpr (polynomial_degree > 64)
            return _mm_or_si128(_mm_srli_epi64(_mm_alignr_epi8(in.hi,in.lo,div),rem),_mm_slli_epi64(in.hi,rem2));
        else
            return _mm_or_si128(_mm_srli_epi64(in.lo,rem),_mm_slli_epi64(_mm_alignr_epi8(in.hi,in.lo,div+8),rem2));
    }
}


static constexpr __m128i mpoly_inv = tovector(gf2x_divide(shl(u256c(1), 2*polynomial_degree), POLY));
static constexpr __m128i mpoly = tovector(POLY);
// This calculates in(x) mod P(x) with just multiplications and additions
// P(x) is assumed to be constant and that's why this can optimize it
// TODO: This maybe can be made faster if we make multiplications which are shift-aware
PUREFUN static inline __m128i gf2x_reduce(gfmm256 in)
{
    return _mm_xor_si128(in.lo,gf2x_mullo(pdgshr(gf2x_mul(pdgshr(in), mpoly_inv)), mpoly));
}

PUREFUN static inline uint128_t gf2x_multiply_fast(uint128_t a, uint128_t b)
{
    return fromvector(gf2x_reduce(gf2x_mul(tovector(b), tovector(a))));
}

PUREFUN static inline uint128_t gf2x_square_fast(uint128_t a)
{
    return fromvector(gf2x_reduce(gf2x_square(tovector(a))));
}


template<int size>
struct exp2vals{
  __m128i arr[size];
  constexpr exp2vals():arr(){
    uint128_t last = 2;
    if (size > 0) arr[0] = tovector(last);
    for(int i = 1; i < size; i++) {
        last = gf2x_multiply(last, last);
        arr[i] = tovector(last);
    }
  }
};

static constexpr exp2vals arre2 = exp2vals<64>();

PUREFUN static inline uint128_t gf_exp2(uint64_t n)
{
    __m128i y = tovector64(1);
    #pragma GCC unroll 64
    for (uint32_t i = 0; i < 64; ++i)
    {
        if (unlikely(n==0))
            return fromvector(y);
        else if (unlikely(n % 2 == 1))
            y = gf2x_reduce(gf2x_mul(y,arre2.arr[i]));
        n>>=1;
    }
    return fromvector(y);
}


#else
PUREFUN static inline uint128_t gf2x_multiply_fast(uint128_t a, uint128_t b)
{
    uint128_t r = 0;

#pragma GCC unroll 128
    for (uint32_t i = 0; i < polynomial_degree+1; ++i)
    {
        if (a&(((uint128_t)1) << i))
            r ^= b;
        b = next_monomial(b);
    }
    return r & pmask;
}
PUREFUN static inline uint128_t gf2x_square_fast(uint128_t a)
{
    return gf2x_multiply_fast(a,a);
}

template<int size1,int size>
struct mul2vals{
  uint128_t arr[size1][size];

  // 'constexpr' constructor:
  constexpr mul2vals():arr(){
    for(int i = 0; i < size1; i++) {
        if (i == 0) arr[i][0] = 1;
        else if (i == 1) arr[i][0] = 2;
        else arr[i][0] = gf2x_multiply(arr[i-1][0], arr[i-1][0]);
        for(int j = 1; j < size; j++)
            arr[i][j] = next_monomial(arr[i][j-1]);
    }
  }
};

__attribute__((always_inline)) PUREFUN constexpr static inline uint128_t gf_mul2(uint128_t a,uint32_t iv) {
    static const mul2vals arrm2 = mul2vals<64,polynomial_degree+1>();
    uint128_t r = 0;

#pragma GCC unroll 128
    for (uint32_t i = 0; i < polynomial_degree+1; ++i)
    {
        if (unlikely(a == 0))
            return r;
        else if (a&1)
            r ^= arrm2.arr[iv][i];
        a>>=1;
    }
    return r;
}


__attribute__ ((noinline)) static uint128_t gf_exp2(uint64_t n)
{
    uint128_t y = 1;
    #pragma GCC unroll 64
    for (uint32_t i = 0; i < 64; ++i)
    {
        if (unlikely(n==0))
            return y;
        else if (unlikely(n % 2 == 1))
            y = gf_mul2(y,i);
        n>>=1;
    }
    return y;
}
#endif
