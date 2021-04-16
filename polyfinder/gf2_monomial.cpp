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
PUREFUN static constexpr inline uint128_t gf2x_multiply(uint128_t b, uint128_t a)
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


template<int size>
struct mul2vals{
  uint128_t arr[size][size];

  // 'constexpr' constructor:
  constexpr mul2vals():arr(){
    for(int i = 0; i < size; i++) {
        if (i == 0) arr[i][0] = 1;
        else if (i == 1) arr[i][0] = 2;
        else arr[i][0] = gf2x_multiply(arr[i-1][0], arr[i-1][0]);
        for(int j = 1; j < size; j++)
            arr[i][j] = next_monomial(arr[i][j-1]);
    }
  }
};

__attribute__((always_inline)) PUREFUN static inline uint128_t gf_mul2(uint128_t a,uint32_t iv) {
    static const mul2vals arrm2 = mul2vals<sizeof(uint128_t)*8+1>();
    uint128_t r = 0;

#pragma GCC unroll 128
    for (uint32_t i = 0; i < polynomial_degree+1; ++i)
    {
        if (a&(((uint128_t)1) << i))
            r ^= arrm2.arr[iv][i];
    }
    return r;
}


__attribute__ ((noinline)) static uint128_t gf_exp2(uint64_t n)
{
    if (unlikely(n==0))
        return 1;
    uint128_t y = 1;
    #pragma GCC unroll 64
    for (uint32_t i = 0; i < 64; ++i)
    {
        if (n & 1) // n %2  == 1
            y = gf_mul2(y,i+1);
        if (unlikely((n>>=1) == 0))
            return y;
    }
    return y;
}
