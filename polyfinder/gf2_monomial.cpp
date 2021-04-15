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
struct exp2vals{
  uint128_t arr[size];

  // 'constexpr' constructor:
  constexpr exp2vals():arr(){
    if (size > 1) arr[0] = 1;
    if (size > 2) arr[1] = 2;
    for(int i = 2; i < size; i++)
      arr[i] = gf2x_multiply(arr[i-1], arr[i-1]);
  }
};
static constexpr exp2vals arre2 = exp2vals<sizeof(uint128_t)*8+1>();

template<int size>
struct mul2vals{
  uint128_t arr[size][size];

  // 'constexpr' constructor:
  constexpr mul2vals():arr(){
    for(int i = 0; i < size; i++) {
        arr[i][0] = arre2.arr[i];
        for(int j = 1; j < size; j++)
            arr[i][j] = next_monomial(arr[i][j-1]);
    }
  }
};
static constexpr mul2vals arrm2 = mul2vals<sizeof(uint128_t)*8+1>();

__attribute__((always_inline)) PUREFUN static constexpr inline uint128_t gf_mul2(uint128_t a,uint32_t iv) {
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
