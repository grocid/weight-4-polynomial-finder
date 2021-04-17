#include <cstdint>
#include <cstdio>
#include "utils.cpp"
#include <emmintrin.h>

#define POLY ((((uint128_t)257) << 64) | ((uint128_t)299733420677929))
static constexpr uint32_t polynomial_degree = get_degree(POLY);
static constexpr uint128_t pmask = polynomial_degree == 127? ((uint128_t)-1) : (((uint128_t)1) << (polynomial_degree+1))-((uint128_t)1);

void constexpr printi128(uint128_t r) {
    for (uint32_t i = 0; i < 128; i+=8) {
        printf("%02x",(unsigned char)(r>>(120-i)));
    }
}

void constexpr printi256(uint256_t r) {
    printi128(r.hi);
    printi128(r.lo);
}

void constexpr printm128(__m128i t) {
    printi128(fromvector(t));
}

#include "gf2_monomial.cpp"


int main(int argc, char const *argv[])
{
    __m128i a = tovector((((uint128_t)0x00000000000001e6) << 64) | ((uint128_t)0xc3fbd3bcc67fa5cb));

    __m128i b = tovector64(1);
    __m128i c = tovector64(2);
    
    uint128_t px = 1;
    uint128_t r;
    __m128i t;
    
    
    printi256(shl(u256c(1), 2*polynomial_degree));
    printf("\n");
    printi128(gf2x_divide(shl(u256c(1), 2*polynomial_degree), POLY));
    printf("\n");
    printf("%d\n",2*polynomial_degree);
    printm128(mpoly_inv);
    printf("\n");
    printm128(mpoly);
    printf("\n");
    
    for (uint32_t i = 0; i < 256; i++) {
        t = gf2x_mullo(b,c);
        r = fromvector(t);
        printf("%d ",i);
        for (uint32_t i = 0; i < 128; i+=8) {
            printf("%02x",(unsigned char)(r>>(120-i)));
        }
        printf("\n");
        t = gf2x_mul(b,c).lo;
        r = fromvector(t);
        printf("%d ",i);
        for (uint32_t i = 0; i < 128; i+=8) {
            printf("%02x",(unsigned char)(r>>(120-i)));
        }
        printf("\n");
        b = gf2x_reduce(gf2x_mul(b, c));
        px = next_monomial(px);
        r = fromvector(b);
        printf("%d ",i);
        for (uint32_t i = 0; i < 128; i+=8) {
            printf("%02x",(unsigned char)(r>>(120-i)));
        }
        printf("\n");
        printf("%d ",i);
        for (uint32_t i = 0; i < 128; i+=8) {
            printf("%02x",(unsigned char)(px>>(120-i)));
        }
        printf("\n");
    }

    return 0;
}
