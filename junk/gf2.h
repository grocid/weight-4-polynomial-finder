#include <stdint.h>
#include <inttypes.h>

typedef unsigned __int128 uint128_t;

uint128_t multiply(uint128_t a, uint128_t b, uint128_t f, uint128_t poly, int degree);
uint128_t gf2_exp(uint128_t a, uint128_t b, uint128_t f, uint128_t poly, int degree);