#pragma once

//This no longer needs to be prime :D
#define THREADS 23
#define BETA    1
#define ALPHA   1
//#define SEED    (time(NULL)+1)
#define SEED    (1618440951)

//#define DEBUG_MESSAGES 1
#define IN_MEM_GENERATION 1

// LILI-128
//#define POLY ((((uint128_t) 257) << 64) + (uint128_t)(299733420677929))
// Bluetooth E0 x^1243366916+x^210123252+x^139501803+1
//#define POLY ((((uint128_t) 0x2090000) << 64) + 0xa0048000000003)
// Bluetooth E0
//#define POLY ((((uint128_t) 0x212011231) << 64) + 0x2fea1a9dc028a229)
#define POLY ((((uint128_t)1) << 80)|(0xa195150d15*2+1))

// Some constants do not touch after this line
typedef unsigned __int128 uint128_t;
#include "gf2_monomial.cpp"

static constexpr uint32_t polynomial_degree = get_degree(POLY);
static constexpr uint64_t masklen = polynomial_degree/3-2;
static constexpr uint64_t imasklen = (polynomial_degree+1)-masklen;
static constexpr uint64_t masklenb = (masklen+7)/8;
static constexpr uint64_t imasklenb = (imasklen+7)/8;
