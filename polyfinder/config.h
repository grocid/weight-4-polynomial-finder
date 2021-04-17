#pragma once

//This no longer needs to be prime :D
//TODO use this to distribute the thread local searches (and reduce the cmap size)
#define THREADS 8
#define BUCKETS (65536/THREADS)
#define COLL_BUCKETS (65536/THREADS)
#define BETA    1
#define ALPHA   1
//#define SEED    (time(NULL)+1)
#define SEED    (1618440951)
// You may want to set these to 1 when 2^x exponentiation is fast

// Drop the inverse mask from the data stored in the hash
// this gains imasklenb bytes per entry (up to max(2^masklen,total_map_size))
// at the expense of one 2^x exponentiation per generated binomial
#define CMAP_DROP_IMASK 1

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
#include <cstdint>
#include "utils.cpp"

static constexpr uint32_t polynomial_degree = get_degree(POLY);
static constexpr uint32_t log_mdrop = get_degree(BUCKETS*THREADS);
static constexpr uint64_t masklen = polynomial_degree/3-2;
static constexpr uint64_t imasklen = (polynomial_degree+1)-masklen;
static constexpr uint64_t masklenb = (masklen-log_mdrop+7)/8;
static constexpr uint64_t imasklenb = (imasklen+7)/8;
static constexpr uint128_t pmask = polynomial_degree == 127? ((uint128_t)-1) : (((uint128_t)1) << (polynomial_degree+1))-((uint128_t)1);
static constexpr uint64_t total_map_size = optimize_map_size((uint64_t)(((uint64_t)1 << (polynomial_degree/3 + BETA)) * ALPHA));
//Items go from 1 to total_map_size-1 the 0 would mean we have to take away one but we already did so with the optimization
static constexpr uint64_t exp_bit_len = get_degree(total_map_size-1) + 1;
static constexpr uint64_t exp_len = (get_degree(total_map_size-1) + 8) / 8;
static constexpr uint64_t coll_set_size = total_map_size/(THREADS*THREADS*BUCKETS);
static constexpr uint64_t base_coll_set_size = total_map_size/(THREADS*BUCKETS);
static constexpr uint64_t stage1_table_len = (1 << masklen) / THREADS;
