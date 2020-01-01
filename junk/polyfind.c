#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <inttypes.h>

#include "gf2.h"

typedef unsigned __int128 uint128_t;

typedef struct {
  uint128_t key;
  uint64_t exp1, exp2;
} __attribute__((packed)) entry;

#define CLAMP 0
#define DEGREE (89-CLAMP)
#define THREADS 14

#define SIZE (uint64_t)(((uint64_t)1 << (DEGREE/3 + 1)) * 1.5)
#define LINES_PER_THREAD (SIZE/THREADS)
#define POLY ((((uint128_t) 257) << 64) + (uint128_t)(299733420677929))
//#define POLY ((((uint128_t) 0x212011231) << 64) + 0x2fea1a9dc028a229)
//#define POLY ((((uint128_t) 0x2090000) << 64) + 0xa0048000000003)

static uint128_t mask, clamp_mask;
static uint128_t f = (((uint128_t) 1) << (DEGREE+CLAMP));
static entry* table;

static pthread_t pids[THREADS];
static uint64_t args[THREADS];

int sorting_function(const void *a, const void *b) 
{ 
    entry *ia = (entry *)a;
    entry *ib = (entry *)b;

    if ((ia->key & mask) > (ib->key & mask)){
        return 1;
    } else if ((ia->key & mask)  < (ib->key & mask) ){
        return -1;
    } else {
        return 0;
    }
}

int sorting_function2(const void *a, const void *b) 
{ 
    entry *ia = (entry *)a;
    entry *ib = (entry *)b;

    if (ia->key > ib->key){
        return 1;
    } else if (ia->key < ib->key){
        return -1;
    } else {
        return 0;
    }
}

void print_binary(uint128_t mask) {
    for(int i = 127; i >= 0; i--) {
        printf("%d", (int)((mask >> i) & 0x1));
    }
    printf("\n");
}


/*
    generates a random mask of weight degree/3 + 2
*/
uint128_t random_mask(uint128_t mask, int weight) 
{
    uint64_t i1 = 0;
    uint128_t new_mask = 0;//(1 << (DEGREE/3)) - 1;
    srand(time(NULL)+1);
    int t;
    for(int i = 0; i < weight; ++i)
    {
        t = 0;
        do
        {
            i1 = rand() % DEGREE;
            t = ((new_mask >> i1) & (uint128_t)0x1);
            t |= ((mask >> i1) & (uint128_t)0x1);
        }
        while (t);

        new_mask = new_mask | ((uint128_t)1 << i1);
    }

    return (uint128_t)new_mask;
}

static void generate(uint64_t *line)
{
    uint64_t i = *line, j = i;
    uint128_t px;
    px = gf2_exp(2, i, f, POLY, DEGREE);
    while(j < *line + LINES_PER_THREAD)
    {
        if ((px & clamp_mask) == 0)
        {
            table[j].key = px;
            table[j].exp1 = i;
            table[j].exp2 = 0;
            j++;
        }

        i++;

        px <<= 1;

        if ((px & f) != 0){
            px ^= POLY;
        }
    }
}

static void match(uint128_t *line) 
{
    uint128_t i, max;
    i = *line;
    max = i + LINES_PER_THREAD;

    // check if previous match this one
    // it means we are in the middle of a group
    // if so, skip this one.

    uint128_t current = table[i].key;
    uint128_t current_m = current & mask;
    uint128_t c_exp = table[i].exp1;

    for (i=i+1; i < max; i++)
    {
        if ((table[i].key & mask) == current_m)
        {
            table[i].key ^= current;
            table[i].exp2 = c_exp;
        }
        else
        {
            current = table[i].key;
            current_m = current & mask;
            c_exp = table[i].exp1;
        }
    }

    // continue until group ends
}

static void sort(uint128_t *line) {
    qsort(table + (*line), LINES_PER_THREAD, sizeof(entry), sorting_function);
}


static void init_array() {
    int i;
    for ( i = 0; i < THREADS; i++) {
        args[i] = i * LINES_PER_THREAD;
        pthread_create(pids + i, NULL, (void *)(&generate), args + i);
    }
}

static void sort_blocks() {
    int i;
    for ( i = 0; i < THREADS; i++) {
        args[i] = i * LINES_PER_THREAD;
        pthread_create(pids + i, NULL, (void *)(&sort), args + i);
    }
}

static void first_match_array() {
    int i;
    for (i = 0; i < THREADS; i++) {
        args[i] = i * LINES_PER_THREAD;
        pthread_create(pids + i, NULL, (void *)(&match), args + i);
    }
}

static void wait_all() {
    for (int i = 0; i < THREADS; i++) {
        pthread_join(pids[i], NULL);
    }
}


void print_poly(uint128_t mask) {
    for(int i = 127; i >= 1; i--) {
        if ((mask >> i) & 0x1) {
            printf("x^%d+", i);    
        }
    }
    printf("1\n");
}

void dump_table() {
    for(uint64_t i = 1; i < SIZE; i++) {
        printf("%016" PRIx64 "%016" PRIx64 "\n", (uint64_t)(table[i].key >> 64), (uint64_t)table[i].key);
    }
}

int main(int argc, char *argv[])
{
    
    table = malloc(sizeof(entry)*(SIZE));
    if (table == NULL) {
        fprintf(stderr, "Fatal: failed to allocate %" PRIu64 " bytes.\n", SIZE);
        abort();
    }
    mask = random_mask(0, DEGREE/3);
    clamp_mask = random_mask(mask, CLAMP);
    printf("Polynominal: "); print_poly(POLY);
    printf("Threads:     %d\n", THREADS);
    printf("Mask:        0x%016" PRIx64 "%016" PRIx64 "\n", (uint64_t)(mask >> 64), (uint64_t)mask);
    printf("Clamp mask:  0x%016" PRIx64 "%016" PRIx64 "\n\n", (uint64_t)(clamp_mask >> 64), (uint64_t)clamp_mask);

    // for each thread, compute x^{size/threads} mod P(x)
    // and generate from there.
    printf("Generating %" PRIu64 " initial entries... ", SIZE); fflush(stdout);
    clock_t start = clock(), diff;
    init_array();
    wait_all();
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("(%d s)\n", msec/1000);

    printf("Sorting...  "); fflush(stdout);
    start = clock(), diff;
    sort_blocks();
    wait_all();
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("(%d s)\n", msec/1000);

/*
    // TODO: parallelize
    printf("Sorting...  "); fflush(stdout);
    start = clock();
    qsort(table, SIZE, sizeof(entry), sorting_function);
    //parallel_sort(table, N, n, sorting_function_asc, sorting_function_desc);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("(%d s)\n", msec/1000);

    printf("Matching... "); fflush(stdout);
    start = clock();
    first_match_array();
    wait_all();
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("(%d s)\n", msec/1000);

    // TODO: parallelize
    printf("Sorting...  "); fflush(stdout);
    start = clock();
    qsort(table, SIZE, sizeof(entry), sorting_function2);
    //parallel_sort(table, N, n, sorting_function_asc, sorting_function_desc);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("(%d s)\n", msec/1000);


    // TODO: parallelize
    // this has determinstic size, so just exclude pivots
    // and start each thread from an equidistant set of indices
    // i.e. exp(2, n/nthreads*thread, poly)
    // consideration: instead of pairing with pivot, we
    // can generate all pairs --> non-determinstic size

    printf("Matching... \n\n");
    uint128_t current;
    uint64_t current_i = 0;
    for(uint64_t i = 1; i < SIZE; i++)
    {
        if ((table[i].key == current))
        {
            printf("x^%" PRIu64 " + x^%" PRIu64 " + x^%" PRIu64 " + x^%" PRIu64 "\n", 
                table[i].exp1, table[i].exp2, 
                table[current_i].exp1, 
                table[current_i].exp2);
        }
        else
        {
            current = table[i].key;
            current_i = i;
        }
    }*/

    free(table);
}