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
} entry;

#define CLAMP 0
#define DEGREE (72-CLAMP)
#define THREADS 11
#define SUBLISTS (THREADS * THREADS)
#define SIZE ((uint64_t)(((((uint64_t)1 << (DEGREE/3 + 1)) * 1.5) / SUBLISTS) * SUBLISTS))
#define SUBLIST_SIZE (SIZE/SUBLISTS)
#define LINES_PER_THREAD (SIZE/THREADS)

typedef struct {
  uint64_t u[THREADS];
  uint64_t l[THREADS];
  entry entries[SUBLIST_SIZE];
} sublist;

#define POLY ((((uint128_t) 257) << 64) + (uint128_t)(299733420677929))
//#define POLY ((((uint128_t) 0x212011231) << 64) + 0x2fea1a9dc028a229)
//#define POLY ((((uint128_t) 0x2090000) << 64) + 0xa0048000000003)

static uint128_t mask, clamp_mask;
static uint128_t f = (((uint128_t) 1) << (DEGREE+CLAMP));
static sublist* table;

static pthread_t pids[THREADS];
static int args[THREADS];

int sorting_function(const void *a, const void *b) 
{ 
    entry *ia = (entry *)a;
    entry *ib = (entry *)b;

    if ((ia->key & mask) > (ib->key & mask)){
        return -1;
    } else if ((ia->key & mask) < (ib->key & mask) ){
        return 1;
    } else {
        return 0;
    }
}

int sorting_function2(const void *a, const void *b) 
{ 
    entry *ia = (entry *)a;
    entry *ib = (entry *)b;

    if (ia->key > ib->key){
        return -1;
    } else if (ia->key < ib->key){
        return 1;
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

inline int phi(uint128_t x) {
    return (x) % THREADS;
}

static void generate(int *line)
{
    uint64_t i = *line, j = 0;
    uint64_t s, t, offset = SUBLIST_SIZE/THREADS * (*line);
    uint128_t px;
    px = gf2_exp(2, i * LINES_PER_THREAD, f, POLY, DEGREE);
    while(j < LINES_PER_THREAD)
    {
        if ((px & clamp_mask) == 0)
        {
            // this is the sublist
            t = phi(px & mask);
            if (table[t].u[*line] < (SUBLIST_SIZE/THREADS)-1) {
                s = table[t].u[*line] + offset;
                table[t].u[*line]++;
                table[t].entries[s].key = px;
                table[t].entries[s].exp1 = i;
                //table[t].entries[s].exp2 = 1;
            }
            j++;
        }

        i++;
        px <<= 1;
        if ((px & f) != 0){
            px ^= POLY;
        }
    }
}

static void reduce(int *line)
{
    /*int64_t j = 0;  // this is the error if unsigned is used
    for(int x = 0; x < THREADS; x++) {
        j += table[*line].u[x];
    }
    int64_t k = j - 1;
    uint128_t cur = table[*line].entries[k].key;
    uint64_t cur_exp = table[*line].entries[k].exp1;

    while(j >= 0)
    {
        if (table[*line].entries[j].key == cur) {
            table[*line].entries[j].key ^= cur;
            table[*line].entries[j].exp2 = cur_exp;
        } else {
            // clear pivot
            table[*line].entries[k].key = 0;
            table[*line].entries[k].exp1 = 0;
            // set new pivot
            k = j;
            cur = table[*line].entries[j].key;
            cur_exp = table[*line].entries[j].exp1;
        }
        j--; // this is the error
    }*/
}


static void sort(int *line) {
    uint64_t size = 0;
    for(int x = 0; x < THREADS; x++) {
        size += table[*line].u[x];
    }
    qsort(table[*line].entries, size, sizeof(entry), sorting_function);
}

static void wait_all() {
    for (int i = 0; i < THREADS; i++) {
        pthread_join(pids[i], NULL);
    }
}

static void init_array() {
    int i;
    for ( i = 0; i < THREADS; i++) {
        args[i] = i;
        pthread_create(pids + i, NULL, (void *)(&generate), args + i);
    }
}

static void sort_blocks() {
    int i;
    for ( i = 0; i < THREADS; i++) {
        args[i] = i;
        pthread_create(pids + i, NULL, (void *)(&sort), args + i);
    }
}

static void reduce_blocks() {
    int i;
    for ( i = 0; i < THREADS; i++) {
        args[i] = i;
        pthread_create(pids + i, NULL, (void *)(&reduce), args + i);
    }
}

/*
static void first_match_array() {
    int i;
    for (i = 0; i < THREADS; i++) {
        args[i] = i * LINES_PER_THREAD;
        pthread_create(pids + i, NULL, (void *)(&match), args + i);
    }
}*/



void clear() {
    for(int i = 0; i < SUBLISTS; i++) {
        for(int j= 0; j < THREADS; j++) {
            table[i].u[j] = 1;
            table[i].l[j] = 1;
        }
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
    for(int j= 0; j < SUBLISTS; j++) {
        printf("sublist: %d %d %d\n", j, table[j].u[0], table[j].l[0]);
        for(uint64_t i = 0; i < SUBLIST_SIZE; i++) {
            printf("%016" PRIx64 "%016" PRIx64 " %d %d\n", (uint64_t)(table[j].entries[i].key >> 64) , (uint64_t)(table[j].entries[i].key >> 0), table[j].entries[i].exp1, table[j].entries[i].exp2);
        }
        printf("----\n");
    }
}

int main(int argc, char *argv[])
{
    
    table = calloc((SUBLISTS), sizeof(sublist));
    if (table == NULL) {
        fprintf(stderr, "Fatal: failed to allocate %" PRIu64 " bytes.\n", SIZE);
        abort();
    }
    mask = random_mask(0, DEGREE/3);
    clamp_mask = random_mask(mask, CLAMP);
    printf("Polynominal: "); print_poly(POLY);
    printf("Threads:     %d\n", THREADS);
    printf("Sublists:    %d\n", SUBLISTS);
    printf("Mask:        0x%016" PRIx64 "%016" PRIx64 "\n", (uint64_t)(mask >> 64), (uint64_t)mask);
    printf("Clamp mask:  0x%016" PRIx64 "%016" PRIx64 "\n\n", (uint64_t)(clamp_mask >> 64), (uint64_t)clamp_mask);

    // for each thread, compute x^{size/threads} mod P(x)
    // and generate from there.
    //clear();
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
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("(%d s)\n", msec/1000);

    //dump_table();
    printf("Reducing...  "); fflush(stdout);
    start = clock(), diff;
    reduce_blocks();
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