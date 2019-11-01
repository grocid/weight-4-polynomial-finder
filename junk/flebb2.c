#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <pthread.h>

typedef unsigned int uint128_t __attribute__((mode(TI)));
typedef __int128 int128_t;

typedef struct {
  uint128_t key;
  uint64_t exp1;
  uint64_t exp2;
} __attribute__((packed)) entry;

#define DEGREE 72
#define THREADS 6
#define LINES_PER_THREAD (size/THREADS)
#define size (uint64_t)(((uint64_t)1 << (DEGREE/3 + 1)) * 1.5)
#define POLY ((((uint128_t) 257) << 64) + (uint128_t)(299733420677929))

static uint128_t mask;
static uint128_t f = (((uint128_t) 1) << DEGREE);
static entry* table;

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
        printf("%u", (mask >> i) & 0x1);
    }
    printf("\n");
}



/*
    GF(2) polynomial multiplication mod P(x)
*/
uint128_t multiply(uint128_t a, uint128_t b)
{
    uint128_t r = 0;
    uint128_t q = b;
    for(int i = 0; i < DEGREE; i++)
    {
        if((a & ((uint128_t)1 << i)) != 0)
        {
            r ^= q;
        }
        q = (q << 1);
        if ((q & f) != 0)
        {
            q ^= POLY;
        }
    }
    return r;
}

/*
    GF(2) polynomial exponentiation mod P(x)
*/
uint128_t gf2_exp(uint128_t x, uint128_t n)
{
    if(n == 0)
    {
        return 1;
    }
    uint128_t y = 1;
    while(n > 1)
    {
        if(n % 2 == 0)
        {
            x = multiply(x, x);
            n = n / 2;
        }
        else
        {
            y = multiply(x, y);
            x = multiply(x, x);
            n = (n - 1) / 2;
        }
    }
    return multiply(x, y);
}

/*
    generates a random mask of weight degree/3 + 2
*/
uint128_t random_mask() 
{
    int i1 = 0, i2 = 0, b;
    int128_t mask = 0;//(1 << (DEGREE/3)) - 1;
    srand(0*time(NULL)+1);

    for(int i = 0; i < DEGREE/3; ++i)
    {
        do {
            i1 = rand() % 31;
        }
        while ((mask >> i1) & (uint128_t)0x1);

        mask = mask | ((uint128_t)1 << i1);
    }

    return (uint128_t)mask;
}


static pthread_t pids[THREADS];
static int args[THREADS];

static void generate(int *line)
{
    int i, max;
    uint128_t px = 1;
    i = *line;
    max = i + LINES_PER_THREAD;
    px = gf2_exp(2, i);
    for (; i < max; i++)
    {
        px <<= 1;
        if ((px & f) != 0){
            px ^= POLY;
        }
        table[i].key = px;
        table[i].exp1 = i;
    }
}

static void match(int *line) 
{
    int i, max;
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
/*
            while ((table[i].key & 1) == 0)
            {
                table[i].key >>= 1;
                table[i].exp1--;
                table[i].exp2--;
            }*/
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


static void init_array() {
     int i;
    for ( i = 0; i < THREADS; i++) {
        args[i] = i* LINES_PER_THREAD;
        pthread_create(pids + i, NULL, &generate, args + i);;
    }
}

static void first_match_array() {
     int i;
    for ( i = 0; i < THREADS; i++) {
        args[i] = i* LINES_PER_THREAD;
        pthread_create(pids + i, NULL, &match, args + i);;
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
            printf("x^%u+", i);    
        }
    }
    printf("1\n");
}

int main(int argc, char *argv[])
{
    
    table = malloc(sizeof(entry)*(size));  
    mask = random_mask();
    printf("Polynominal: ");
    print_poly(POLY);
    printf("Threads:     %d\n", THREADS);
    printf("Mask:        %lx%lx\n\n", mask);

    printf("Generating %d initial entries... ", size); fflush(stdout);

    // TODO: parallelize
    // for each thread, compute x^{size/threads} mod P(x)
    // and generate from there.
    clock_t start = clock(), diff;
    init_array();
    wait_all();
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("(%d s)\n", msec/1000);

    // TODO: parallelize
    printf("Sorting...  "); fflush(stdout);
    start = clock();
    qsort(table, size, sizeof(entry), sorting_function);
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
    qsort(table, size, sizeof(entry), sorting_function2);
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
    uint32_t current_i = 0;
    for(uint32_t i = 1; i < size; i++)
    {
        if (table[i].key == current)
        {
            if ((table[i].exp1 != table[current_i].exp1) && (table[i].exp1 != table[current_i].exp2)) {
            printf("x^%u + x^%u + x^%u + x^%u\n", 
                table[i].exp1, table[i].exp2, 
                table[current_i].exp1, 
                table[current_i].exp2);}
        }
        else
        {
            current = table[i].key;
            current_i = i;
        }
    }

    free(table);
}