//#include "sort.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>

typedef unsigned int uint128_t __attribute__((mode(TI)));
typedef __int128 int128_t;

typedef struct {
  uint128_t key;
  uint64_t exp1, exp2;
} __attribute__((packed)) entry;

struct timeval startwtime, endwtime;
double seq_time;

uint64_t N, n; 
entry* a; 

int threadlayers;

static uint128_t mask;

const int ASCENDING = 1;
const int DESCENDING = 0;

void init(void);
void print(void);
void test(void);
inline void exchange(uint64_t i, uint64_t j);
void compare(uint64_t i, uint64_t j, int dir);
void bitonic_merge(uint64_t lo, uint64_t cnt, int dir);
void parallel_sort(void);
void* sort(void* arg);
void* merge(void* arg);

int sorting_function_desc(const void *a, const void *b) 
{ 
    entry *ia = (entry *)a;
    entry *ib = (entry *)b;

    if ((ia->key ) > (ib->key ))
    {
        return -1;
    }
    else if ((ia->key )  < (ib->key ) )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int sorting_function_asc(const void *a, const void *b) 
{ 
    entry *ia = (entry *)a;
    entry *ib = (entry *)b;

    if ((ia->key ) > (ib->key ))
    {
        return 1;
    }
    else if ((ia->key )  < (ib->key) )
    {
        return -1;
    }
    else
    {
        return 0;
    }
}


void test() {
  int pass = 1;
  uint64_t i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1].key <= a[i].key);
  }

  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}

/** the main program **/
int main(int argc, char** argv)
{

    if (argc != 3 || atoi(argv[2]) > 256)
    {
        printf("Usage: %s q t\n  where n=2^q is problem size (power of two), and t is the number of threads, <=256, to use.\n", argv[0]);
        exit(1);
    }

    N = 1 << atoi(argv[1]);
    n = atoi(argv[2]);

    threadlayers = atoi(argv[2]);

    if (threadlayers != 0 && threadlayers != 1)
    {
        --threadlayers;
    }

    a = (entry*)malloc(N * sizeof(entry));

    init();
    gettimeofday(&startwtime, NULL);
    parallel_sort();
    gettimeofday(&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    printf("Bitonic parallel recursive with qsort and %i threads wall clock time = %f\n", 1 << atoi(argv[2]), seq_time);
    test();

}

void init()
{
    uint64_t i;

    for (i = 0; i < N; i++)
    {
        a[i].key = rand() % N;
    }
}

inline void exchange(uint64_t i, uint64_t j)
{
    entry t;
    t = a[i];
    a[i] = a[j];
    a[j] = t;
}


inline void compare(uint64_t i, uint64_t j, int dir)
{
    if (dir == (a[i].key > a[j].key))
    {
        exchange(i, j);
    }
}

void bitonic_merge(uint64_t lo, uint64_t cnt, int dir)
{
    if (cnt > 1) 
    {
        uint64_t i, k = cnt / 2;

        for (i = lo; i < lo + k; i++)
        {
            compare(i, i + k, dir);
        }

        bitonic_merge(lo, k, dir);
        bitonic_merge(lo + k, k, dir);
    }
}

typedef struct {
    uint64_t lo, cnt, layer;
    int dir;
    const void *asc, *desc;
} sarg;


void* merge(void* arg)
{
    uint64_t lo = ((sarg*)arg)->lo;
    uint64_t cnt = ((sarg*)arg)->cnt;
    uint64_t layer = ((sarg*)arg)->layer;

    int dir = ((sarg*)arg)->dir;

    if (cnt > 1)
    {
        uint64_t i, k = cnt / 2;

        for (i = lo; i < lo + k; ++i)
        {
            compare(i, i + k, dir);
        }

        if (layer <= 0)
        {
            bitonic_merge(lo, k, dir);
            bitonic_merge(lo + k, k, dir);
        
            return 0;
        }

        sarg arg1, arg2;
        pthread_t thread1, thread2;

        arg1.lo = lo;
        arg1.cnt = k;
        arg1.dir = dir;
        arg1.layer = layer - 1;
        arg2.lo = lo + k;
        arg2.cnt = k;
        arg2.dir = dir;
        arg2.layer = layer - 1;

        pthread_create(&thread1, NULL, merge, &arg1);
        pthread_create(&thread2, NULL, merge, &arg2);

        pthread_join(thread1, NULL);
        pthread_join(thread2, NULL);
    }

    return 0;
}


void* sort(void* arg)
{
    uint64_t lo = ((sarg*)arg)->lo;
    uint64_t cnt = ((sarg*)arg)->cnt;
    uint64_t layer = ((sarg*)arg)->layer;

    int dir = ((sarg*)arg)->dir;

    if (cnt > 1)
    {
        uint64_t k = cnt / 2;

        if (layer >= threadlayers)
        {
            qsort(a + lo, k, sizeof(entry), ((sarg*)arg)->asc);
            qsort(a + (lo + k), k, sizeof(entry), ((sarg*)arg)->desc);
        }
        else
        {
            sarg arg1, arg2;
            pthread_t thread1, thread2;

            arg1.lo = lo;
            arg1.cnt = k;
            arg1.dir = ASCENDING;
            arg1.layer = layer + 1;
            arg1.asc = sorting_function_asc;
            arg1.desc = sorting_function_desc;

            arg2.lo = lo + k;
            arg2.cnt = k;
            arg2.dir = DESCENDING;
            arg2.layer = layer + 1;
            arg2.asc = sorting_function_asc;
            arg2.desc = sorting_function_desc;

            pthread_create(&thread1, NULL, sort, &arg1);
            pthread_create(&thread2, NULL, sort, &arg2);

            pthread_join(thread1, NULL);
            pthread_join(thread2, NULL);
        }
        
        sarg arg3;
        
        arg3.lo = lo;
        arg3.cnt = cnt;
        arg3.dir = dir;
        arg3.layer = threadlayers - layer;

        merge(&arg3);
    }

    return 0;
}

void parallel_sort()
{
    sarg arg;
    arg.lo = 0;
    arg.cnt = N;
    arg.dir = ASCENDING;
    arg.layer = 0;
    arg.asc = sorting_function_asc;
    arg.desc = sorting_function_desc;

    sort(&arg);
}