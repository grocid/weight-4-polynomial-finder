//#include "sort.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>

typedef unsigned int uint128_t __attribute__((mode(TI)));
typedef struct {
  uint128_t key;
  uint64_t exp1, exp2;
} __attribute__((packed)) entry;
uint64_t N, n; 
entry* a; 

int threadlayers;

const int ASCENDING = 1;
const int DESCENDING = 0;

void bitonic_merge(uint64_t lo, uint64_t cnt, int dir);
void* sort(void* arg);
void* merge(void* arg);

inline void exchange(uint64_t i, uint64_t j);
inline void compare(uint64_t i, uint64_t j, int dir);

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
            arg1.asc = ((sarg*)arg)->asc;
            arg1.desc = ((sarg*)arg)->desc;

            arg2.lo = lo + k;
            arg2.cnt = k;
            arg2.dir = DESCENDING;
            arg2.layer = layer + 1;
            arg2.asc = ((sarg*)arg)->asc;
            arg2.desc = ((sarg*)arg)->desc;

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

void parallel_sort(entry* a, uint64_t size, int threads, int (*)(const void *, const void *) sorting_function_asc, int (*)(const void *, const void *)sorting_function_desc)
{
    threadlayers = threads;
    N = size;

    if (threadlayers != 0 && threadlayers != 1)
    {
        --threadlayers;
    }

    sarg arg;
    arg.lo = 0;
    arg.cnt = N;
    arg.dir = ASCENDING;
    arg.layer = 0;
    arg.asc = sorting_function_asc;
    arg.desc = sorting_function_desc;

    sort(&arg);
}