typedef unsigned __int128 uint128_t;

uint128_t random_mask(uint128_t mask, int weight, int degree, unsigned int seed) 
{
    uint64_t index = 0;
    uint128_t new_mask = 0;
    srand(seed);

    for(int i = 0; i < weight; ++i)
    {
        int t = 0;
        do
        {
            index = rand() % degree;
            t = ((new_mask >> index) & (uint128_t)0x1);
            t |= ((mask >> index) & (uint128_t)0x1);
        }
        while (t);

        new_mask = new_mask | ((uint128_t)1 << index);
    }

    return (uint128_t)new_mask;
}
