#include "Abstract_Classes.h"
#include <stdlib.h>


// Uniformly distributed bond strengths
class uniform_Strength: public Strength
{
private:
    // Random number generator seed
    unsigned long long rdtsc(){
        unsigned int lo,hi;
        __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
        return ((unsigned long long)hi << 32) | lo;
    }

public:
    uniform_Strength()
    {
        srand48(rdtsc());
    }

    double get_new_strength(void)
    {
        return drand48();
    }
};

// Class with sequentially increasing strengths for testing
class testing_Strength: public Strength
{
private:
    int next_strength;

public:
    testing_Strength()
    {
        next_strength=0;
    }

    double get_new_strength()
    {
        return next_strength++;
    }

};
