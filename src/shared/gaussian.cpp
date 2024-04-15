#include <random>
#include <iostream>

// fixed random seed
static std::default_random_engine e(11);

float gaussian(float mean,float sigma)
{
    std::normal_distribution<float> n(mean,sigma);
    return n(e);
}