#ifndef RNG_H
#define RNG_H
#include <random>

struct RNG {
    // Mersenne Twister 19937 random number generator
    std::mt19937 gen;

    explicit RNG(const uint32_t seed = std::random_device{}()) : gen(seed) {}

    double normal(const double mean, const double stddev) {
        std::normal_distribution<double> dist(mean, stddev);
        return dist(gen);
    }
};

#endif //RNG_H
