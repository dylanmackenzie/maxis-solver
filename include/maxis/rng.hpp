#include <random>

#ifndef MAXIS_RNG_HPP
#define MAXIS_RNG_HPP

// Mixin for seeding and generating random doubles between 0 and 1
class RNG {
public:
    RNG() : engine{0}, probability{0, 1} {};

    template<typename SeedType>
    void seed_random(SeedType seed) { engine.seed(seed); };

    inline double random_prob() { return probability(engine); };
    inline size_t random_index(size_t start, size_t end) {
        std::uniform_int_distribution<size_t> range{start, end};
        return range(engine);
    };

private:
    std::mt19937 engine;
    std::uniform_real_distribution<double> probability;
};

#endif
