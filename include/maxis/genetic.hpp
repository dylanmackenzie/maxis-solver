#include <random>
#include <memory>

#include "graph.hpp"

#ifndef MAXIS_GENETIC_H
#define MAXIS_GENETIC_H

namespace maxis {
namespace genetic {

struct Phenotype {
    Phenotype() : chromosome{nullptr}, fitness{0.0} {}
    Phenotype(size_t n) : chromosome{new BitVector(n)}, fitness{0.0} {}

    std::unique_ptr<BitVector> chromosome;
    double fitness;

    inline bool operator<(const Phenotype &cmp) const {
        return fitness < cmp.fitness;
    }
};

// Mixin for seeding and generating random doubles between 0 and 1
class RNG {
public:
    RNG() : engine{0}, probability{0, 1} {};

    template<typename SeedType>
    void seed_random(SeedType seed) { engine.seed(seed); };

    double random_prob();
    size_t random_index(size_t, size_t);

private:
    std::mt19937 engine;
    std::uniform_real_distribution<double> probability;
};

// Abstract Base Classes

class Generator {

public:
    virtual ~Generator() {};
    virtual void initialize(Phenotype&) =0;
};

class Selector {
public:
    virtual ~Selector() {};
    virtual Phenotype& select(std::vector<Phenotype>::iterator, std::vector<Phenotype>::iterator, double, double) =0;
};

class Recombinator {
public:
    virtual ~Recombinator() {};
    virtual void crossover(Phenotype&, Phenotype&, Phenotype&) =0;
};

class Mutator {
public:
    virtual ~Mutator() {};
    virtual void mutate(Phenotype&) =0;
};

// Generator Implementations

class MaxisHeuristicGenerator : public virtual Generator {
public:
    MaxisHeuristicGenerator();

    virtual void initialize(Phenotype &);

    const Graph *graph;
    RNG rng;
};

// Selector implementations

class TournamentSelector : public virtual Selector {
public:
    TournamentSelector(size_t size) : size{size} {};
    Phenotype& select(std::vector<Phenotype>::iterator, std::vector<Phenotype>::iterator, double, double);
    RNG rng;

private:
    size_t size;
};

class RouletteSelector : public virtual Selector {
public:
    Phenotype& select(std::vector<Phenotype>::iterator, std::vector<Phenotype>::iterator, double, double);
    RNG rng;
private:
};

// Recombinator implementations

// class SinglePointRecombinator : public virtual Recombinator {
// public:
    // void crossover(Phenotype&, Phenotype&);
    // RNG rng;
// };

// class DoublePointRecombinator : public virtual Recombinator {
// public:
    // void crossover(Phenotype&, Phenotype&);
    // RNG rng;
// };

class BlendingRecombinator : public virtual Recombinator {
public:
    void crossover(Phenotype&, Phenotype&, Phenotype&);
    RNG rng;
};

// Mutator implementations

class SinglePointMutator : public virtual Mutator {
public:
    SinglePointMutator(double mul);
    void mutate(Phenotype&);

    RNG rng;
    double probability_multiplier;
};

class MaxisHeuristicMutator : public virtual Mutator {
public:
    MaxisHeuristicMutator();
    void mutate(Phenotype&);

    const Graph *graph;
    RNG rng;
};

} // namespace genetic
} // namespace maxis

#endif // MAXIS_GENETIC_H
