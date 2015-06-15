#include <algorithm>
#include <memory>

#include "graph.hpp"
#include "rng.hpp"

#ifndef MAXIS_GENETIC_H
#define MAXIS_GENETIC_H

namespace genetic {

using maxis::BitVector;

struct Phenotype {
    Phenotype(BitVector *ch) : chromosome{ch}, fitness{0.0} {}

    BitVector * chromosome;
    double fitness;

    inline bool operator<(const Phenotype &cmp) const {
        return fitness < cmp.fitness;
    }
};

struct AlgorithmState {
    AlgorithmState() :
        min_fitness{0}, max_fitness{0}, total_fitness{0},
        adjusted_fitness{0}, iterations{0} {};

    double min_fitness;
    double max_fitness;
    double total_fitness;
    double adjusted_fitness; // sum(fitness - min_fitness) over the population

    unsigned int iterations;
};

// Abstract Base Classes

class Selector {
public:
    virtual ~Selector() {};
    virtual Phenotype& select(
            const AlgorithmState&,
            std::vector<Phenotype>::iterator,
            std::vector<Phenotype>::iterator) =0;
};

class Recombinator {
public:
    virtual ~Recombinator() {};
    virtual void breed(
            const AlgorithmState&,
            const Phenotype&, const Phenotype&, Phenotype&) =0;
};

class Mutator {
public:
    virtual ~Mutator() {};
    virtual void mutate(const AlgorithmState&, Phenotype&) =0;
};

// Selector implementations

class TournamentSelector : public virtual Selector {
public:
    TournamentSelector(size_t size) : size{size} {};
    Phenotype& select(
            const AlgorithmState&,
            std::vector<Phenotype>::iterator,
            std::vector<Phenotype>::iterator);
    RNG rng;

private:
    size_t size;
};

class RouletteSelector : public virtual Selector {
public:
    Phenotype& select(
            const AlgorithmState&,
            std::vector<Phenotype>::iterator,
            std::vector<Phenotype>::iterator);
    RNG rng;
};

// Recombinator implementations

class BlendingRecombinator : public virtual Recombinator {
public:
    void breed(const AlgorithmState&, const Phenotype&, const Phenotype&, Phenotype&);
    RNG rng;
};

// Mutator implementations

class SimpleMutator : public virtual Mutator {
public:
    SimpleMutator(double mul) : probability_multiplier{mul} {};
    void mutate(const AlgorithmState&, Phenotype&);

    RNG rng;

private:
    double probability_multiplier;
};

class VariableRateMutator : public virtual Mutator {
public:
    VariableRateMutator(double ss, size_t halfway, double gradient)
        : steady_state{ss}, halfway{halfway}, gradient{gradient} {};

    void mutate(const AlgorithmState&, Phenotype&);

    RNG rng;

private:
    double steady_state;
    size_t halfway;
    double gradient;
};

} // namespace genetic

#endif // MAXIS_GENETIC_H
