#include <algorithm>
#include <memory>

#include "graph.hpp"
#include "rng.hpp"

#ifndef MAXIS_GENETIC_H
#define MAXIS_GENETIC_H

namespace maxis {
namespace genetic {

// Phenotype is a container which holds a chromosome and its associated
// fitness.
struct Phenotype {
    Phenotype(BitVector *ch) : chromosome{ch}, fitness{0.0} {}

    inline bool operator<(const Phenotype &cmp) const {
        return fitness < cmp.fitness;
    }

    // chromosome should be a unique_ptr, but the hash table for
    // checking duplicate chromosomes is more easily implemented with
    // pointers to chromosomes, and the overhead of a shared_ptr is
    // unnecessary.
    BitVector *chromosome;
    double fitness;
};

// AlgorithmState is passed to the genetic operators on every
// invocation. It is kept up to date by the solver.
struct AlgorithmState {
    AlgorithmState();
    void initialize(std::vector<Phenotype>::iterator, std::vector<Phenotype>::iterator);

    double min_fitness;
    double max_fitness;
    double total_fitness;
    double adjusted_fitness; // sum(fitness - min_fitness) over the population

    unsigned int iterations;
    size_t size; // current size of population
};

// Genetic Operator Interfaces
// ===========================
// These interfaces allow us to plug different genetic operators into
// the solver for testing. If performance is critical, the solver
// function can be templatized and these interface classes can be
// removed.

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

// AlgorithmStrategy packages a group of genetic operators
struct AlgorithmStrategy {
    AlgorithmStrategy(Selector &sel, Recombinator &rec, Mutator &mut)
        : selector{sel}, recombinator{rec}, mutator{mut} {};

    Selector &selector;
    Recombinator &recombinator;
    Mutator &mutator;
};

// Selector implementations

// The tournament selector picks `size` cantidates out of the population
// at random and returns the best one.
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

// The roulette selector picks a population member with a probability
// proportional to the adjusted fitness (fitness - min_fitness) of the
// member. Note that if the genetic algorithm is near convergence and
// the population is equally fit with the exception of one fitter
// outlier, that outlier will be selected every time.
class RouletteSelector : public virtual Selector {
public:
    Phenotype& select(
            const AlgorithmState&,
            std::vector<Phenotype>::iterator,
            std::vector<Phenotype>::iterator);
    RNG rng;
};

// Recombinator implementations

// The blending recombinator implements the genetic-fusion operator
// described by Beasley and Chu. It randomly selects individul genes
// from each parent with a probability proportional to the fitness of
// the parent.
class BlendingRecombinator : public virtual Recombinator {
public:
    void breed(const AlgorithmState&, const Phenotype&, const Phenotype&, Phenotype&);
    RNG rng;
};

// Mutator implementations

// The simple mutator maintains a constant average mutation rate of
// {mul} bits.
class SimpleMutator : public virtual Mutator {
public:
    SimpleMutator(double mul) : probability_multiplier{mul} {};
    void mutate(const AlgorithmState&, Phenotype&);

    RNG rng;

private:
    double probability_multiplier;
};

// The variable rate mutator uses a sigmoid function to set the mutation
// rate. {ss} is the mutation rate as the number of iterations
// approaches infinity. {halfway} is the number of iterations before the
// mutation rate reaches half of its steady state value. {gradient} is the
// slope at the halfway point.
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
} // namespace maxis

#endif // MAXIS_GENETIC_H
