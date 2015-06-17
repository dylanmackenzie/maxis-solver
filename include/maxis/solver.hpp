#include "graph.hpp"
#include "genetic.hpp"

#ifndef MAXIS_SOLVER_H
#define MAXIS_SOLVER_H

namespace maxis {

// Solver interface
class MaxisSolver {
public:
    virtual ~MaxisSolver() {};

    virtual BitVector operator()() =0;
};

// Solver implementations
// ======================

// The brute force solver is a naive O(2^n) algorithm. It is included
// only as a baseline for other solver implementations.
class BruteForceMaxisSolver : public virtual MaxisSolver {
public:
    BruteForceMaxisSolver(const Graph&);
    virtual BitVector operator()();

private:
    const Graph &graph;
};

// The genetic solver returns an approximation of the maximum weighted
// indepenent set using a genetic algorithm based on the one by Beasley
// and Chu. It accepts various implementations of each genetic operator.
class GeneticMaxisSolver : public virtual MaxisSolver {
public:
    GeneticMaxisSolver(
            const Graph&,
            genetic::Selector&,
            genetic::Recombinator&,
            genetic::Mutator&);
    virtual BitVector operator()();

    static void initialize_set(size_t, BitVector&, BitVector&);
    static void heuristic_feasibility(size_t, BitVector&, BitVector&);

    // constraint is the fitness goal for the genetic algorithm.
    // The solver will return immediately when a solution with the given
    // fitness(weight) is found.
    double constraint;

    // size is the size of the population
    size_t size;

protected:
    Graph graph;
    genetic::Selector &selector;
    genetic::Recombinator &recombinator;
    genetic::Mutator &mutator;

    // Stores the new order of the vertices after they sorted by weight
    // and degree
    std::vector<size_t> permutation;
};

// The parallel version of the genetic solver. It keeps a segment of the
// population on each core, and performs a traditional genetic algorithm
// on each sub-population. At the end of {migration_period} iterations,
// population members are moved from core to core. No attempt is made to
// eliminate duplicates at the global level, only individual
// sub-populations are guaranteed to be unique.
class ParallelGeneticMaxisSolver : public GeneticMaxisSolver {
public:
    ParallelGeneticMaxisSolver(
            const Graph&,
            genetic::Selector&,
            genetic::Recombinator&,
            genetic::Mutator&);
    virtual BitVector operator()();

    size_t migration_period;
};

} // namespace maxis

#endif // MAXIS_SOLVER_H
