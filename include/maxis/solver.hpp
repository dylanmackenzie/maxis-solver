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

class BruteForceMaxisSolver : public virtual MaxisSolver {
public:
    BruteForceMaxisSolver(const Graph&);
    virtual BitVector operator()();

private:
    const Graph &graph;
};

class BranchBoundMaxisSolver : public virtual MaxisSolver {
public:
    BranchBoundMaxisSolver(const Graph&);
    virtual BitVector operator()();

private:
    const Graph &graph;
};

class GeneticMaxisSolver : public virtual MaxisSolver {
public:
    GeneticMaxisSolver(
            const Graph&,
            genetic::Selector&,
            genetic::Recombinator&,
            genetic::Mutator&);
    virtual BitVector operator()();

    size_t size;
    double constraint;

private:
    Graph graph;
    genetic::Selector &selector;
    genetic::Recombinator &recombinator;
    genetic::Mutator &mutator;

    void compute_fitness(genetic::Phenotype&) const;

    // Stores the permutation for reordering the vertices after sorting
    std::vector<size_t> permutation;
};

} // namespace maxis

#endif // MAXIS_SOLVER_H
