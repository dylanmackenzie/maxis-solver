#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <unordered_set>

#include "graph.hpp"
#include "genetic.hpp"
#include "hash_table.hpp"

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
    using PopIter = std::vector<genetic::Phenotype>::iterator;

    GeneticMaxisSolver(const Graph&, genetic::AlgorithmStrategy);

    virtual BitVector operator()();

    static void iterate(
            const Graph&, PopIter, PopIter,
            genetic::AlgorithmState&,
            genetic::AlgorithmStrategy&,
            BitVectorHashTable&);

    static void initialize_set(const Graph&, BitVector&);
    static void heuristic_feasibility(const Graph&, BitVector&);

    // constraint is the fitness goal for the genetic algorithm.
    // The solver will return immediately when a solution with the given
    // fitness(weight) is found.
    double constraint;

    // size is the size of the population
    size_t size;

protected:
    Graph graph;
    genetic::AlgorithmStrategy strategy;

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
    ParallelGeneticMaxisSolver(const Graph&, genetic::AlgorithmStrategy);
    virtual BitVector operator()();

    unsigned int migration_period;
};

class WorkerSynchronizer {
public:
    // The number of workers that will be executing is passed to the
    // constructor. The number of threads attempting to use the
    // synchronizer SHOULD BE EXACTLY THIS NUMBER.
    WorkerSynchronizer(unsigned int);

    // wait_on_cycle blocks until all workers are complete. Once it
    // returns, no new handles can be created until next_cycle is
    // called.
    void wait_on_cycle();

    // next_cycle is called to start a new work cycle.
    void next_cycle();

    // This class allows us to use RAII to synchronize the workers. A
    // worker constructs a handle for the current work cycle and
    // destroys it when it is done with that work cycle. The worker then
    // increments the number of cycles it has completed, and constructs
    // a handle for the next work cycle. That constructor will block
    // until the remaining threads complete the current work cycle and
    // next_cycle is called.
    class Handle {
    public:
        Handle(WorkerSynchronizer&, unsigned int);
        ~Handle();
        Handle(const Handle&) =delete;
        Handle& operator=(const Handle&) =delete;
        // Move constructors will be disabled automatically

    private:
        WorkerSynchronizer &sync;
    };

    struct SynchronizationError : public std::runtime_error {
        SynchronizationError(const char* s) : std::runtime_error(s) {};
    };

private:
    // available_handles is the total number of handles available for a
    // given work cycle
    const unsigned int available_handles;
    std::atomic<unsigned int> utilized_handles;  // number of successfully constructed handles
    std::atomic<unsigned int> discarded_handles; // number of successfully destructed handles
    std::atomic<unsigned int> current_cycle;

    std::condition_variable worker_cv;
    std::mutex worker_lock;
    std::condition_variable manager_cv;
    std::mutex manager_lock;
};

class ParallelGeneticWorker {
public:
    using PopIter = std::vector<genetic::Phenotype>::iterator;
    ParallelGeneticWorker(
            const Graph&, unsigned int,
            PopIter, PopIter, genetic::AlgorithmState&,
            genetic::AlgorithmStrategy, WorkerSynchronizer&);

    virtual void operator()();

private:
    const Graph &graph;
    unsigned int migration_period;
    PopIter b;
    PopIter e;
    genetic::AlgorithmState &state;
    genetic::AlgorithmStrategy strategy;
    WorkerSynchronizer &sync;
};

} // namespace maxis

#endif // MAXIS_SOLVER_H
