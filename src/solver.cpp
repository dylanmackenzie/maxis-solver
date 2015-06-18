#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <csignal>
#include <iostream>
#include <limits>
#include <mutex>
#include <thread>
#include <unordered_set>
#include <vector>

#include "maxis/graph.hpp"
#include "maxis/genetic.hpp"
#include "maxis/solver.hpp"
#include "maxis/rng.hpp"

using std::vector;
using std::begin; using std::end;
using namespace maxis::genetic;

static volatile bool sigint_flag;

void sigint_handler(int signum) {
    sigint_flag = true;
}

namespace maxis {

// Brute Force Solver
// ==================

BruteForceMaxisSolver::BruteForceMaxisSolver(const Graph &g) : graph{g} {}

// Ripple carry increment for a collection of bits.
// Returns false on overflow.
bool
increment_bit_vector(BitVector &bv) {

    bool carry = true;
    for (auto it = begin(bv); it != end(bv); ++it) {
        if (*it == 0){
            *it = 1;
            carry = false;
            break;
        } else {
            *it = 0;
        }
    }

    return !carry;
}

BitVector
BruteForceMaxisSolver::operator()() {
    // max_weight tracks the current maximum weight
    int max_weight = 0;
    BitVector max_set;

    auto adjacency_matrix = graph.adjacency_matrix();

    BitVector bv(graph.order(), 0);
    do {
        if (graph.is_independent_set(bv)) {
            auto weight = graph.weighted_total(bv);
            if (weight > max_weight) {
                max_weight = weight;
                max_set = bv;
            }
        }
    } while (increment_bit_vector(bv));

    return max_set;
}

// Genetic Maxis Solver
// ====================

GeneticMaxisSolver::GeneticMaxisSolver(
            const Graph &g,
            genetic::Selector &sel,
            genetic::Recombinator &rec,
            genetic::Mutator &mut
        ) : constraint {std::numeric_limits<double>::infinity()},
            size{50},
            graph{g},
            selector{sel},
            recombinator{rec},
            mutator{mut} {

    // Reorder graph by sorting vertices by their weight divided by
    // their degree. Vertices with a large weight and few neighbors have
    // the highest probability of appearing in the solution.
    permutation.resize(graph.order());
    std::iota(begin(permutation), end(permutation), 0);
    auto &weights = graph.weights;
    auto &adj = graph.adjacency_list;
    std::sort(begin(permutation), end(permutation), [&weights, &adj](size_t i, size_t j) {
        return weights[i] / adj[i].size() > weights[j] / adj[j].size();
    });

    graph.reorder(permutation);
}

// The heuristic feasability operator ensures that genotypes created
// from mutation and breeding are valid independent sets. It is
// essentially a greedy solver and was described by Beasley and Chu. It
// makes a higher mutation rate necessary because some of the mutated
// bits will be immediately cancelled out to make the set independent.
//
// Currently it is the performance bottleneck for the genetic algorithm as it
// runs in O(|V|^2) time where |V| is the order of the graph.
void
GeneticMaxisSolver::heuristic_feasibility(size_t order, BitVector &adj, BitVector &chromosome) {
    // Traverse the chromosome, marking vertices which are most likely
    // to appear in an optimal solution as kept. Delete all neighbors of
    // the kept vertices from the chromosome.
    for (auto it = std::make_pair(begin(adj), begin(chromosome)); it.second != end(chromosome); it.first += order, ++it.second) {
        // If vertex is not selected, there is no conflict
        if (!*it.second) {
            continue;
        }

        // Delete all neighbors of the vertex. Because we are assuming
        // that the graph is bidirectional, we only need to traverse
        // half of the adjacency matrix.
        for (auto jt = std::make_pair(std::next(it.first, std::distance(begin(chromosome), it.second)), it.second);
                    jt.second != chromosome.end(); ++jt.first, ++jt.second) {

            if (*jt.first) {
                *jt.second = 0;
            }
        }
    }

    // Add back vertices where we can
    for (auto it = std::make_pair(begin(adj), begin(chromosome)); it.second != end(chromosome); it.first += order, ++it.second) {
        if (*it.second) {
            continue;
        }
        if(std::inner_product(begin(chromosome), end(chromosome), it.first, 0) == 0) {
            *it.second = 1;
        }
    }
}

// initialize_set randomly creates an independent set on the graph with
// adjacency matrix adj
void
GeneticMaxisSolver::initialize_set(size_t order, BitVector &adj, BitVector &chromosome) {
    static RNG rng;

    BitVector cover(order, 0);
    std::fill(begin(chromosome), end(chromosome), 0);
    auto uncovered_cnt = order;

    // Randomly select a set of vertices that completely cover the
    // graph
    while (uncovered_cnt > 0) {
        // Pick an uncovered vertex at random
        size_t skips = rng.random_index(0, uncovered_cnt - 1);
        size_t dist = 0;
        auto it = begin(cover);
        for (; it != end(cover); ++it)  {
            if (*it == 0 && skips-- == 0) {
                dist = std::distance(begin(cover), it);
                *it = 1;
                break;
            }
        }

        // Select it and mark all of its neighbors as covered
        chromosome[dist] = 1;
        --uncovered_cnt;
        for (auto i = begin(cover), j = begin(adj) + dist*order; i != end(cover); ++i, ++j) {
            if (*j != 0 && *i == 0) {
                *i = 1;
                if (--uncovered_cnt == 0) {
                    break;
                }
            }
        }
    }
}

void
GeneticMaxisSolver::iterate(AlgorithmState &state, size_t order,
        BitVector &adj, PopIter b, PopIter e, HashTable dupes) {

    // Select the weakest member of the population to be replaced
    dupes.erase(b->chromosome);

    // Select two parents for breeding and store the result into the
    // child, while ensuring that the newly created child is unique.
    // If the selector returns the same result every time, and the
    // mutator does not create a sufficiently unique solution, this
    // could loop forever.
    do {
        recombinator.breed(
            state,
            selector.select(state, b, e),
            selector.select(state, b, e),
            *b
        );
        mutator.mutate(state, *b);
        GeneticMaxisSolver::heuristic_feasibility(order, adj, *b->chromosome);

    } while (dupes.insert(b->chromosome) == false);

    // Calculate fitness of the new phenotype
    auto old_fitness = b->fitness;
    auto new_fitness = b->fitness = graph.weighted_total(*b->chromosome);

    // Use a single pass of bubble sort to put the new child in the
    // correct order. Phenotypes with the same fitness are ordered
    // by freshness.
    for (auto it = std::next(b); it != e && it->fitness <= new_fitness; ++it) {
        std::iter_swap(it, std::prev(it));
    }

    // Update state information
    state.total_fitness -= old_fitness;
    state.total_fitness += new_fitness;

    if (new_fitness > state.max_fitness) {
        state.max_fitness = new_fitness;
    }

    state.min_fitness = min_fitness;
    state.adjusted_fitness = state.total_fitness - state.min_fitness * state.size;
    ++state.iterations;
}

BitVector
permute_chromosome(const BitVector &c, vector<size_t> perm) {
    auto order = c.size();
    BitVector result(order);

    for (size_t i = 0; i < order; ++i) {
        result[perm[i]] = c[i];
    }

    return result;
}

void
initialize_algorithm_state(
        genetic::AlgorithmState &state,
        PopIter b,
        PopIter e) {

    state.min_fitness = b->fitness;
    state.max_fitness = std::prev(e)->fitness;
    state.total_fitness = std::accumulate(
        b, e, 0.0,
        [](double acc, genetic::Phenotype &ph) {
            return acc + ph.fitness;
        }
    );
    state.adjusted_fitness = state.total_fitness - state.min_fitness * state.size;
}

BitVector
GeneticMaxisSolver::operator()() {
    auto order = graph.order();
    auto adj = graph.adjacency_matrix();
    genetic::AlgorithmState state{size};

    // Keep a hash table to check for duplicate chromosomes
    BitVectorHashTable dupes{size};

    // Generate the initial population
    vector<BitVector> chromosomes;
    vector<genetic::Phenotype> pop;
    chromosomes.reserve(size);
    pop.reserve(size);

    // If the population size is larger than the set of possible
    // unique solutions, this will loop forever.
    for (decltype(size) i = 0; i < size; ++i) {
        chromosomes.emplace_back(order);
        pop.emplace_back(&chromosomes[i]);
        do {
            initialize_set(order, adj, chromosomes[i]);
        } while (dupes.insert(&chromosomes[i]) == false);

        pop[i].fitness = graph.weighted_total(*pop[i].chromosome);
    }

    // Sort population by fitness and initialize state information
    auto b = begin(pop), e = end(pop);
    initialize_algorithm_state(state, b, e);

    // Start handling SIGINT
    std::signal(SIGINT, sigint_handler);

    while (state.max_fitness < constraint && !sigint_flag) {
        iterate(state, order, adj, b, e, dupes) {
        if ((state.iterations & 0x3ff) == 0) {
            std::cout << "Iterations: " << state.iterations
                      << "; Average Fitness: " << state.total_fitness / size << std::endl;
        }
    }

    if (sigint_flag) {
        std::cout << "Interrupted..." << std::endl;
    }

    // Rearrange the most fit chromosome into the format of the
    // original, unsorted graph
    return permute_chromosome(*std::prev(end(pop))->chromosome, permutation);
}

// Parallel Genetic Solver
// =======================

WorkerSynchronizer::WorkerSynchronizer(
        unsigned int n
    ) : available_handles{n},
        current_cycle{0},
        utilized_handles{0},
        discarded_handles{0} {}

WorkerSynchronizer::Handle::Handle(WorkerSynchronizer &sync, unsigned int desired_cycle)
        : sync{sync} {

    std::unique_lock<std::mutex> l(sync.worker_lock);
    while(desired_cycle != sync.current_cycle) {
        sync.worker_cv.wait(l)
    }

    if (++sync.utilized_handles > sync.available_handles) {
        throw SynchronizationError("Too many workers attempting to reserve handles");
    }
}

WorkerSynchronizer::Handle::~Handle() {
    ++sync.discarded_handles;
    sync.migrator_cv.notify_all();
}

void
WorkerSynchronizer::manager_resume() {
    std::unique_lock<std::mutex> l(manager_lock);
    while(discarded_handles != available_handles) {
        manager_cv.wait(l);
    }
}

void
WorkerSynchronizer::manager_finish() {
    utilized_handles = 0;
    discarded_handles = 0;
    ++current_cycle;
    worker_cv.notify_all();
}

ParallelGeneticMaxisSolver::ParallelGeneticMaxisSolver(
            const Graph &g,
            genetic::AlgorithmStrategy &strat,
        ) : GeneticMaxisSolver(g, strat),
            migration_period{1000} {}

BitVector
ParallelGeneticMaxisSolver::operator()() {
    auto order = graph.order();
    auto adj = graph.adjacency_matrix();
    auto num_threads = std::thread::hardware_concurrency();

    // Ensure population size is evenly divisible by num_threads
    size += size % num_threads;

    // Initialize worker synchronizer
    WorkerSynchronizer sync{num_threads};

    // Initialize population
    vector<BitVector> chromosomes;
    vector<genetic::Phenotype> pop;
    chromosomes.reserve(size);
    pop.reserve(size);
    for (decltype(size) i = 0; i < size; ++i) {
        chromosomes.emplace_back(order);
        pop.emplace_back(&chromosomes[i]);
        GeneticMaxisSolver::initialize_set(order, adj, chromosomes[i]);
        pop[i].fitness = graph.weighted_total(*pop[i].chromosome);
    }

    // Assign subsets of the population to each thread and spawn them
    for (decltype(num_threads) i = 0; i < num_threads; ++i) {
        Worker w{
            graph,
            migration_period,
            b + i * size / num_threads,
            b + (i+1) * size / num_threads,
            strat, sync
        };
        std::thread t{w};
        t.detach();
    }

    // Wait for all worker threads to finish, then migrate phenotypes
    // between threads
    std::default_random_engine engine;
    while (1) {
        sync.wait_on_cycle();

        if (std::max_element(b, e)->fitness >= constraint) {
            break;
        }

        // Shuffle the population
        std::shuffle(begin(pop), end(pop), engine);

        sync.next_cycle();
    }

    return permute_chromosome(*std::max_element(b, e)->chromosome, permutation);
}

ParallelGeneticWorker::ParallelGeneticWorker(
        const Graph &g, unsigned int period,
        PopIter b, PopIter e,
        genetic::AlgorithmStrategy strat, WorkerSynchronizer &sync
    ) : graph{g}, migration_period{period}
        b{b}, e{e},
        strat{strat}, sync{sync} {}

void
ParallelGeneticWorker::operator()() {
    auto size = static_cast<size_t>(std::distance(b, e));
    genetic::AlgorithmState state{size};
    auto order = graph.order();
    auto adj = graph.adjacency_matrix();

    BitVectorHashTable dupes{size};

    for (unsigned int cycles = 0; true; ++cycles) {
        WorkerSynchronizer::Handle h{sync, cycles};

        // Ensure newly migrated population is unique
        dupes.clear();
        for (auto it = b; it != e; ++it) {
            auto is_changed = false;
            while(dupes.insert(it->chromosome) == false) {
                GeneticMaxisSolver::initialize_set(order, adj, *it->chromosome);
                is_changed = true;
            }
            if (is_changed) {
                it->fitness = graph.weighted_total(*it->chromosome);
            }
        }

        // Sort population and calculate fitness information
        std::sort(b, e);
        initialize_algorithm_state(state, b, e);

        for (decltype(migration_period) i = 0; i < migration_period; ++i) {
            GeneticMaxisSolver::iterate(state, strat, adj, b, e, dupes);
        }
    }
}


} // namespace maxis
