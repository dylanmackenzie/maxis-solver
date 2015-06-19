#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <csignal>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <thread>
#include <unordered_set>
#include <vector>

#include "maxis/graph.hpp"
#include "maxis/genetic.hpp"
#include "maxis/solver.hpp"
#include "maxis/synchronization.hpp"
#include "maxis/rng.hpp"

using std::vector;
using std::begin; using std::end;
using namespace maxis::genetic;
using PopIter = std::vector<Phenotype>::iterator;

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
            genetic::AlgorithmStrategy strat
        ) : constraint {std::numeric_limits<double>::infinity()},
            size{50}, graph{g}, strategy{strat} {

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
GeneticMaxisSolver::heuristic_feasibility(const Graph &graph, BitVector &chromosome) {
    auto order = graph.order();
    auto adj = graph.adjacency_matrix();

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
GeneticMaxisSolver::initialize_set(const Graph &graph, BitVector &chromosome) {
    static thread_local RNG rng;

    auto order = graph.order();
    auto adj = graph.adjacency_matrix();

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
GeneticMaxisSolver::iterate(const Graph &graph, PopIter b, PopIter e,
        AlgorithmState &state, AlgorithmStrategy &strat, BitVectorHashTable<BitVector> &dupes) {

    // Select the weakest member of the population to be replaced
    dupes.erase(b->chromosome);

    // Select two parents for breeding and store the result into the
    // child, while ensuring that the newly created child is unique.
    // If the selector returns the same result every time, and the
    // mutator does not create a sufficiently unique solution, this
    // could loop forever.
    do {
        strat.recombinator.breed(
            state,
            strat.selector.select(state, b, e),
            strat.selector.select(state, b, e),
            *b
        );
        strat.mutator.mutate(state, *b);
        GeneticMaxisSolver::heuristic_feasibility(graph, *b->chromosome);

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

    state.min_fitness = b->fitness;
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
debug_heading() {
    using std::cout; using std::endl; using std::setw;
    const int cw = 15;

    cout << setw(cw) << "Iteration"
         << setw(8) << "Thread"
         << setw(cw) << "Maximum"
         << setw(cw) << "Average\n";

    cout << std::string(cw * 4, '=') << endl;
}

void
debug_iteration(const genetic::AlgorithmState &state, unsigned int n) {
    using std::cout; using std::endl; using std::setw;
    const int cw = 15;

    if (n == 0) {
        cout << setw(cw) << state.iterations;
    } else {
        cout << setw(cw) << "";
    }

    cout << setw(6) << n
         << setw(cw) << state.max_fitness
         << setw(cw) << state.total_fitness / state.size
         << endl;
}

BitVector
GeneticMaxisSolver::operator()() {
    auto order = graph.order();
    auto adj = graph.adjacency_matrix();
    genetic::AlgorithmState state{};

    // Keep a hash table to check for duplicate chromosomes
    BitVectorHashTable<BitVector> dupes{size};

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
            initialize_set(graph, chromosomes[i]);
        } while (dupes.insert(&chromosomes[i]) == false);

        pop[i].fitness = graph.weighted_total(*pop[i].chromosome);
    }

    // Sort population by fitness and initialize state information
    auto b = begin(pop), e = end(pop);
    std::sort(b, e);
    state.initialize(b, e);

    // Start handling SIGINT and begin logging to stdout
    std::signal(SIGINT, sigint_handler);

#ifdef DEBUG
    debug_heading();
#endif

    while (state.max_fitness < constraint && !sigint_flag) {
        // Iterate the genetic algorithm
        iterate(graph, b, e, state, strategy, dupes);

#ifdef DEBUG
        // Log every 1024 iterations
        if ((state.iterations & 0x3ff) == 0) {
            debug_iteration(state, 0);
        }
#endif
    }

    if (sigint_flag) {
        std::cout << "\nInterrupted..." << std::endl;
    }

    // Rearrange the most fit chromosome into the format of the
    // original, unsorted graph
    return permute_chromosome(*std::prev(end(pop))->chromosome, permutation);
}

// Parallel Genetic Solver
// =======================

ParallelGeneticMaxisSolver::ParallelGeneticMaxisSolver(
            const Graph &g,
            genetic::AlgorithmStrategy strat
        ) : GeneticMaxisSolver(g, strat),
            migration_period{512} {}

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
    vector<genetic::AlgorithmState> worker_states(num_threads);
    chromosomes.reserve(size);
    pop.reserve(size);
    for (decltype(size) i = 0; i < size; ++i) {
        chromosomes.emplace_back(order);
        pop.emplace_back(&chromosomes[i]);
        GeneticMaxisSolver::initialize_set(graph, chromosomes[i]);
        pop[i].fitness = graph.weighted_total(*pop[i].chromosome);
    }

    // Assign subsets of the population to each thread and spawn them
    auto b = begin(pop), e = end(pop);
    for (decltype(num_threads) i = 0; i < num_threads; ++i) {
        ParallelGeneticWorker w{
            graph,
            migration_period,
            std::next(b, i * size / num_threads),
            std::next(b, (i+1) * size / num_threads),
            worker_states[i],
            strategy, sync
        };
        std::thread t{w};
        t.detach();
    }

    // Start handling SIGINT
    std::signal(SIGINT, sigint_handler);

#ifdef DEBUG
    debug_heading();
#endif

    // Wait for all worker threads to finish, then migrate phenotypes
    // between threads
    std::default_random_engine engine;
    while (1) {
        sync.wait_on_cycle();

        if (std::max_element(b, e)->fitness >= constraint || sigint_flag) {
            break;
        }

#ifdef DEBUG
        for (decltype(num_threads) i = 0; i < num_threads; ++i) {
            debug_iteration(worker_states[i], i);
        }
#endif

        // Shuffle the population
        std::shuffle(b, e, engine);

        sync.next_cycle();
    }

    if (sigint_flag) {
        std::cout << "\nInterrupted..." << std::endl;
    }

    auto max = std::max_element(b, e)->chromosome;
    return permute_chromosome(*max, permutation);
}

ParallelGeneticWorker::ParallelGeneticWorker(
        const Graph &g, unsigned int period,
        PopIter b, PopIter e, genetic::AlgorithmState &state,
        genetic::AlgorithmStrategy strat, WorkerSynchronizer &sync
    ) : graph{g}, migration_period{period},
        b{b}, e{e}, state{state},
        strategy{strat}, sync{sync} {}

void
ParallelGeneticWorker::operator()() {
    auto size = static_cast<size_t>(std::distance(b, e));

    BitVectorHashTable<BitVector> dupes{size};

    for (unsigned int cycles = 0; true; ++cycles) {
        WorkerSynchronizer::Handle h{sync, cycles};

        // Ensure newly migrated population is unique
        dupes.clear();
        for (auto it = b; it != e; ++it) {
            auto is_changed = false;
            while(dupes.insert(it->chromosome) == false) {
                GeneticMaxisSolver::initialize_set(graph, *it->chromosome);
                is_changed = true;
            }
            if (is_changed) {
                it->fitness = graph.weighted_total(*it->chromosome);
            }
        }

        // Sort population and calculate fitness information
        std::sort(b, e);
        state.initialize(b, e);

        for (decltype(migration_period) i = 0; i < migration_period; ++i) {
            GeneticMaxisSolver::iterate(graph, b, e, state, strategy, dupes);
        }
    }
}


} // namespace maxis
