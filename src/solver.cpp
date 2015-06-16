#include <algorithm>
#include <csignal>
#include <iostream>
#include <limits>
#include <unordered_set>

#include "maxis/graph.hpp"
#include "maxis/genetic.hpp"
#include "maxis/solver.hpp"
#include "maxis/rng.hpp"

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
    using std::begin; using std::end;

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

    using std::begin; using std::end;

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
// Currently it is the bottleneck for a single iteration as it
// runs in O(|V|^2) time where |V| is the order of the graph.
void
GeneticMaxisSolver::heuristic_feasibility(size_t order, BitVector &adj, BitVector &chromosome) {
    using std::begin; using std::end;

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
    using std::begin; using std::end;

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

BitVector
GeneticMaxisSolver::operator()() {
    using std::begin; using std::end;

    auto order = graph.order();
    auto adj = graph.adjacency_matrix();
    genetic::AlgorithmState state;

    // Keep a hash table to check for duplicate chromosomes
    auto hash_func = [](BitVector *chromosome) {
        return std::hash<std::vector<bool>>()(*chromosome);
    };
    auto eq_func = [](BitVector *lhs, BitVector *rhs) {
        return std::equal(begin(*lhs), end(*lhs), begin(*rhs));
    };
    std::unordered_set<BitVector*, decltype(hash_func), decltype(eq_func)> dupes(size, hash_func, eq_func);

    // Generate the initial population
    std::vector<BitVector> chromosomes;
    std::vector<genetic::Phenotype> pop;
    pop.reserve(size);
    chromosomes.reserve(size);

    // If the population size is larger than the set of possible
    // unique solutions, this will loop forever.
    for (decltype(size) i = 0; i < size; ++i) {
        chromosomes.emplace_back(order);
        pop.emplace_back(&chromosomes[i]);
        do {
            initialize_set(order, adj, chromosomes[i]);
        } while (dupes.insert(&chromosomes[i]).second == false);

        pop[i].fitness = graph.weighted_total(*pop[i].chromosome);
    }

    // Sort population by fitness and initialize state information
    std::sort(begin(pop), end(pop));
    state.min_fitness = begin(pop)->fitness;
    state.max_fitness = end(pop)->fitness;
    state.total_fitness = std::accumulate(
        begin(pop), end(pop), 0.0,
        [](double acc, genetic::Phenotype &ph) {
            return acc + ph.fitness;
        }
    );
    state.adjusted_fitness = state.total_fitness - state.min_fitness * size;

    // Start handling SIGINT
    std::signal(SIGINT, sigint_handler);

    while (state.max_fitness < constraint && !sigint_flag) {

        // Select the weakest member of the population to be replaced
        auto &child = *begin(pop);
        dupes.erase(child.chromosome);

        // Select two parents for breeding and store the result into the
        // child, while ensuring that the newly created child is unique.
        // If the selector returns the same result every time, and the
        // mutator does create a sufficiently unique solution, this
        // could loop forever.
        do {
            recombinator.breed(
                state,
                selector.select(state, begin(pop), end(pop)),
                selector.select(state, begin(pop), end(pop)),
                child
            );
            mutator.mutate(state, child);
            heuristic_feasibility(order, adj, *child.chromosome);

        } while (dupes.insert(child.chromosome).second == false);

        // Calculate fitness of the new phenotype and update the
        // total fitness
        state.total_fitness -= child.fitness;
        auto child_fitness = child.fitness = graph.weighted_total(*child.chromosome);
        state.total_fitness += child_fitness;

        // Log whenever we have an increase in the maximum fitness
        if (child_fitness > state.max_fitness) {
            state.max_fitness = child_fitness;
            std::cout << "Best independent set: " << state.max_fitness << std::endl;
        }

        // Use a single pass of bubble sort to put the new child in the
        // correct order. Phenotypes with the same fitness are ordered
        // by freshness.
        for (auto it = std::next(begin(pop)); it != end(pop) && it->fitness <= child_fitness; ++it) {
            std::swap(*it, *std::prev(it));
        }

        // Update the rest of the state information
        state.min_fitness = begin(pop)->fitness;
        state.adjusted_fitness = state.total_fitness - state.min_fitness * size;
        ++state.iterations;
        if ((state.iterations & 0x3ff) == 0) {
            std::cout << "Iterations: " << state.iterations
                      << "; Average Fitness: " << state.total_fitness / size << std::endl;
        }
    }

    auto max = std::prev(end(pop));

    if (sigint_flag) {
        std::cout << "Interrupted..." << std::endl;
    }

    // Rearrange the chromosome into the format of the original,
    // unsorted graph
    BitVector result(order);
    for (size_t i = 0; i < order; ++i) {
        result[permutation[i]] = (*max->chromosome)[i];
    }
    return result;
}

} // namespace maxis
