#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <limits>

#include "maxis/graph.hpp"
#include "maxis/genetic.hpp"
#include "maxis/solver.hpp"
#include "maxis/rng.hpp"

namespace maxis {

// Brute Force Solver

BruteForceMaxisSolver::BruteForceMaxisSolver(const Graph &g) : graph{g} {}

// Ripple carry increment for a collection of bits.
// Returns true on overflow.
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

    return carry;
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
    } while (increment_bit_vector(bv) == false);

    return max_set;
}

// Genetic Maxis Solver

// The heuristic feasability operator ensures that genotypes created
// from mutation and breeding are valid independent sets
// TODO: use decltype as adj need not be a BitVector
void
heuristic_feasibility(size_t order, BitVector &adj, BitVector &chromosome) {
    using std::begin; using std::end;

    // Remove vertices until we have an independent set
    // it.first is the iterator over the adjacency matrix
    // it.second is the iterator over the set of vertices
    for (auto it = std::make_pair(begin(adj), begin(chromosome)); it.second != end(chromosome); it.first += order, ++it.second) {
        // If vertex is not selected, there is no conflict
        if (!*it.second) {
            continue;
        }

        // Delete all vertices that are in the set and neighbors of the
        // vertex which is currently being processed
        // TODO: change to std::rbegin() when gcc supports it
        auto cover = std::inner_product(begin(chromosome), end(chromosome), it.first, 0);
        BitVector::reverse_iterator rit{std::next(it.first, order)};
        for (auto jt = std::make_pair(rit, chromosome.rbegin()); cover != 0 && jt.second != chromosome.rend(); ++jt.first, ++jt.second) {
            if (*jt.first && *jt.second) {
                *jt.second = 0;
                --cover;
            }
        }
    }

    // Add back vertices to fill empty spaces
    for (auto it = std::make_pair(begin(adj), begin(chromosome)); it.second != end(chromosome); it.first += order, ++it.second) {
        if (*it.second) {
            continue;
        }
        if(std::inner_product(begin(chromosome), end(chromosome), it.first, 0) == 0) {
            *it.second = 1;
        }
    }
}

// initialize_set is used to generate the initial population with valid
// independent sets on the graph
void
initialize_set(size_t order, BitVector &adj, BitVector &chromosome) {
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

GeneticMaxisSolver::GeneticMaxisSolver(
            const Graph &g,
            genetic::Selector &sel,
            genetic::Recombinator &rec,
            genetic::Mutator &mut
        ) : graph{g},
            selector{sel},
            recombinator{rec},
            mutator{mut}
        {

    using std::begin; using std::end;

    size = graph.order() / 2;

    // Reorder graph by sorting vertices by (weight / degree)
    permutation.resize(graph.order());
    std::iota(begin(permutation), end(permutation), 0);
    auto &weights = graph.weights;
    auto &adj = graph.adjacency_list;
    std::sort(begin(permutation), end(permutation), [&weights, &adj](size_t i, size_t j) {
        return weights[i] / adj[i].size() > weights[j] / adj[j].size();
    });

    graph.reorder(permutation);
}

BitVector
GeneticMaxisSolver::operator()() {
    using std::begin; using std::end;
    using genetic::Phenotype;

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
    std::vector<Phenotype> pop;
    pop.reserve(size);
    chromosomes.reserve(size);
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
        [](double acc, Phenotype &ph) {
            return acc + ph.fitness;
        }
    );
    state.adjusted_fitness = state.total_fitness - state.min_fitness * size;

    // TODO: parallelize
    while (state.max_fitness < constraint) {

        // Select the weakest member of the population to be replaced
        auto &child = *begin(pop);
        dupes.erase(child.chromosome);

        do {
            // Select two parents for breeding and store the result into the
            // child, while ensuring that the newly created child is unique
            recombinator.breed(
                state,
                selector.select(state, begin(pop), end(pop)),
                selector.select(state, begin(pop), end(pop)),
                child
            );
            mutator.mutate(state, child);
            heuristic_feasibility(order, adj, *child.chromosome);

        } while (dupes.insert(child.chromosome).second == false);

        // Calculate fitness of new population member and update the
        // total_fitness
        state.total_fitness -= child.fitness;
        child.fitness = graph.weighted_total(*child.chromosome);
        state.total_fitness += child.fitness;

        // Log fitness info
        if (child.fitness > state.max_fitness) {
            state.max_fitness = child.fitness;
            std::cout << "Best independent set (" << state.max_fitness
                      << "): "  << std::endl;
        }

        // Insert new child into sorted position and update algorithm state
        auto sorted_pos = std::lower_bound(std::next(begin(pop)), end(pop), child);
        for (auto it = begin(pop); it != std::prev(sorted_pos); ++it) {
            std::swap(*it, *std::next(it));
        }
        state.min_fitness = begin(pop)->fitness;
        state.adjusted_fitness = state.total_fitness - state.min_fitness * size;
        ++state.iterations;
    }

    auto max = std::prev(end(pop));

    BitVector result(order);
    for (size_t i = 0; i < order; ++i) {
        result[permutation[i]] = (*max->chromosome)[i];
    }
    return result;
}

} // namespace maxis
