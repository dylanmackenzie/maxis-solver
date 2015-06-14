#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <limits>

#include "maxis/graph.hpp"
#include "maxis/genetic.hpp"
#include "maxis/solver.hpp"

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

void
heuristic_feasibility(const Graph &graph, genetic::Phenotype &ph) {
    using std::begin; using std::end;

    auto adj = graph.adjacency_matrix();
    auto order = graph.order();

    // Remove vertices until we have an independent set
    // it.first is the iterator over the adjacency matrix
    // it.second is the iterator over the set of vertices
    for (auto it = std::make_pair(begin(adj), begin(*ph.chromosome)); it.second != end(*ph.chromosome); it.first += order, ++it.second) {
        // If vertex is not selected, there is no conflict
        if (!*it.second) {
            continue;
        }

        // Delete all vertices that are in the set and neighbors of the
        // vertex which is currently being processed
        // TODO: change to std::rbegin() when gcc supports it
        auto cover = std::inner_product(begin(*ph.chromosome), end(*ph.chromosome), it.first, 0);
        decltype(adj)::reverse_iterator rit{it.first + order};
        for (auto jt = std::make_pair(rit, ph.chromosome->rbegin()); cover != 0 && jt.second != ph.chromosome->rend(); ++jt.first, ++jt.second) {
            if (*jt.first && *jt.second) {
                *jt.second = 0;
                --cover;
            }
        }
    }

    // Add back vertices to fill empty spaces
    for (auto it = std::make_pair(begin(adj), begin(*ph.chromosome)); it.second != end(*ph.chromosome); it.first += order, ++it.second) {
        if (*it.second) {
            continue;
        }
        if(std::inner_product(begin(*ph.chromosome), end(*ph.chromosome), it.first, 0) == 0) {
            *it.second = 1;
        }
    }
}

void
initialize_set(const Graph &g, genetic::Phenotype &ph) {
    using std::begin; using std::end;

    static genetic::RNG rng;

    auto order = g.order();
    auto adj = g.adjacency_matrix();

    BitVector cover(order, 0);
    std::fill(begin(*ph.chromosome), end(*ph.chromosome), 0);
    auto uncovered_cnt = order;

    // Randomly select a set of vertices that completely cover the
    // graph
    while (uncovered_cnt > 0) {
        // Pick an uncovered vertex at random
        size_t skips = rng.random_index(0, uncovered_cnt - 1);
        size_t dist;
        auto it = begin(cover);
        for (; it != end(cover); ++it)  {
            if (*it == 0 && skips-- == 0) {
                dist = std::distance(begin(cover), it);
                *it = 1;
                break;
            }
        }

        // Select it and mark all of its neighbors as covered
        (*ph.chromosome)[dist] = 1;
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

    size =  graph.order() / 2;

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

    // Keep a hash table to check for duplicate chromosomes
    auto hash_func = [](const Phenotype *ph) {
        return std::hash<std::vector<bool>>()(*ph->chromosome);
    };
    auto eq_func = [](const Phenotype *lhs, const Phenotype *rhs) {
        return std::equal(begin(*lhs->chromosome), end(*lhs->chromosome), begin(*rhs->chromosome));
    };
    std::unordered_set<Phenotype*, decltype(hash_func), decltype(eq_func)> dupes(size, hash_func, eq_func);

    // Generate the initial population
    std::vector<Phenotype> pop;
    pop.reserve(size);
    for (decltype(size) i = 0; i < size; ++i) {
        pop.emplace_back(order);
        do {
            initialize_set(graph, pop[i]);
        } while (!(dupes.insert(&pop[i]).second));

        pop[i].fitness = graph.weighted_total(*pop[i].chromosome);
    }

    // Calculate fitness information for the population
    auto min_fitness = std::min_element(begin(pop), end(pop))->fitness;
    auto max_fitness = std::max_element(begin(pop), end(pop))->fitness;
    std::cout << "Best independent set: " << max_fitness << std::endl;
    auto total_fitness = std::accumulate(
        begin(pop), end(pop), 0.0,
        [min_fitness](double acc, Phenotype &ph) {
            return acc + (ph.fitness - min_fitness);
        }
    );

    // TODO: parallelize
    while (max_fitness < constraint) {
        // Select the weakest member of the population to be replaced
        auto &child = *std::find_if(
            begin(pop), end(pop),
            [min_fitness](const Phenotype &ph) { return ph.fitness == min_fitness; }
        );

        // Erase it and update the aggregate fitness information
        dupes.erase(&child);

        do {
            // Select two parents for breeding
            auto &p1 = selector.select(begin(pop), end(pop), min_fitness, total_fitness);
            auto &p2 = selector.select(begin(pop), end(pop), min_fitness, total_fitness);

            // Select two parents for breeding and store the result into the
            // child, while ensuring that the newly created child is unique
            recombinator.crossover(p1, p2, child);
            mutator.mutate(child);
            heuristic_feasibility(graph, child);

        } while (!(dupes.insert(&child).second));

        // Compute fitness of new child
        child.fitness = graph.weighted_total(*child.chromosome);
        if (child.fitness <= min_fitness) {
            min_fitness = child.fitness;
        } else {
            min_fitness = std::min_element(begin(pop), end(pop))->fitness;
        }

        if (child.fitness > max_fitness) {
            max_fitness = child.fitness;
            std::cout << "Best independent set: " << max_fitness << std::endl;
        }

        total_fitness = std::accumulate(
            begin(pop), end(pop), 0.0,
            [min_fitness](double acc, Phenotype &ph) {
                return acc + (ph.fitness - min_fitness);
            }
        );
    }

    auto max = std::max_element(begin(pop), end(pop));

    BitVector result(order);
    for (size_t i = 0; i < order; ++i) {
        result[permutation[i]] = (*max->chromosome)[i];
    }
    return result;
}

} // namespace maxis
