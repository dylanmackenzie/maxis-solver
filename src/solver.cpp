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

GeneticMaxisSolver::GeneticMaxisSolver(
            const Graph &g,
            genetic::MaxisHeuristicGenerator &gen,
            genetic::Selector &sel,
            genetic::Recombinator &rec,
            genetic::Mutator &mut
        ) : graph{g},
            generator{gen},
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

    // TODO: this is a hack for the Maxis specific functions. We need
    // to make this function a template.
    gen.graph = &graph;

    graph.reorder(permutation);
}

// Fitness for a given selection of vertices is just the inner
// product of the weights and a bit vector where a 1 represents a
// selected vertex and a 0 represents an unselected one.
void
GeneticMaxisSolver::compute_fitness(genetic::Phenotype& ph) const {
    ph.fitness = std::inner_product(
        std::begin(graph.weights),
        std::end(graph.weights),
        std::begin(*ph.chromosome),
        0.0
    );
}

BitVector
GeneticMaxisSolver::operator()() {
    using std::begin; using std::end;
    using genetic::Phenotype;

    genetic::MaxisHeuristicMutator heuristic{};
    heuristic.graph = &graph;

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
            generator.initialize(pop[i]);
        } while (!(dupes.insert(&pop[i]).second));

        compute_fitness(pop[i]);
    }

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
            heuristic.mutate(child);

        } while (!(dupes.insert(&child).second));

        // Compute fitness of new child
        compute_fitness(child);
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
