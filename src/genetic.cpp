#include <algorithm>
#include <vector>
#include <unordered_set>
#include <cassert>

#include "maxis/genetic.hpp"

namespace maxis {
namespace genetic {

// Random Number Generator
double
RNG::random_prob() {
    return probability(engine);
}

size_t
RNG::random_index(size_t start, size_t end) {
    std::uniform_int_distribution<size_t> range{start, end};
    return range(engine);
}

// Generators
MaxisHeuristicGenerator::MaxisHeuristicGenerator() {}

void
MaxisHeuristicGenerator::initialize(Phenotype &ph) {
    using std::begin; using std::end;

    auto order = graph->order();
    auto adj = graph->adjacency_matrix();

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

// Selectors

// TODO: implement a tournament selector
Phenotype&
TournamentSelector::select(std::vector<Phenotype>::iterator b, std::vector<Phenotype>::iterator e, double min, double total) {
    using std::begin; using std::end;

    Phenotype *selection = &*(b + rng.random_index(0, std::distance(b, e) - 1));
    double best_fitness = selection->fitness;
    for (size_t i = 1; i < size; ++i) {
        auto idx = rng.random_index(0, std::distance(b, e) - 1);
        if ((b + idx)->fitness > best_fitness) {
            selection = &*(b + idx);
        }

    }

    return *selection;
}

Phenotype&
RouletteSelector::select(std::vector<Phenotype>::iterator b, std::vector<Phenotype>::iterator e, double min, double total) {
    using std::begin; using std::end;

    // TODO: ensure correctness in face of floating point rounding
    auto rn = rng.random_prob();
    double sum = 0;
    for (auto it = b; it != e; ++it) {
        sum += (it->fitness - min) / total;
        if (sum >= rn) {
            return *it;
        }
    }

    // Loop should always return
    return *(e - 1);
}

// Recombinators

// void
// SinglePointRecombinator::crossover(Phenotype &p1, Phenotype &p2, Phenotype &c1, Phenotype &c2) {
    // using std::begin; using std::end;

    // size_t cp = rng.random_index(0, (*p1.chromosome).size());
    // std::copy(begin(*p1.chromosome), begin(*p1.chromosome) + cp, begin(c1.chromosome));
    // std::copy(begin(*p2.chromosome), begin(*p2.chromosome) + cp, begin(c2.chromosome));
    // std::copy(begin(*p1.chromosome) + cp, end(*p1.chromosome), begin(c2.chromosome) + cp);
    // std::copy(begin(*p2.chromosome) + cp, end(*p2.chromosome), begin(c1.chromosome) + cp);
// }

// void
// DoublePointRecombinator::crossover(Phenotype &p1, Phenotype &p2, Phenotype &c1, Phenotype &c2) {
    // using std::begin; using std::end;

    // size_t size = (*p1.chromosome).size();
    // size_t cp1 = rng.random_index(0, size);
    // size_t cp2 = rng.random_index(0, size);
    // if (cp1 > cp2) {
        // std::swap(cp1, cp2);
    // }

    // size_t cp = rng.random_index(0, (*p1.chromosome).size());
    // std::copy(begin(*p1.chromosome), begin(*p1.chromosome) + cp1, begin(c1.chromosome));
    // std::copy(begin(*p2.chromosome), begin(*p2.chromosome) + cp1, begin(c2.chromosome));
    // std::copy(begin(*p1.chromosome) + cp1, begin(*p1.chromosome) + cp2, begin(c2.chromosome) + cp1);
    // std::copy(begin(*p2.chromosome) + cp1, begin(*p2.chromosome) + cp2, begin(c1.chromosome) + cp1);
    // std::copy(begin(*p1.chromosome) + cp2, end(*p1.chromosome), begin(c1.chromosome) + cp2);
    // std::copy(begin(*p2.chromosome) + cp2, end(*p2.chromosome), begin(c2.chromosome) + cp2);

// }

// TODO: ensure probabilities are accurate
void
BlendingRecombinator::crossover(Phenotype &p1, Phenotype &p2, Phenotype &c) {
    using std::begin; using std::end;

    auto f1 = p1.fitness;
    auto f2 = p2.fitness;

    for (auto i = begin(*p1.chromosome), j = begin(*p2.chromosome), k = begin(*c.chromosome);
            i != end(*p1.chromosome); ++i, ++j, ++k) {

        if (*i == *j) {
            *k = *i;
        }

        double prob = f1 / (f1 + f2);

        // Inherit bit proportional to fitness
        if (rng.random_prob() < prob) {
            *k = *i;
        } else {
            *k = *j;
        }
    }
}

// Mutators

SinglePointMutator::SinglePointMutator(double mul) : probability_multiplier{mul} {}

void
SinglePointMutator::mutate(Phenotype &ph) {
    using std::begin; using std::end;

    auto l = ph.chromosome->size();
    double prob = probability_multiplier / l;

    for (auto it = begin(*ph.chromosome); it != end(*ph.chromosome); ++it) {
        if (rng.random_prob() < prob) {
            *it = !*it;
        }
    }
}

MaxisHeuristicMutator::MaxisHeuristicMutator() {}

void
MaxisHeuristicMutator::mutate(Phenotype &ph) {
    using std::begin; using std::end;

    auto adj = graph->adjacency_matrix();
    auto order = graph->order();

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

    assert(graph->is_independent_set(*ph.chromosome));
}

} // namespace genetic
} // namespace maxis
