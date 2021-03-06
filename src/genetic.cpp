#include <algorithm>
#include <iostream>
#include <vector>

#include "maxis/bit_vector.hpp"
#include "maxis/genetic.hpp"

namespace maxis {
namespace genetic {

using std::begin; using std::end;
using std::cbegin; using std::cend;
using maxis::begin; using maxis::end;
using maxis::cbegin; using maxis::cend;

AlgorithmState::AlgorithmState() :
        min_fitness{0}, max_fitness{0}, total_fitness{0},
        adjusted_fitness{0}, iterations{0}, size{0} {};

// initialize assumes the list is in sorted order
void
AlgorithmState::initialize(std::vector<Phenotype>::iterator b, std::vector<Phenotype>::iterator e) {
    size = static_cast<size_t>(std::distance(b, e));
    min_fitness = b->fitness;
    max_fitness = std::prev(e)->fitness;
    total_fitness = std::accumulate(
        b, e, 0.0,
        [](double acc, genetic::Phenotype &ph) {
            return acc + ph.fitness;
        }
    );
    adjusted_fitness = total_fitness - min_fitness * size;
}

std::ostream&
operator<<(std::ostream& os, const AlgorithmState& state) {
    os << state.iterations << ": " << state.max_fitness;
    return os;
}


// Selectors
// =========

Phenotype&
TournamentSelector::select(const AlgorithmState &state, std::vector<Phenotype>::iterator b, std::vector<Phenotype>::iterator e) {

    auto selection = std::next(b, rng.random_index(0, std::distance(b, e) - 1));
    double best_fitness = selection->fitness;
    for (size_t i = 1; i < size; ++i) {
        auto n = rng.random_index(0, std::distance(b, e) - 1);
        if (std::next(b, n)->fitness > best_fitness) {
            selection = std::next(b, n);
            best_fitness = selection->fitness;
        }

    }

    return *selection;
}

Phenotype&
RouletteSelector::select(const AlgorithmState &state, std::vector<Phenotype>::iterator b, std::vector<Phenotype>::iterator e) {
    // TODO: ensure correctness in face of floating point rounding
    auto rn = rng.random_prob();
    double sum = 0;
    for (auto it = b; it != e; ++it) {
        sum += (it->fitness - state.min_fitness) / state.adjusted_fitness;
        if (sum >= rn) {
            return *it;
        }
    }

    // select should always return from loop
    return *(e - 1);
}

// Recombinators
// =============

void
BlendingRecombinator::breed(const AlgorithmState &state, const Phenotype &p1, const Phenotype &p2, Phenotype &c) {
    auto f1 = p1.fitness;
    auto f2 = p2.fitness;

    auto k = begin(*c.chromosome);
    for (auto i = cbegin(*p1.chromosome), j = cbegin(*p2.chromosome); i != cend(*p1.chromosome); ++i, ++j, ++k) {
        if (*i == *j) {
            *k = *i;
        }

        double prob = f1 / (f1 + f2);

        // Inherit bit with probability proportional to fitness
        if (rng.random_prob() < prob) {
            *k = *i;
        } else {
            *k = *j;
        }
    }
}

// Mutators
// ========

void
SimpleMutator::mutate(const AlgorithmState &state, Phenotype &ph) {
    double rate = probability_multiplier / ph.chromosome->size();

    for (auto it = begin(*ph.chromosome); it != end(*ph.chromosome); ++it) {
        if (rng.random_prob() < rate) {
            *it = !*it;
        }
    }
}

void
VariableRateMutator::mutate(const AlgorithmState &state, Phenotype &ph) {
    double rate = exp(-4 * gradient * (state.iterations - halfway) / steady_state);
    rate = steady_state / (1 + rate);
    double prob = rate / ph.chromosome->size();

    for (auto it = begin(*ph.chromosome); it != end(*ph.chromosome); ++it) {
        if (rng.random_prob() < prob) {
            *it = !*it;
        }
    }
}

} // namespace genetic
} // namespace maxis
