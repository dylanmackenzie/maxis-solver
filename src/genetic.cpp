#include <algorithm>
#include <vector>

#include "maxis/genetic.hpp"

namespace genetic {

// Selectors

// TODO: implement a tournament selector
Phenotype&
TournamentSelector::select(const AlgorithmState &state, std::vector<Phenotype>::iterator b, std::vector<Phenotype>::iterator e) {
    using std::begin; using std::end;

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
    using std::begin; using std::end;

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

// TODO: ensure probabilities are accurate
void
BlendingRecombinator::breed(const AlgorithmState &state, const Phenotype &p1, const Phenotype &p2, Phenotype &c) {
    using std::begin; using std::end;

    auto f1 = p1.fitness;
    auto f2 = p2.fitness;

    for (auto i = begin(*p1.chromosome), j = begin(*p2.chromosome), k = begin(*c.chromosome);
            i != end(*p1.chromosome); ++i, ++j, ++k) {

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

void
SimpleMutator::mutate(const AlgorithmState &state, Phenotype &ph) {
    using std::begin; using std::end;

    double rate = probability_multiplier / ph.chromosome->size();

    for (auto it = begin(*ph.chromosome); it != end(*ph.chromosome); ++it) {
        if (rng.random_prob() < rate) {
            *it = !*it;
        }
    }
}

void
VariableRateMutator::mutate(const AlgorithmState &state, Phenotype &ph) {
    using std::begin; using std::end;

    double rate = exp(-4 * gradient * (state.iterations - halfway) / steady_state);
    rate = steady_state / (1 + rate);
    rate /= ph.chromosome->size();

    for (auto it = begin(*ph.chromosome); it != end(*ph.chromosome); ++it) {
        if (rng.random_prob() < rate) {
            *it = !*it;
        }
    }
}

} // namespace genetic
