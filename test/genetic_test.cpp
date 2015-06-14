#include <algorithm>
#include <fstream>
#include <gtest/gtest.h>

#include "maxis/graph.hpp"
#include "maxis/genetic.hpp"

using namespace maxis;

class GeneticTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        std::ifstream naaru_stream("test/graphs/naaru24.mis");
        std::ifstream wiki_stream("test/graphs/wiki6.mis");
        naaru = Graph::from_ascii_dimacs(naaru_stream);
        wiki = Graph::from_ascii_dimacs(wiki_stream);
    }

    Graph naaru;
    Graph wiki;
};

TEST_F(GeneticTest, MaxisHeuristicGenerator) {
    using namespace genetic;

    auto gen = MaxisHeuristicGenerator{};
    gen.graph = &naaru;
    auto pop = gen.initial_population(100);

    for (auto &ph : pop) {
        EXPECT_TRUE(naaru.is_independent_set(*ph.chromosome));
    }
}

TEST_F(GeneticTest, MaxisHeuristicMutator) {
    using namespace genetic;
    auto mut = MaxisHeuristicMutator{};
    mut.graph = &naaru;
    Phenotype pn{naaru.order()};
    std::fill(std::begin(*pn.chromosome), std::end(*pn.chromosome), 1);
    mut.mutate(pn);
    EXPECT_TRUE(naaru.is_independent_set(*pn.chromosome));

    mut.graph = &wiki;
    Phenotype pw{wiki.order()};
    std::fill(std::begin(*pw.chromosome), std::end(*pw.chromosome), 1);
    mut.mutate(pw);
    EXPECT_TRUE(wiki.is_independent_set(*pw.chromosome));
}
