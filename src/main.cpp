#include <iostream>
#include <fstream>
#include <sstream>

#include "maxis/graph.hpp"
#include "maxis/solver.hpp"
#include "maxis/genetic.hpp"

maxis::BitVector
solve_maxis(maxis::Graph &g, double constraint) {
    using namespace maxis::genetic;

    MaxisHeuristicGenerator gen{};
    TournamentSelector sel{2};
    BlendingRecombinator rec{};
    SinglePointMutator mut{4};

    maxis::GeneticMaxisSolver solver{g, gen, sel, rec, mut};
    solver.constraint = constraint;

    return solver();
}

int
main(int argc, char *argv[]) {
    using std::cout;
    using std::endl;

    if (argc < 2) {
        cout << "Need input file" << endl;
        return 1;
    }
    std::ifstream is(argv[1]);

    auto test = maxis::Graph::from_ascii_dimacs(is);
    auto coloring = solve_maxis(test, std::stoi(argv[2]));

    cout << endl << "Results" << endl << "=======" << endl;
    cout << test.is_independent_set(coloring) << endl;
    cout << "Vertices: ";

    for (size_t i = 0; i < test.order(); ++i) {
        if (coloring[i]) {
            cout << i << " ";
        }
    }

    cout << endl;
    cout << "Weight: " << test.weighted_total(coloring) << endl;
}
