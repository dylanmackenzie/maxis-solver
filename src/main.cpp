#include <iostream>
#include <fstream>
#include <sstream>

#include "boost/program_options.hpp"

#include "maxis/graph.hpp"
#include "maxis/solver.hpp"
#include "maxis/genetic.hpp"

using namespace maxis;

int
main(int argc, char *argv[]) {
    namespace po = boost::program_options;
    using std::cout;
    using std::endl;

    // Parse arguments
    po::options_description opts{"Options"};
    opts.add_options()
        ("help,h", "Display usage information")
        ("file", po::value<std::string>(), "Input file")
        ("population,p", po::value<size_t>(), "Population size")
        ("target,w", po::value<double>(), "Target fitness value")
        ("mutation,m", po::value<double>()->default_value(9), "Mutation rate")
        ("mutation-start", po::value<size_t>()->default_value(5000),
             "(Variable Rate Mutator) Point at which mutation reaches half of its final rate")
        ("mutation-gradient", po::value<double>()->default_value(0.05, "0.05"),
             "(Variable Rate Mutator) Mutation growth rate")
        ("selector,s", po::value<size_t>()->default_value(2),
             "(Tournament Selector) Tournament Size");

    po::positional_options_description p;
    p.add("file", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(opts).positional(p).run(), vm);
    po::notify(vm);
    if (vm.count("help")) {
        cout << opts << endl;
        return 0;
    }
    if (!vm.count("file")) {
        cout << "No input file specified" << endl;
        return 1;
    }

    // Parse and initialize graph
    std::ifstream is(vm["file"].as<std::string>());
    auto graph = maxis::Graph::from_ascii_dimacs(is);

    // Initialize genetic solver
    auto pop_size = graph.order() / 2;
    if (vm.count("population")) {
        pop_size = vm["population"].as<size_t>();
    }

    // genetic::RouletteSelector sel{};
    genetic::TournamentSelector sel{vm["selector"].as<size_t>()};

    genetic::BlendingRecombinator rec{};

    // genetic::SimpleMutator mut{vm["mutation"].as<double>()};
    genetic::VariableRateMutator mut{
        vm["mutation"].as<double>(),
        vm["mutation-start"].as<size_t>(),
        vm["mutation-gradient"].as<double>()
    };

    maxis::ParallelGeneticMaxisSolver solver{graph, genetic::AlgorithmStrategy(sel, rec, mut)};
    solver.constraint = vm["target"].as<double>();
    solver.size = pop_size;
    solver.migration_period = 2048;

    // Solve and print results
    auto result = graph.weighted_maxis(solver);

    auto is_valid_result = graph.is_independent_set(result);
    cout << endl << "Results" << endl << "=======" << endl;
    if (!is_valid_result) {
        cout << "ERROR Invalid Independent Set: ";
    } else {
        cout << "Vertices: ";
    }

    for (auto it = begin(result); it != end(result); ++it) {
        if (*it) {
            cout << std::distance(begin(result), it) << " ";
        }
    }
    cout << endl;

    if (is_valid_result) {
        cout << "Weight: " << graph.weighted_total(result) << endl;
    }
}
