#include <iostream>
#include <stdexcept>
#include <vector>

#ifndef MAXIS_GRAPH_H
#define MAXIS_GRAPH_H

namespace maxis {

// TODO: modify code to allow use of vector<bool> or vector<char>
using BitVector = std::vector<bool>;

// Forward declaration of solvers
class MaxisSolver;
class BruteForceMaxisSolver;
class GeneticMaxisSolver;

// Error thrown when parsing a text file into a graph
// TODO: use proper class declaration
using ParseError = std::runtime_error;

class Graph {
    // TODO: remove friendship dependency
    friend BruteForceMaxisSolver;
    friend GeneticMaxisSolver;

public:
    // Constructors

    Graph();
    // Construct a graph directly from an adjacency matrix
    Graph(const BitVector&);
    // Construct a graph as a subset of another graph
    Graph(const Graph&, const BitVector&);

    // Factory methods

    static Graph from_ascii_dimacs(std::istream&);

    // Methods

    // computes weighted maxis with the given algorithm
    BitVector weighted_maxis() const;
    BitVector weighted_maxis(MaxisSolver&) const;

    // Gets the weighted total of a given coloring
    double weighted_total(const BitVector&) const;

    // Returns true if a bit vector forms an independent set on the
    // graph
    bool is_independent_set(const BitVector&) const;

    // Returns the adjacency matrix of the graph
    // TODO: No need for the adjacency matrix to be the same type as a graph coloring
    BitVector adjacency_matrix() const;

    // Finds the independent subgraphs
    std::vector<BitVector> independent_subgraphs() const;

    // Reorders the graph vertices according to the given permutation
    void reorder(std::vector<size_t>&);

    // Getters

    unsigned int order() const { return _order; };
    unsigned int size() const { return _size; };

private:
    std::vector<double> weights;   // vertex weights
    std::vector<std::vector<size_t>> adjacency_list; // adjacency list for vertices

    // These are prefixed because the have a getter with the same name
    // (_lowercase_name is not illegal when not in global scope).
    unsigned int _order; // number of vertices
    unsigned int _size; // number of edges
};

} // end namespace maxis

#endif // MAXIS_GRAPH_H
