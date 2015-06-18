#include <algorithm>
#include <sstream>
#include <set>
#include <stack>
#include <vector>

#include "maxis/graph.hpp"
#include "maxis/solver.hpp"

using std::begin; using std::end;

namespace maxis {

Graph::Graph() {

}

// TODO: implement constructor from adjacency_matrix
Graph::Graph(const BitVector &adjacency_matrix) {

}

Graph
Graph::from_ascii_dimacs(std::istream &input) {
    std::string line, token;
    char c;
    size_t n1, n2;
    bool is_initialized = false;

    Graph g;

    // TODO: handle 0 and 1 indexing in dimacs
    while (std::getline(input, line)) {
        if(line.empty()) {
            continue;
        }

        std::istringstream ss{line};

        ss >> c;
        switch (c) {
        case 'p': // p FORMAT NODES EDGES (FORMAT=edge)
            ss >> token;
            if (token != "edge") throw ParseError("'p' must have FORMAT of 'edge'");
            ss >> n1;
            g.weights.resize(n1);
            std::fill(std::begin(g.weights), std::end(g.weights), 1);
            g.adjacency_list.resize(n1);
            g.order_ = n1;
            ss >> n2;
            g.size_ = n1;
            is_initialized = true;
            break;
        case 'e': // e W V (edge connecting w and v)
            if (!is_initialized) throw ParseError("'p' must be defined before edge definitions");;
            ss >> n1;
            ss >> n2;
            g.adjacency_list.at(n1-1).push_back(n2-1);
            g.adjacency_list.at(n2-1).push_back(n1-1);
            break;
        case 'n': // n ID VALUE (weight of given node)
            if (!is_initialized) throw ParseError("'p' must be defined before node definitions");;
            ss >> n1;
            ss >> g.weights.at(n1-1);
            break;
        case 'c':
        default:
            continue;
        }
    }

    if (!is_initialized) throw ParseError("'p' definition not found");

    g.adjacency_matrix_.resize(g.order_ * g.order_, 0);

    for (auto it = begin(g.adjacency_list); it != end(g.adjacency_list); ++it) {
        auto offset = g.order_ * std::distance(begin(g.adjacency_list), it);
        for (auto jt = begin(*it); jt != end(*it); ++jt) {
            g.adjacency_matrix_[offset + *jt] = 1;
        }
    }

    return g;
}


// Finds the disconnected subgraphs
std::vector<BitVector>
Graph::independent_subgraphs() const {
    std::vector<BitVector> subgraphs{};
    std::set<size_t> visited{};
    std::stack<size_t> search{};

    // Enumerate all vertices connected to a given start vertex, marking
    // all nodes we see as visited and putting them into a subgraph.
    // Once all these vertex have been visited. Pick a new, unvisited
    // vertex and repeat the search, using a new subgraph. Once all
    // nodes have been visited, we are done.
    while (visited.size() != order_) {
        subgraphs.emplace_back(order_, 0);
        auto &sub = subgraphs.back();

        // Pick the lowest-indexed, unvisited vertex as the search root.
        // This relies on visited being sorted.
        auto it = begin(visited);
        while (*it != static_cast<decltype(*it)>(std::distance(begin(visited), it))) {
            ++it;
        }
        search.push(*it);

        // Perform a DFS.
        while (!search.empty()) {
            auto curr = search.top();
            search.pop();

            // if vertex has already been visited, skip it.
            if (std::find(begin(visited), end(visited), curr) != end(visited)) {
                continue;
            }

            // Iterate over the neighbors of the current vertex, adding
            // them to the visited list, the current subgraph, and the
            // search if they have not yet been visited.
            for (auto &n: adjacency_list[curr]) {
                auto unique = visited.insert(n).second;
                if (unique) {
                    sub[n] = 1;
                    search.push(n);
                }
            }
        }
    }

    return subgraphs;
}

const BitVector&
Graph::adjacency_matrix() const {
    return adjacency_matrix_;
}

bool
Graph::is_independent_set(const BitVector &bv) const {

    auto adj = adjacency_matrix();

    for (decltype(order_) i = 0; i < order_; ++i) {
        auto offset = i * order_;
        for (decltype(order_) j = 0; j < order_; ++j) {
            if (adj[offset + j] == 0) {
                continue;
            }

            if (bv[i] + bv[j] == 2) {
                return false;
            }
        }
    }

    return true;
}

double
Graph::weighted_total(const BitVector &c) const {
    return std::inner_product(std::begin(c), std::end(c), begin(weights), 0.0);
}

BitVector
Graph::weighted_maxis() const {
    BruteForceMaxisSolver solve{*this};

    return solve();
}

BitVector
Graph::weighted_maxis(MaxisSolver &solve) const {
    return solve();
}

void
Graph::reorder(std::vector<size_t> &perm) {

    // Apply permutation to adjacency list and weights
    decltype(weights) tmp_weights(order_);
    decltype(adjacency_list) tmp_list(order_);
    std::transform(begin(perm), end(perm), begin(tmp_list), [this, &perm](size_t i) {
        return adjacency_list[i];
    });
    std::transform(begin(perm), end(perm), begin(tmp_weights), [this](size_t i) {
        return weights[i];
    });
    for (auto &list : tmp_list) {
        for (auto &index : list) {
            index = std::distance(begin(perm), std::find(begin(perm), end(perm), index));
        }
    }
    adjacency_list = std::move(tmp_list);
    weights = std::move(tmp_weights);

    // Recompute adjacency matrix
    std::fill(begin(adjacency_matrix_), end(adjacency_matrix_), 0);
    for (auto it = begin(adjacency_list); it != end(adjacency_list); ++it) {
        auto offset = order_ * std::distance(begin(adjacency_list), it);
        for (auto jt = begin(*it); jt != end(*it); ++jt) {
            adjacency_matrix_[offset + *jt] = 1;
        }
    }
}

} // namespace maxis
