#include "boost/dynamic_bitset.hpp"
#include "boost/iterator/iterator_facade.hpp"

#ifndef MAXIS_BIT_VECTOR_H
#define MAXIS_BIT_VECTOR_H

namespace maxis {

// Could be any random access container including vector<bool>
using BitVector = boost::dynamic_bitset<>;

class bitset_iterator :
    public boost::iterator_facade<
        bitset_iterator,
        bool,
        boost::random_access_traversal_tag,
        boost::dynamic_bitset<>::reference,
        long int >
    {

public:
    bitset_iterator() : b{nullptr}, pos{0} {}
    bitset_iterator(boost::dynamic_bitset<>* b, boost::dynamic_bitset<>::size_type i = 0) : b{b}, pos{i} {}

private:
    friend class boost::iterator_core_access;

    boost::dynamic_bitset<> *b;
    boost::dynamic_bitset<>::size_type pos;

    boost::dynamic_bitset<>::reference dereference() const { return (*b)[pos]; }
    bool equal(bitset_iterator const& other) const {
        return this->b == other.b && this->pos == other.pos;
    }
    void increment() { ++pos; }
    void decrement() { --pos; }
    void advance(long int n) { pos += n; }
    long int distance_to(boost::dynamic_bitset<>::size_type z) const { return pos - z; }
};

// Overload begin and end for a boost::dynamic_bitset in our local
// namespace. Becasue most code uses using declarations for std::begin
// and std::end, these will be called instead of their standard library
// counterparts.

auto begin(boost::dynamic_bitset<> &b) -> maxis::bitset_iterator;
auto end(boost::dynamic_bitset<> &b) -> maxis::bitset_iterator;

} // namespace maxis

#endif // MAXIS_BIT_VECTOR_H
