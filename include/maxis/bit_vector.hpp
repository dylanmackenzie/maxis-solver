#include "boost/dynamic_bitset.hpp"
#include "boost/iterator/iterator_facade.hpp"

#ifndef MAXIS_BIT_VECTOR_H
#define MAXIS_BIT_VECTOR_H

namespace maxis {

// Could be any random access container including vector<bool>
using BitVector = boost::dynamic_bitset<>;

namespace detail {
    // TODO: long int may not be large enough to hold a distance value
    // in all cases. Find a suitable replacement.
    template <typename ReferenceType>
    class bitset_iterator :
        public boost::iterator_facade<
            bitset_iterator<ReferenceType>,
            bool,
            boost::random_access_traversal_tag,
            ReferenceType,
            long int>
        {

    public:
        bitset_iterator() : b{nullptr}, pos{0} {}
        bitset_iterator(boost::dynamic_bitset<>* b, boost::dynamic_bitset<>::size_type i = 0)
            : b{b}, pos{i} {}

        template <typename OtherReferenceType>
        bitset_iterator(const bitset_iterator<OtherReferenceType> &other)
            : b{other.b}, pos{other.pos} {}

    private:
        friend class boost::iterator_core_access;

        boost::dynamic_bitset<> *b;
        boost::dynamic_bitset<>::size_type pos;

        ReferenceType dereference() const { return (*b)[pos]; }
        void increment() { ++pos; }
        void decrement() { --pos; }
        void advance(long int n) { pos += n; }

        template <typename OtherReferenceType>
        bool equal(const bitset_iterator<OtherReferenceType> &other) const {
            return this->b == other.b && this->pos == other.pos;
        }

        template <typename OtherReferenceType>
        long int distance_to(const bitset_iterator<OtherReferenceType> &other) const {
            return this->pos - other.pos;
        }
    };
} // namespace detail

using bitset_iterator = detail::bitset_iterator<boost::dynamic_bitset<>::reference>;
using const_bitset_iterator = detail::bitset_iterator<boost::dynamic_bitset<>::const_reference>;

// Overload begin and end for a boost::dynamic_bitset in our local
// namespace. Becasue most code uses using declarations for std::begin
// and std::end, these will be called instead of their standard library
// counterparts.
auto begin(boost::dynamic_bitset<> &b) -> maxis::bitset_iterator;
auto end(boost::dynamic_bitset<> &b) -> maxis::bitset_iterator;

} // namespace maxis

#endif // MAXIS_BIT_VECTOR_H
