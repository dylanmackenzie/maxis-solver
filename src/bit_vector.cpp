#include "maxis/bit_vector.hpp"

namespace maxis {

bitset_iterator
begin(boost::dynamic_bitset<> &b) {
    return bitset_iterator(&b);
}

const_bitset_iterator
begin(const boost::dynamic_bitset<> &b) {
    return const_bitset_iterator(&b);
}

bitset_iterator
end(boost::dynamic_bitset<> &b) {
    return bitset_iterator(&b, b.size());
}

const_bitset_iterator
end(const boost::dynamic_bitset<> &b) {
    return const_bitset_iterator(&b, b.size());
}

const_bitset_iterator
cbegin(const boost::dynamic_bitset<> &b) {
    return const_bitset_iterator(&b);
}

const_bitset_iterator
cend(const boost::dynamic_bitset<> &b) {
    return const_bitset_iterator(&b, b.size());
}

} // namespace std
