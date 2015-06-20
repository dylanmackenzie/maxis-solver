#include "maxis/bit_vector.hpp"

namespace std {

maxis::bitset_iterator
begin(boost::dynamic_bitset<> &b) {
    return maxis::bitset_iterator(&b, 0);
}

maxis::bitset_iterator
end(boost::dynamic_bitset<> &b) {
    return maxis::bitset_iterator(&b, b.size());
}

} // namespace std
