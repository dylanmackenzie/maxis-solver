#include <unordered_set>

#include "maxis/hash_table.hpp"

namespace maxis{

size_t
BitVectorHashTable::HashFunctor::operator()(BitVector *chromosome) const {
    return std::hash<BitVector>()(*chromosome);
}

bool
BitVectorHashTable::EqFunctor::operator()(BitVector *lhs, BitVector *rhs) const {
    return std::equal(begin(*lhs), end(*lhs), begin(*rhs));
}

// TODO: determine optimal hash table size
BitVectorHashTable::BitVectorHashTable(size_t size)
    : table(6*size + 1, HashFunctor(), EqFunctor()) {}

bool
BitVectorHashTable::insert(BitVector* bv) {
    return table.insert(bv).second;
}

bool
BitVectorHashTable::erase(BitVector* bv) {
    return table.erase(bv);
}

void
BitVectorHashTable::clear() {
    table.clear();
}

} // namespace maxis
