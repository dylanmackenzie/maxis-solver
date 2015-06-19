#include <vector>

#include "maxis/hash_table.hpp"

namespace maxis {

// TODO: find optimal hash table size
BitVectorHashTable<std::vector<bool>>::BitVectorHashTable(size_t size)
    : table(6*size + 1, HashFunctor(), EqFunctor()) {}

bool
BitVectorHashTable<std::vector<bool>>::insert(std::vector<bool>* bv) {
    return table.insert(bv).second;
}

bool
BitVectorHashTable<std::vector<bool>>::erase(std::vector<bool>* bv) {
    return table.erase(bv);
}

void
BitVectorHashTable<std::vector<bool>>::clear() {
    table.clear();
}

} // namespace maxis
