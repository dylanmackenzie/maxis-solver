#include <unordered_set>

#ifndef MAXIS_HASH_TABLE_H
#define MAXIS_HASH_TABLE_H

namespace maxis {

class BitVectorHashTable {
public:
    BitVectorHashTable(size_t size);

    bool insert(BitVector* bv);
    bool erase(BitVector* bv);
    void clear();

private:
    struct HashFunctor {
        size_t operator()(BitVector*) const;
    };

    struct EqFunctor {
        bool operator()(BitVector*, BitVector*) const;
    };

    // Using functors here avoids the overhead of std::function
    std::unordered_set<BitVector*,
        BitVectorHashTable::HashFunctor,
        BitVectorHashTable::EqFunctor> table;
};

} // namespace maxis

#endif
