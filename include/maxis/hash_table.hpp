#include <cstddef>
#include <unordered_set>
#include <set>

#ifndef MAXIS_HASH_TABLE_H
#define MAXIS_HASH_TABLE_H

namespace maxis {

template<typename Bv>
class BitVectorHashTable {
public:
    BitVectorHashTable(size_t);

    bool insert(Bv* bv);
    bool erase(Bv* bv);
    void clear();

private:
    struct CmpFunctor {
        bool operator()(Bv*, Bv*) const;
    };

    std::set<Bv*, CmpFunctor> table;
};

template<typename Bv>
BitVectorHashTable<Bv>::BitVectorHashTable(size_t size) : table(CmpFunctor()) {}

template<typename Bv>
bool
BitVectorHashTable<Bv>::CmpFunctor::operator()(Bv *lhs, Bv *rhs) const {
    return *lhs < *rhs;
}

template<typename Bv>
bool
BitVectorHashTable<Bv>::insert(Bv* bv) {
    return table.insert(bv).second;
}

template<typename Bv>
bool
BitVectorHashTable<Bv>::erase(Bv* bv) {
    return table.erase(bv);
}

template<typename Bv>
void
BitVectorHashTable<Bv>::clear() {
    table.clear();
}


// Specialization
// ==============

// Use an unordered map for vector<bool> so we can take advantage of the
// hash function
template<>
class BitVectorHashTable<std::vector<bool>> {
public:
    BitVectorHashTable(size_t);

    bool insert(std::vector<bool> *bv);
    bool erase(std::vector<bool> *bv);
    void clear();

private:
    struct HashFunctor {
        size_t operator()(std::vector<bool> *chromosome) const {
            return std::hash<std::vector<bool>>()(*chromosome);
        }
    };

    struct EqFunctor {
        bool operator()(std::vector<bool>* lhs, std::vector<bool>* rhs) const {
            return std::equal(begin(*lhs), end(*lhs), begin(*rhs));
        }
    };

    std::unordered_set<std::vector<bool>*, HashFunctor, EqFunctor> table;
};

// BitVectorHashTable<std::vector<bool>>::BitVectorHashTable(size_t size);
// bool BitVectorHashTable<std::vector<bool>>::insert(std::vector<bool>* bv);
// bool BitVectorHashTable<std::vector<bool>>::erase(std::vector<bool>* bv);
// void BitVectorHashTable<std::vector<bool>>::clear();

} // namespace maxis

#endif
