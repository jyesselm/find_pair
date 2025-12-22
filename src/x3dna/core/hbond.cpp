/**
 * @file hbond.cpp
 * @brief Hydrogen bond implementation
 */

#include <x3dna/core/hbond.hpp>

namespace x3dna {
namespace core {

int HBond::legacy_linkage_type() const {
    // Legacy linkage_type mapping from idx2 array sum
    // idx2[k][0] tracks donor sharing, idx2[k][1] tracks acceptor sharing
    // linkage_type = idx2[k][0] + idx2[k][1]
    switch (conflict_state) {
        case ConflictState::NO_CONFLICT:
            return 0;
        case ConflictState::SHARES_DONOR_WITH_WINNER:
            return 1; // idx2[m][0]=1, idx2[m][1]=0
        case ConflictState::SHARES_ACCEPTOR_WITH_WINNER:
            return 1; // idx2[m][0]=0, idx2[m][1]=1
        case ConflictState::SHARES_BOTH_WITH_WINNER:
            return 2; // idx2[m][0]=1, idx2[m][1]=1
        case ConflictState::IS_CONFLICT_WINNER:
            return 18; // idx2[k][0]=9, idx2[k][1]=9 -> 9+9=18
        default:
            return 0;
    }
}

} // namespace core
} // namespace x3dna
