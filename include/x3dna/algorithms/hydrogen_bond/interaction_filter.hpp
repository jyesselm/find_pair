/**
 * @file interaction_filter.hpp
 * @brief Filter H-bonds by interaction type using bitwise flags
 */

#pragma once

#include <vector>
#include <x3dna/algorithms/hydrogen_bond/hbond.hpp>
#include <x3dna/algorithms/hydrogen_bond/hbond_types.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

/**
 * @brief Filter H-bonds by interaction type
 *
 * Uses HBondInteractionType bitwise flags to select which types of
 * interactions to include. Allows combining multiple types with |.
 *
 * Example:
 * @code
 * // Filter for only base-base interactions
 * auto filtered = InteractionFilter::filter(hbonds, HBondInteractionType::BASE_BASE);
 *
 * // Filter for base-base OR protein-protein
 * auto filtered = InteractionFilter::filter(hbonds,
 *     HBondInteractionType::BASE_BASE | HBondInteractionType::PROTEIN_PROTEIN);
 *
 * // Include all RNA internal interactions
 * auto filtered = InteractionFilter::filter(hbonds, HBondInteractionType::RNA_INTERNAL);
 * @endcode
 */
class InteractionFilter {
public:
    /**
     * @brief Filter H-bonds to only include specified interaction types
     * @param hbonds Input H-bonds
     * @param allowed_types Bitwise combination of allowed HBondInteractionType flags
     * @return Vector of H-bonds matching the allowed types
     */
    [[nodiscard]] static std::vector<core::HBond> filter(const std::vector<core::HBond>& hbonds,
                                                          core::HBondInteractionType allowed_types);

    /**
     * @brief Check if an H-bond matches the allowed interaction types
     * @param hbond H-bond to check
     * @param allowed_types Bitwise combination of allowed HBondInteractionType flags
     * @return true if the H-bond's context matches allowed types
     */
    [[nodiscard]] static bool matches(const core::HBond& hbond, core::HBondInteractionType allowed_types);

    /**
     * @brief Convert HBondContext to HBondInteractionType
     * @param context H-bond context
     * @return Corresponding interaction type flag
     */
    [[nodiscard]] static core::HBondInteractionType context_to_interaction_type(core::HBondContext context);
};

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
