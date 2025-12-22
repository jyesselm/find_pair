/**
 * @file interaction_filter.cpp
 * @brief Implementation of H-bond interaction filter
 */

#include <x3dna/algorithms/hydrogen_bond/interaction_filter.hpp>
#include <algorithm>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

using namespace x3dna::core;

std::vector<HBond> InteractionFilter::filter(const std::vector<HBond>& hbonds,
                                              HBondInteractionType allowed_types) {
    std::vector<HBond> result;
    result.reserve(hbonds.size());

    for (const auto& hbond : hbonds) {
        if (matches(hbond, allowed_types)) {
            result.push_back(hbond);
        }
    }

    return result;
}

bool InteractionFilter::matches(const HBond& hbond, HBondInteractionType allowed_types) {
    // ANY allows all interactions
    if (allowed_types == HBondInteractionType::ANY) {
        return true;
    }

    const HBondInteractionType interaction = context_to_interaction_type(hbond.context);
    return interaction & allowed_types;
}

HBondInteractionType InteractionFilter::context_to_interaction_type(HBondContext context) {
    switch (context) {
        case HBondContext::BASE_BASE:
            return HBondInteractionType::BASE_BASE;

        case HBondContext::BASE_BACKBONE:
            return HBondInteractionType::BASE_BACKBONE;

        case HBondContext::BASE_SUGAR:
            return HBondInteractionType::BASE_SUGAR;

        case HBondContext::BACKBONE_BACKBONE:
        case HBondContext::SUGAR_SUGAR:
            return HBondInteractionType::RNA_INTERNAL;

        case HBondContext::BASE_PROTEIN:
        case HBondContext::SUGAR_PROTEIN:
        case HBondContext::BACKBONE_PROTEIN:
            return HBondInteractionType::BASE_PROTEIN;

        case HBondContext::BASE_LIGAND:
            return HBondInteractionType::BASE_LIGAND;

        case HBondContext::PROTEIN_MAINCHAIN:
        case HBondContext::PROTEIN_SIDECHAIN:
            return HBondInteractionType::PROTEIN_PROTEIN;

        case HBondContext::PROTEIN_LIGAND:
            return HBondInteractionType::PROTEIN_LIGAND;

        case HBondContext::LIGAND_LIGAND:
            return HBondInteractionType::ANY; // No specific flag for ligand-ligand

        case HBondContext::UNKNOWN:
        default:
            return HBondInteractionType::ANY; // Unknown defaults to match ANY
    }
}

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
