/**
 * @file pair_geometry_helper.cpp
 * @brief Implementation of base pair geometry calculations
 */

#include <x3dna/algorithms/helix/pair_geometry_helper.hpp>

namespace x3dna::algorithms::helix {

geometry::Vector3D PairGeometryHelper::get_pair_origin(const core::BasePair& pair) {
    // Legacy uses the average of both base origins (morg) for pair distance calculations
    // See refs_right_left() in cmn_fncs.c: morg[j] = (org[base1][j] + org[base2][j]) / 2
    // Note: organize() validates all pairs have both frames
    const auto& o1 = pair.frame1()->origin();
    const auto& o2 = pair.frame2()->origin();
    return geometry::Vector3D((o1.x() + o2.x()) / 2.0, (o1.y() + o2.y()) / 2.0, (o1.z() + o2.z()) / 2.0);
}

geometry::Vector3D PairGeometryHelper::get_pair_z_axis(const core::BasePair& pair) {
    // Legacy uses average z-axis from both frames (or difference if they point opposite)
    // This matches bp_context: (d <= 0.0) ? ddxyz(z2, z1, zave) : sumxyz(z2, z1, zave)
    // Note: organize() validates all pairs have both frames
    auto z1 = pair.frame1()->z_axis();
    auto z2 = pair.frame2()->z_axis();
    double d = z1.dot(z2);

    geometry::Vector3D zave = (d <= 0.0) ? (z2 - z1) : (z2 + z1);
    zave.normalize();
    return zave;
}

geometry::Vector3D PairGeometryHelper::get_frame_z(const core::BasePair& pair, bool swapped) {
    // Note: organize() validates all pairs have both frames
    return swapped ? pair.frame2()->z_axis() : pair.frame1()->z_axis();
}

StrandResidues PairGeometryHelper::get_strand_residues(const core::BasePair& pair, bool swapped) {
    // BasePair stores 0-based indices normalized to (smaller, larger).
    // finding_order_swapped() indicates if original finding order was (larger, smaller).
    // Legacy code uses the ORIGINAL finding order for strand assignments.
    //
    // To match legacy:
    // - Start with original finding order (apply finding_order_swapped to restore it)
    // - Then apply the five2three swap flag
    // This is equivalent to XOR of finding_order_swapped and swapped.
    bool use_reversed_order = (pair.finding_order_swapped() != swapped);

    if (use_reversed_order) {
        return {pair.residue_idx2() + 1, pair.residue_idx1() + 1};
    }
    return {pair.residue_idx1() + 1, pair.residue_idx2() + 1};
}

} // namespace x3dna::algorithms::helix
