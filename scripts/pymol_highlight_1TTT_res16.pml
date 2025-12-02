# PyMOL Script to Highlight Chain D, ResSeq 16 (H2U) in 1TTT
# This is the residue that modern recognizes but legacy does not
# ============================================================

# Load the PDB file
load data/pdb/1TTT.pdb, 1TTT

# Show the structure
show cartoon, 1TTT
show sticks, 1TTT and resn H2U

# Highlight the problematic residue (Chain D, ResSeq 16)
# This is the ONLY H2U that legacy does NOT recognize
select problematic_res, 1TTT and chain D and resi 16
show sticks, problematic_res
show spheres, problematic_res
color red, problematic_res

# Also show the other H2U residues for comparison
# These are all recognized by both legacy and modern
select other_h2u, 1TTT and resn H2U and not (chain D and resi 16)
show sticks, other_h2u
color yellow, other_h2u

# Show residue 59 (the partner that modern pairs with residue 16)
select res59, 1TTT and chain D and resi 59
show sticks, res59
show spheres, res59
color blue, res59

# Center view on the problematic residue
center problematic_res
zoom problematic_res, 10

# Label the residue
label problematic_res and name CA, "%s-%s%s" % (resn, chain, resi)

# Print summary
print("=" * 60)
print("1TTT Residue Highlighting")
print("=" * 60)
print("RED: Chain D, ResSeq 16 (H2U) - NOT recognized by legacy")
print("YELLOW: Other H2U residues - recognized by both")
print("BLUE: Chain D, ResSeq 59 - paired with residue 16 by modern")
print("=" * 60)

