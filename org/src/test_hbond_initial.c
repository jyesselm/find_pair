/**
 * @file test_hbond_initial.c
 * @brief Test ONLY the initial H-bond detection (good_hbatoms + within_limits)
 * 
 * This tool isolates the initial H-bond detection step:
 * - Which atoms are checked (seidx range)
 * - Which H-bonds pass good_hbatoms
 * - Which H-bonds pass within_limits
 * - Final initial H-bonds found
 * 
 * Does NOT do conflict resolution or validation - just initial detection
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "x3dna.h"

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <pdb_file> <residue_i> <residue_j>\n", argv[0]);
        fprintf(stderr, "Example: %s data/pdb/3G8T.pdb 946 947\n", argv[0]);
        return 1;
    }
    
    char *pdb_file = argv[1];
    long residue_i = atol(argv[2]);
    long residue_j = atol(argv[3]);
    
    // Initialize globals
    set_my_globals(argv[0]);
    
    // Count atoms
    long num = number_of_atoms(pdb_file, TRUE, "*");
    if (num <= 0) {
        fprintf(stderr, "Error: No atoms found in %s\n", pdb_file);
        return 1;
    }
    
    // Allocate arrays
    char **AtomName, **ResName, *ChainID, **Miscs;
    long *ResSeq, *idx;
    double **xyz;
    
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    
    // Read PDB
    read_pdb(pdb_file, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    
    // Get atom indices
    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);
    
    // Get residue indices
    long num_residue;
    long **seidx;
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    
    if (residue_i < 1 || residue_i > num_residue || residue_j < 1 || residue_j > num_residue) {
        fprintf(stderr, "Error: Residue indices out of range (1-%ld)\n", num_residue);
        return 1;
    }
    
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "Initial H-bond Detection Test\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "Pair: (%ld, %ld)\n", residue_i, residue_j);
    fprintf(stdout, "Residue i: %s (chain %c, seq %ld)\n", 
            ResName[seidx[residue_i][1]], ChainID[seidx[residue_i][1]], ResSeq[seidx[residue_i][1]]);
    fprintf(stdout, "Residue j: %s (chain %c, seq %ld)\n", 
            ResName[seidx[residue_j][1]], ChainID[seidx[residue_j][1]], ResSeq[seidx[residue_j][1]]);
    
    // Print seidx ranges
    fprintf(stdout, "\nAtom ranges (seidx):\n");
    fprintf(stdout, "  Residue %ld: atoms [%ld, %ld] (total: %ld atoms)\n", 
            residue_i, seidx[residue_i][1], seidx[residue_i][2], 
            seidx[residue_i][2] - seidx[residue_i][1] + 1);
    fprintf(stdout, "  Residue %ld: atoms [%ld, %ld] (total: %ld atoms)\n", 
            residue_j, seidx[residue_j][1], seidx[residue_j][2],
            seidx[residue_j][2] - seidx[residue_j][1] + 1);
    
    // Print atoms in range
    fprintf(stdout, "\nAtoms in residue %ld:\n", residue_i);
    for (long m = seidx[residue_i][1]; m <= seidx[residue_i][2]; m++) {
        fprintf(stdout, "  [%ld] %s (idx=%ld)\n", m, AtomName[m], idx[m]);
    }
    
    fprintf(stdout, "\nAtoms in residue %ld:\n", residue_j);
    for (long n = seidx[residue_j][1]; n <= seidx[residue_j][2]; n++) {
        fprintf(stdout, "  [%ld] %s (idx=%ld)\n", n, AtomName[n], idx[n]);
    }
    
    // Now do initial H-bond detection (matches get_hbond_ij initial loop)
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "Initial H-bond Detection\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "Checking: good_hbatoms() && within_limits()\n");
    fprintf(stdout, "Distance range: [%.3f, %.3f]\n", 
            Gvars.misc_pars.hb_lower, Gvars.misc_pars.hb_dist1);
    
    long num_hbonds = 0;
    char hb_atom1[BUF512][5], hb_atom2[BUF512][5];
    double hb_dist[BUF512];
    
    for (long m = seidx[residue_i][1]; m <= seidx[residue_i][2]; m++) {
        for (long n = seidx[residue_j][1]; n <= seidx[residue_j][2]; n++) {
            // Check distance first
            double dist = p1p2_dist(xyz[n], xyz[m]);
            int in_range = within_limits(xyz[n], xyz[m], Gvars.misc_pars.hb_lower, Gvars.misc_pars.hb_dist1);
            
            // Check good_hbatoms
            int good_atoms = good_hbatoms(&Gvars.misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]);
            
            if (good_atoms && in_range) {
                num_hbonds++;
                strcpy(hb_atom1[num_hbonds], AtomName[m]);
                strcpy(hb_atom2[num_hbonds], AtomName[n]);
                hb_dist[num_hbonds] = dist;
                
                fprintf(stdout, "\nH-bond #%ld:\n", num_hbonds);
                fprintf(stdout, "  %s [%ld] -> %s [%ld]\n", AtomName[m], m, AtomName[n], n);
                fprintf(stdout, "  Distance: %.6f\n", dist);
                fprintf(stdout, "  idx[m]=%ld, idx[n]=%ld\n", idx[m], idx[n]);
            } else {
                // Debug: show why it was rejected
                if (dist >= Gvars.misc_pars.hb_lower && dist <= Gvars.misc_pars.hb_dist1 && !good_atoms) {
                    fprintf(stdout, "  REJECTED: %s -> %s (dist=%.3f) - good_hbatoms failed\n", 
                            AtomName[m], AtomName[n], dist);
                } else if (good_atoms && !in_range) {
                    fprintf(stdout, "  REJECTED: %s -> %s (dist=%.3f) - out of range\n", 
                            AtomName[m], AtomName[n], dist);
                }
            }
        }
    }
    
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "Summary\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "Total initial H-bonds found: %ld\n", num_hbonds);
    fprintf(stdout, "\nInitial H-bonds:\n");
    for (long k = 1; k <= num_hbonds; k++) {
        fprintf(stdout, "  %ld. %s -> %s, dist=%.6f\n", k, hb_atom1[k], hb_atom2[k], hb_dist[k]);
    }
    
    // Output JSON for easy parsing
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "JSON Output\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "{\n");
    fprintf(stdout, "  \"residue_i\": %ld,\n", residue_i);
    fprintf(stdout, "  \"residue_j\": %ld,\n", residue_j);
    fprintf(stdout, "  \"num_initial_hbonds\": %ld,\n", num_hbonds);
    fprintf(stdout, "  \"seidx_i\": [%ld, %ld],\n", seidx[residue_i][1], seidx[residue_i][2]);
    fprintf(stdout, "  \"seidx_j\": [%ld, %ld],\n", seidx[residue_j][1], seidx[residue_j][2]);
    fprintf(stdout, "  \"hbonds\": [\n");
    for (long k = 1; k <= num_hbonds; k++) {
        fprintf(stdout, "    {\n");
        fprintf(stdout, "      \"hbond_idx\": %ld,\n", k);
        fprintf(stdout, "      \"donor_atom\": \"%s\",\n", hb_atom1[k]);
        fprintf(stdout, "      \"acceptor_atom\": \"%s\",\n", hb_atom2[k]);
        fprintf(stdout, "      \"distance\": %.6f\n", hb_dist[k]);
        fprintf(stdout, "    }%s\n", k < num_hbonds ? "," : "");
    }
    fprintf(stdout, "  ]\n");
    fprintf(stdout, "}\n");
    
    return 0;
}

