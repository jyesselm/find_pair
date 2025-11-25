/**
 * @file test_hbond_conflict.c
 * @brief Test ONLY conflict resolution (hb_atompair)
 * 
 * This tool isolates the conflict resolution step:
 * - Takes initial H-bonds as input
 * - Applies hb_atompair algorithm
 * - Shows which H-bonds are marked as conflicts
 * - Shows linkage types assigned
 * 
 * Does NOT do validation - just conflict resolution
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
    fprintf(stdout, "H-bond Conflict Resolution Test\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "Pair: (%ld, %ld)\n", residue_i, residue_j);
    
    // Step 1: Get initial H-bonds (same as test_hbond_initial)
    long num_hbonds = 0;
    char **hb_atom1, **hb_atom2;
    double *hb_dist;
    hb_atom1 = cmatrix(1, BUF512, 0, 4);
    hb_atom2 = cmatrix(1, BUF512, 0, 4);
    hb_dist = dvector(1, BUF512);
    
    fprintf(stdout, "\nStep 1: Finding initial H-bonds...\n");
    for (long m = seidx[residue_i][1]; m <= seidx[residue_i][2]; m++) {
        for (long n = seidx[residue_j][1]; n <= seidx[residue_j][2]; n++) {
            if (good_hbatoms(&Gvars.misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]) &&
                within_limits(xyz[n], xyz[m], Gvars.misc_pars.hb_lower, Gvars.misc_pars.hb_dist1)) {
                num_hbonds++;
                strcpy(hb_atom1[num_hbonds], AtomName[m]);
                strcpy(hb_atom2[num_hbonds], AtomName[n]);
                hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m]);
            }
        }
    }
    
    fprintf(stdout, "Found %ld initial H-bonds\n", num_hbonds);
    for (long k = 1; k <= num_hbonds; k++) {
        fprintf(stdout, "  %ld. %s -> %s, dist=%.6f\n", k, hb_atom1[k], hb_atom2[k], hb_dist[k]);
    }
    
    if (num_hbonds == 0) {
        fprintf(stdout, "\nNo H-bonds to resolve conflicts for.\n");
        return 0;
    }
    
    // Step 2: Apply conflict resolution (hb_atompair)
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "Step 2: Conflict Resolution (hb_atompair)\n");
    fprintf(stdout, "========================================\n");
    
    long *lkg_type;
    lkg_type = lvector(1, num_hbonds);
    
    // Save original distances
    double *hb_dist_original;
    hb_dist_original = dvector(1, num_hbonds);
    for (long k = 1; k <= num_hbonds; k++) {
        hb_dist_original[k] = hb_dist[k];
    }
    
    fprintf(stdout, "Calling hb_atompair...\n");
    hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, &Gvars.misc_pars);
    
    fprintf(stdout, "\nAfter conflict resolution:\n");
    long num_conflicts = 0;
    long num_kept = 0;
    for (long k = 1; k <= num_hbonds; k++) {
        int is_conflict = (hb_dist[k] < 0);
        if (is_conflict) {
            num_conflicts++;
            fprintf(stdout, "  %ld. %s -> %s, dist=%.6f (NEGATED - conflict), lkg=%ld\n", 
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k], lkg_type[k]);
        } else {
            num_kept++;
            fprintf(stdout, "  %ld. %s -> %s, dist=%.6f (kept), lkg=%ld\n", 
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k], lkg_type[k]);
        }
    }
    
    fprintf(stdout, "\nSummary:\n");
    fprintf(stdout, "  Total H-bonds: %ld\n", num_hbonds);
    fprintf(stdout, "  Kept (positive distance): %ld\n", num_kept);
    fprintf(stdout, "  Conflicts (negative distance): %ld\n", num_conflicts);
    
    // Output JSON
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "JSON Output\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "{\n");
    fprintf(stdout, "  \"residue_i\": %ld,\n", residue_i);
    fprintf(stdout, "  \"residue_j\": %ld,\n", residue_j);
    fprintf(stdout, "  \"num_initial_hbonds\": %ld,\n", num_hbonds);
    fprintf(stdout, "  \"num_kept\": %ld,\n", num_kept);
    fprintf(stdout, "  \"num_conflicts\": %ld,\n", num_conflicts);
    fprintf(stdout, "  \"hbonds\": [\n");
    for (long k = 1; k <= num_hbonds; k++) {
        fprintf(stdout, "    {\n");
        fprintf(stdout, "      \"hbond_idx\": %ld,\n", k);
        fprintf(stdout, "      \"donor_atom\": \"%s\",\n", hb_atom1[k]);
        fprintf(stdout, "      \"acceptor_atom\": \"%s\",\n", hb_atom2[k]);
        fprintf(stdout, "      \"distance_original\": %.6f,\n", hb_dist_original[k]);
        fprintf(stdout, "      \"distance_after_conflict\": %.6f,\n", hb_dist[k]);
        fprintf(stdout, "      \"is_conflict\": %s,\n", hb_dist[k] < 0 ? "true" : "false");
        fprintf(stdout, "      \"linkage_type\": %ld\n", lkg_type[k]);
        fprintf(stdout, "    }%s\n", k < num_hbonds ? "," : "");
    }
    fprintf(stdout, "  ]\n");
    fprintf(stdout, "}\n");
    
    free_cmatrix(hb_atom1, 1, BUF512, 0, 4);
    free_cmatrix(hb_atom2, 1, BUF512, 0, 4);
    free_dvector(hb_dist, 1, BUF512);
    free_dvector(hb_dist_original, 1, num_hbonds);
    free_lvector(lkg_type, 1, num_hbonds);
    
    return 0;
}

