/**
 * @file test_hbond_validation.c
 * @brief Test ONLY validation (validate_hbonds)
 * 
 * This tool isolates the validation step:
 * - Takes H-bonds after conflict resolution as input
 * - Applies validate_hbonds
 * - Shows which H-bonds get which types (' ', '-', '*')
 * 
 * Does NOT do conflict resolution - assumes it's already done
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
    long *RY;
    char *bseq;
    
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
    
    // Get base types
    RY = lvector(1, num_residue);
    bseq = cvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);
    
    char basei = bseq[residue_i];
    char basej = bseq[residue_j];
    
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "H-bond Validation Test\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "Pair: (%ld, %ld)\n", residue_i, residue_j);
    fprintf(stdout, "Base i: %c, Base j: %c\n", basei, basej);
    
    // Step 1: Get initial H-bonds
    long num_hbonds = 0;
    char **hb_atom1, **hb_atom2;
    double *hb_dist;
    long *lkg_type;
    hb_atom1 = cmatrix(1, BUF512, 0, 4);
    hb_atom2 = cmatrix(1, BUF512, 0, 4);
    hb_dist = dvector(1, BUF512);
    lkg_type = lvector(1, BUF512);
    
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
    
    if (num_hbonds == 0) {
        fprintf(stdout, "\nNo H-bonds to validate.\n");
        return 0;
    }
    
    // Step 2: Conflict resolution
    fprintf(stdout, "\nStep 2: Conflict resolution (hb_atompair)...\n");
    hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, &Gvars.misc_pars);
    
    long num_after_conflict = 0;
    for (long k = 1; k <= num_hbonds; k++) {
        if (hb_dist[k] > 0) {
            num_after_conflict++;
        }
    }
    fprintf(stdout, "After conflict resolution: %ld H-bonds (positive distance)\n", num_after_conflict);
    
    // Step 3: Validation
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "Step 3: Validation (validate_hbonds)\n");
    fprintf(stdout, "========================================\n");
    
    char *hb_type;
    hb_type = cvector(1, num_hbonds);
    
    long num_validated = validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej,
                                         hb_atom1, hb_atom2);
    
    fprintf(stdout, "\nAfter validation:\n");
    long count_type_space = 0, count_type_dash = 0, count_type_star = 0;
    for (long k = 1; k <= num_hbonds; k++) {
        if (hb_dist[k] > 0) {  // Only show non-conflicts
            char type_char = hb_type[k];
            if (type_char == ' ') count_type_space++;
            else if (type_char == '-') count_type_dash++;
            else if (type_char == '*') count_type_star++;
            
            fprintf(stdout, "  %ld. %s -> %s, dist=%.6f, type='%c', lkg=%ld\n", 
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k], type_char, lkg_type[k]);
        }
    }
    
    fprintf(stdout, "\nSummary:\n");
    fprintf(stdout, "  Total initial H-bonds: %ld\n", num_hbonds);
    fprintf(stdout, "  After conflict resolution: %ld\n", num_after_conflict);
    fprintf(stdout, "  After validation: %ld\n", num_validated);
    fprintf(stdout, "  Type ' ': %ld\n", count_type_space);
    fprintf(stdout, "  Type '-': %ld\n", count_type_dash);
    fprintf(stdout, "  Type '*': %ld\n", count_type_star);
    
    // Output JSON
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "JSON Output\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "{\n");
    fprintf(stdout, "  \"residue_i\": %ld,\n", residue_i);
    fprintf(stdout, "  \"residue_j\": %ld,\n", residue_j);
    fprintf(stdout, "  \"base_i\": \"%c\",\n", basei);
    fprintf(stdout, "  \"base_j\": \"%c\",\n", basej);
    fprintf(stdout, "  \"num_initial_hbonds\": %ld,\n", num_hbonds);
    fprintf(stdout, "  \"num_after_conflict\": %ld,\n", num_after_conflict);
    fprintf(stdout, "  \"num_after_validation\": %ld,\n", num_validated);
    fprintf(stdout, "  \"hbonds\": [\n");
    long json_count = 0;
    for (long k = 1; k <= num_hbonds; k++) {
        if (hb_dist[k] > 0) {  // Only include non-conflicts
            json_count++;
            fprintf(stdout, "    {\n");
            fprintf(stdout, "      \"hbond_idx\": %ld,\n", k);
            fprintf(stdout, "      \"donor_atom\": \"%s\",\n", hb_atom1[k]);
            fprintf(stdout, "      \"acceptor_atom\": \"%s\",\n", hb_atom2[k]);
            fprintf(stdout, "      \"distance\": %.6f,\n", hb_dist[k]);
            fprintf(stdout, "      \"type\": \"%c\",\n", hb_type[k]);
            fprintf(stdout, "      \"linkage_type\": %ld\n", lkg_type[k]);
            fprintf(stdout, "    }%s\n", json_count < num_after_conflict ? "," : "");
        }
    }
    fprintf(stdout, "  ]\n");
    fprintf(stdout, "}\n");
    
    free_cmatrix(hb_atom1, 1, BUF512, 0, 4);
    free_cmatrix(hb_atom2, 1, BUF512, 0, 4);
    free_dvector(hb_dist, 1, BUF512);
    free_lvector(lkg_type, 1, BUF512);
    free_cvector(hb_type, 1, num_hbonds);
    free_lvector(RY, 1, num_residue);
    free_cvector(bseq, 1, num_residue);
    
    return 0;
}

