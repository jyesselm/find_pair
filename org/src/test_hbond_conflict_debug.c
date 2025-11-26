/**
 * @file test_hbond_conflict_debug.c
 * @brief Enhanced debug tool to trace conflict resolution step-by-step
 * 
 * This tool provides detailed step-by-step tracing of hb_atompair algorithm:
 * - Phase 1: Initial conflict detection (iterative algorithm)
 * - Phase 2: idx2 population
 * - Phase 3: Linkage type calculation and additional conflict marking
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "x3dna.h"

void hb_atompair_debug(long num_hbonds, char **hb_atom1, char **hb_atom2, double *hb_dist,
                       long *lkg_type, miscPars *misc_pars) {
    double dtmp[3];
    long k, m = 0, n, num_iter = 1;
    long ddidx[3], *matched_idx, **idx2;
    
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "PHASE 1: Initial Conflict Detection\n");
    fprintf(stdout, "========================================\n");
    
    if (!num_hbonds)
        return;
    matched_idx = lvector(1, num_hbonds);
    
    long phase1_iter = 0;
    while (1) {
        phase1_iter++;
        if (matched_idx[num_iter]) {
            num_iter++;
            continue;
        }
        
        fprintf(stdout, "\n--- Phase 1 Iteration %ld (num_iter=%ld) ---\n", phase1_iter, num_iter);
        fprintf(stdout, "Current H-bond: %s -> %s (dist=%.6f)\n", 
                hb_atom1[num_iter], hb_atom2[num_iter], hb_dist[num_iter]);
        
        for (k = 1; k <= 2; k++)
            update_hb_idx(k, dtmp, ddidx, hb_dist, num_iter);
        
        fprintf(stdout, "Initial shortest distances:\n");
        fprintf(stdout, "  dtmp[1] = %.6f (donor atom: %s)\n", dtmp[1], hb_atom1[num_iter]);
        fprintf(stdout, "  dtmp[2] = %.6f (acceptor atom: %s)\n", dtmp[2], hb_atom2[num_iter]);
        
        for (n = 1; n <= num_hbonds; n++) {
            if (n == num_iter || matched_idx[n])
                continue;
            if (!strcmp(hb_atom1[n], hb_atom1[num_iter]) && hb_dist[n] < dtmp[1]) {
                fprintf(stdout, "  Found shorter for donor: H-bond %ld (%s -> %s, dist=%.6f)\n",
                        n, hb_atom1[n], hb_atom2[n], hb_dist[n]);
                update_hb_idx(1, dtmp, ddidx, hb_dist, n);
            }
            if (!strcmp(hb_atom2[n], hb_atom2[num_iter]) && hb_dist[n] < dtmp[2]) {
                fprintf(stdout, "  Found shorter for acceptor: H-bond %ld (%s -> %s, dist=%.6f)\n",
                        n, hb_atom1[n], hb_atom2[n], hb_dist[n]);
                update_hb_idx(2, dtmp, ddidx, hb_dist, n);
            }
        }
        
        fprintf(stdout, "Final shortest:\n");
        fprintf(stdout, "  ddidx[1] = %ld, dtmp[1] = %.6f\n", ddidx[1], dtmp[1]);
        fprintf(stdout, "  ddidx[2] = %ld, dtmp[2] = %.6f\n", ddidx[2], dtmp[2]);
        
        if (ddidx[1] == ddidx[2]) {
            k = ddidx[1];
            fprintf(stdout, "  CONFLICT DETECTED! Both point to H-bond %ld\n", k);
            fprintf(stdout, "  Marking H-bond %ld (%s -> %s) as conflict (negating distance)\n",
                    k, hb_atom1[k], hb_atom2[k]);
            fprintf(stdout, "  Original distance: %.6f, new distance: %.6f\n",
                    hb_dist[k], -hb_dist[k]);
            hb_dist[k] = -hb_dist[k];
            
            num_iter = 1;
            for (n = 1; n <= num_hbonds; n++) {
                if (matched_idx[n])
                    continue;
                if (!strcmp(hb_atom1[n], hb_atom1[k]) || !strcmp(hb_atom2[n], hb_atom2[k])) {
                    matched_idx[n] = 1;
                    m++;
                    fprintf(stdout, "  Marking H-bond %ld as matched\n", n);
                }
            }
            if (m >= num_hbonds)
                break;
        } else {
            fprintf(stdout, "  No conflict (donor and acceptor point to different H-bonds)\n");
            num_iter++;
        }
    }
    
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "PHASE 2: idx2 Population\n");
    fprintf(stdout, "========================================\n");
    
    idx2 = lmatrix(1, num_hbonds, 1, 2);
    for (k = 1; k <= num_hbonds; k++) {
        if (hb_dist[k] > 0.0)
            continue;
        fprintf(stdout, "\nProcessing conflicted H-bond %ld: %s -> %s (dist=%.6f)\n",
                k, hb_atom1[k], hb_atom2[k], hb_dist[k]);
        idx2[k][1] = 9;
        idx2[k][2] = 9;
        fprintf(stdout, "  Setting idx2[%ld][1] = 9, idx2[%ld][2] = 9\n", k, k);
        
        for (m = 1; m <= num_hbonds; m++) {
            if (m == k || hb_dist[m] < 0.0)
                continue;
            if (!strcmp(hb_atom1[m], hb_atom1[k])) {
                idx2[m][1] = 1;
                fprintf(stdout, "  H-bond %ld shares atom1 (%s) -> idx2[%ld][1] = 1\n", m, hb_atom1[m], m);
            }
            if (!strcmp(hb_atom2[m], hb_atom2[k])) {
                idx2[m][2] = 1;
                fprintf(stdout, "  H-bond %ld shares atom2 (%s) -> idx2[%ld][2] = 1\n", m, hb_atom2[m], m);
            }
        }
    }
    
    fprintf(stdout, "\nidx2 values after Phase 2:\n");
    for (k = 1; k <= num_hbonds; k++) {
        fprintf(stdout, "  H-bond %ld: idx2[%ld][1] = %ld, idx2[%ld][2] = %ld\n",
                k, k, idx2[k][1], k, idx2[k][2]);
    }
    
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "PHASE 3: Linkage Type & Additional Conflicts\n");
    fprintf(stdout, "========================================\n");
    
    for (k = 1; k <= num_hbonds; k++) {
        m = idx2[k][1] + idx2[k][2];
        lkg_type[k] = m;
        fprintf(stdout, "\nH-bond %ld: %s -> %s\n", k, hb_atom1[k], hb_atom2[k]);
        fprintf(stdout, "  idx2[%ld][1] = %ld, idx2[%ld][2] = %ld\n", k, idx2[k][1], k, idx2[k][2]);
        fprintf(stdout, "  Linkage type = %ld + %ld = %ld\n", idx2[k][1], idx2[k][2], m);
        fprintf(stdout, "  Current distance: %.6f %s\n", hb_dist[k], hb_dist[k] < 0.0 ? "(CONFLICT)" : "(positive)");
        
        if (m != 18 && dval_in_range(hb_dist[k], misc_pars->hb_lower, misc_pars->hb_dist2)) {
            if (hb_dist[k] > 0.0) {
                fprintf(stdout, "  Linkage type != 18 and distance in range [%.2f, %.2f]\n",
                        misc_pars->hb_lower, misc_pars->hb_dist2);
                fprintf(stdout, "  -> Marking as additional conflict (negating distance)\n");
                fprintf(stdout, "  Original distance: %.6f, new distance: %.6f\n",
                        hb_dist[k], -hb_dist[k]);
                hb_dist[k] = -hb_dist[k];
            } else {
                fprintf(stdout, "  Already a conflict (distance < 0)\n");
            }
        } else {
            if (m == 18) {
                fprintf(stdout, "  Linkage type = 18 (no conflicts) -> keeping positive\n");
            } else {
                fprintf(stdout, "  Distance out of range [%.2f, %.2f] -> keeping positive\n",
                        misc_pars->hb_lower, misc_pars->hb_dist2);
            }
        }
    }
    
    free_lmatrix(idx2, 1, num_hbonds, 1, 2);
    free_lvector(matched_idx, 1, num_hbonds);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <pdb_file> <residue_i> <residue_j> [output_file]\n", argv[0]);
        fprintf(stderr, "Example: %s data/pdb/1VBY.pdb 45 62 debug_output.txt\n", argv[0]);
        return 1;
    }
    
    char *pdb_file = argv[1];
    long residue_i = atol(argv[2]);
    long residue_j = atol(argv[3]);
    
    FILE *output_fp = stdout;
    if (argc >= 5) {
        output_fp = fopen(argv[4], "w");
        if (!output_fp) {
            fprintf(stderr, "Error: Cannot open output file: %s\n", argv[4]);
            return 1;
        }
    }
    
    // Redirect stdout to output file
    FILE *original_stdout = stdout;
    stdout = output_fp;
    
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
    fprintf(stdout, "H-bond Conflict Resolution Debug\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "PDB: %s\n", pdb_file);
    fprintf(stdout, "Pair: (%ld, %ld)\n", residue_i, residue_j);
    
    // Get base types
    char *bseq;
    long *RY;
    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);
    char basei = bseq[residue_i];
    char basej = bseq[residue_j];
    fprintf(stdout, "Base types: %c - %c\n", basei, basej);
    
    // Step 1: Get initial H-bonds
    long num_hbonds = 0;
    char **hb_atom1, **hb_atom2;
    double *hb_dist;
    hb_atom1 = cmatrix(1, BUF512, 0, 4);
    hb_atom2 = cmatrix(1, BUF512, 0, 4);
    hb_dist = dvector(1, BUF512);
    
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "INITIAL H-BONDS (before conflict resolution)\n");
    fprintf(stdout, "========================================\n");
    fprintf(stdout, "Checking atoms in range:\n");
    fprintf(stdout, "  Residue %ld: atoms %ld-%ld\n", residue_i, seidx[residue_i][1], seidx[residue_i][2]);
    fprintf(stdout, "  Residue %ld: atoms %ld-%ld\n", residue_j, seidx[residue_j][1], seidx[residue_j][2]);
    
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
    
    fprintf(stdout, "\nFound %ld initial H-bonds:\n", num_hbonds);
    for (long k = 1; k <= num_hbonds; k++) {
        fprintf(stdout, "  %ld. %s -> %s, dist=%.6f\n", k, hb_atom1[k], hb_atom2[k], hb_dist[k]);
    }
    
    if (num_hbonds == 0) {
        fprintf(stdout, "\nNo H-bonds to resolve conflicts for.\n");
        goto cleanup;
    }
    
    // Step 2: Apply conflict resolution with debug output
    long *lkg_type;
    lkg_type = lvector(1, num_hbonds);
    
    // Save original distances
    double *hb_dist_original;
    hb_dist_original = dvector(1, num_hbonds);
    for (long k = 1; k <= num_hbonds; k++) {
        hb_dist_original[k] = hb_dist[k];
    }
    
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "CONFLICT RESOLUTION (hb_atompair_debug)\n");
    fprintf(stdout, "========================================\n");
    hb_atompair_debug(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, &Gvars.misc_pars);
    
    fprintf(stdout, "\n========================================\n");
    fprintf(stdout, "FINAL STATE (after conflict resolution)\n");
    fprintf(stdout, "========================================\n");
    long num_conflicts = 0;
    long num_kept = 0;
    for (long k = 1; k <= num_hbonds; k++) {
        int is_conflict = (hb_dist[k] < 0);
        if (is_conflict) {
            num_conflicts++;
            fprintf(stdout, "  %ld. %s -> %s, dist=%.6f (CONFLICT), lkg=%ld\n", 
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k], lkg_type[k]);
        } else {
            num_kept++;
            fprintf(stdout, "  %ld. %s -> %s, dist=%.6f (positive), lkg=%ld\n", 
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k], lkg_type[k]);
        }
    }
    
    fprintf(stdout, "\nSummary:\n");
    fprintf(stdout, "  Total H-bonds: %ld\n", num_hbonds);
    fprintf(stdout, "  Kept (positive distance): %ld\n", num_kept);
    fprintf(stdout, "  Conflicts (negative distance): %ld\n", num_conflicts);
    
cleanup:
    free_cmatrix(hb_atom1, 1, BUF512, 0, 4);
    free_cmatrix(hb_atom2, 1, BUF512, 0, 4);
    free_dvector(hb_dist, 1, BUF512);
    if (num_hbonds > 0) {
        free_dvector(hb_dist_original, 1, num_hbonds);
        free_lvector(lkg_type, 1, num_hbonds);
    }
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    
    if (output_fp != original_stdout) {
        fclose(output_fp);
    }
    stdout = original_stdout;
    
    if (argc >= 5) {
        fprintf(original_stdout, "Debug output written to: %s\n", argv[4]);
    }
    
    return 0;
}

