#include "x3dna.h"
#include "json_writer.h"

// Detailed debug test for N3->N2 H-bond type issue
// This tests the specific case where N3->N2 should be type='-' but shows as type='*'
// Usage: test_n3_n2_debug <pdb_file> <residue_i> <residue_j>

int main(int argc, char *argv[]) {
    char *pdbfile;
    long i, j, num, num_residue;
    char **AtomName, **ResName, **Miscs;
    char *ChainID;
    long *ResSeq, *idx, *RY, **seidx;
    double **xyz;
    char *bseq;
    char hb_info[BUF512];
    char *hb_type, **hb_atom1, **hb_atom2;
    double *hb_dist;
    long k, m, n, num_hbonds = 0, *lkg_type;
    
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <pdb_file> <residue_i> <residue_j>\n", argv[0]);
        fprintf(stderr, "  residue_i, residue_j: 1-based legacy residue indices\n");
        fprintf(stderr, "Example: %s data/pdb/3G8T.pdb 92 160\n", argv[0]);
        return 1;
    }
    
    pdbfile = argv[1];
    i = atol(argv[2]);
    j = atol(argv[3]);
    
    // Initialize globals
    set_my_globals(argv[0]);
    
    // Count atoms
    num = number_of_atoms(pdbfile, TRUE, "*");
    if (num <= 0) {
        fprintf(stderr, "Error: No atoms found in %s\n", pdbfile);
        return 1;
    }
    
    // Allocate arrays
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    
    // Read PDB
    read_pdb(pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    
    // Get residue information
    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);
    
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    RY = lvector(1, num_residue);
    bseq = cvector(1, num_residue);
    
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);
    
    // Validate residue indices
    if (i < 1 || i > num_residue || j < 1 || j > num_residue) {
        fprintf(stderr, "Error: Residue indices out of range [1, %ld]\n", num_residue);
        return 1;
    }
    
    if (RY[i] < 0 || RY[j] < 0) {
        fprintf(stderr, "Error: One or both residues are not nucleotides\n");
        return 1;
    }
    
    char basei = bseq[i];
    char basej = bseq[j];
    
    fprintf(stdout, "=== N3->N2 H-bond Debug Test ===\n\n");
    fprintf(stdout, "Residue %ld: base=%c, RY=%ld\n", i, basei, RY[i]);
    fprintf(stdout, "Residue %ld: base=%c, RY=%ld\n\n", j, basej, RY[j]);
    
    // Allocate H-bond arrays
    hb_atom1 = cmatrix(1, BUF512, 0, 4);
    hb_atom2 = cmatrix(1, BUF512, 0, 4);
    hb_dist = dvector(1, BUF512);
    
    // Step 1: Find all potential H-bonds (matches get_hbond_ij initial loop)
    fprintf(stdout, "=== Step 1: Finding potential H-bonds ===\n");
    for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
        for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
            if (good_hbatoms(&Gvars.misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]) &&
                within_limits(xyz[n], xyz[m], Gvars.misc_pars.hb_lower, Gvars.misc_pars.hb_dist1)) {
                if (++num_hbonds > BUF512)
                    fatal("Too many possible H-bonds between two bases\n");
                strcpy(hb_atom1[num_hbonds], AtomName[m]);
                strcpy(hb_atom2[num_hbonds], AtomName[n]);
                hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m]);
                
                // Check if this is N3->N2 or N2->N3
                int is_n3_n2 = 0;
                if (strncmp(AtomName[m], " N3 ", 4) == 0 && strncmp(AtomName[n], " N2 ", 4) == 0) {
                    is_n3_n2 = 1;
                } else if (strncmp(AtomName[m], " N2 ", 4) == 0 && strncmp(AtomName[n], " N3 ", 4) == 0) {
                    is_n3_n2 = 1;
                }
                
                if (is_n3_n2) {
                    fprintf(stdout, "  Found N3/N2 H-bond #%ld: %s -> %s, dist=%.6f\n", 
                            num_hbonds, AtomName[m], AtomName[n], hb_dist[num_hbonds]);
                    fprintf(stdout, "    Atom1 (from residue %ld): %s\n", i, AtomName[m]);
                    fprintf(stdout, "    Atom2 (from residue %ld): %s\n", j, AtomName[n]);
                    fprintf(stdout, "    Base1: %c, Base2: %c\n", basei, basej);
                }
            }
        }
    }
    fprintf(stdout, "Total potential H-bonds found: %ld\n\n", num_hbonds);
    
    if (num_hbonds == 0) {
        fprintf(stdout, "No H-bonds found. Exiting.\n");
        goto cleanup;
    }
    
    // Step 2: Resolve conflicts (hb_atompair)
    fprintf(stdout, "=== Step 2: Conflict Resolution (hb_atompair) ===\n");
    lkg_type = lvector(1, num_hbonds);
    
    // Show distances before conflict resolution
    for (k = 1; k <= num_hbonds; k++) {
        int is_n3_n2 = 0;
        if ((strncmp(hb_atom1[k], " N3 ", 4) == 0 && strncmp(hb_atom2[k], " N2 ", 4) == 0) ||
            (strncmp(hb_atom1[k], " N2 ", 4) == 0 && strncmp(hb_atom2[k], " N3 ", 4) == 0)) {
            is_n3_n2 = 1;
        }
        if (is_n3_n2) {
            fprintf(stdout, "  Before conflict resolution: H-bond #%ld: %s -> %s, dist=%.6f\n",
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k]);
        }
    }
    
    hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, &Gvars.misc_pars);
    
    // Show distances after conflict resolution (negative = conflict)
    fprintf(stdout, "After conflict resolution:\n");
    for (k = 1; k <= num_hbonds; k++) {
        int is_n3_n2 = 0;
        if ((strncmp(hb_atom1[k], " N3 ", 4) == 0 && strncmp(hb_atom2[k], " N2 ", 4) == 0) ||
            (strncmp(hb_atom1[k], " N2 ", 4) == 0 && strncmp(hb_atom2[k], " N3 ", 4) == 0)) {
            is_n3_n2 = 1;
        }
        if (is_n3_n2) {
            fprintf(stdout, "  H-bond #%ld: %s -> %s, dist=%.6f, lkg_type=%ld\n",
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k], lkg_type[k]);
            if (hb_dist[k] < 0.0) {
                fprintf(stdout, "    -> Marked as CONFLICT (negative distance)\n");
            } else {
                fprintf(stdout, "    -> NOT a conflict (positive distance)\n");
            }
        }
    }
    fprintf(stdout, "\n");
    
    // Step 3: Validate H-bonds (validate_hbonds)
    fprintf(stdout, "=== Step 3: H-bond Validation (validate_hbonds) ===\n");
    hb_type = cvector(1, num_hbonds);
    
    // Test donor_acceptor directly for N3->N2 before validation
    for (k = 1; k <= num_hbonds; k++) {
        int is_n3_n2 = 0;
        if ((strncmp(hb_atom1[k], " N3 ", 4) == 0 && strncmp(hb_atom2[k], " N2 ", 4) == 0) ||
            (strncmp(hb_atom1[k], " N2 ", 4) == 0 && strncmp(hb_atom2[k], " N3 ", 4) == 0)) {
            is_n3_n2 = 1;
        }
        if (is_n3_n2) {
            fprintf(stdout, "Testing donor_acceptor for H-bond #%ld:\n", k);
            fprintf(stdout, "  Input: basei=%c, basej=%c, atom1=%s, atom2=%s\n",
                    basei, basej, hb_atom1[k], hb_atom2[k]);
            
            // Test both orders
            char type1 = donor_acceptor(basei, basej, hb_atom1[k], hb_atom2[k]);
            char type2 = donor_acceptor(basej, basei, hb_atom2[k], hb_atom1[k]);
            char type3 = donor_acceptor(basei, basej, hb_atom2[k], hb_atom1[k]);
            char type4 = donor_acceptor(basej, basei, hb_atom1[k], hb_atom2[k]);
            
            fprintf(stdout, "  donor_acceptor(%c, %c, \"%s\", \"%s\") = '%c'\n",
                    basei, basej, hb_atom1[k], hb_atom2[k], type1);
            fprintf(stdout, "  donor_acceptor(%c, %c, \"%s\", \"%s\") = '%c'\n",
                    basej, basei, hb_atom2[k], hb_atom1[k], type2);
            fprintf(stdout, "  donor_acceptor(%c, %c, \"%s\", \"%s\") = '%c'\n",
                    basei, basej, hb_atom2[k], hb_atom1[k], type3);
            fprintf(stdout, "  donor_acceptor(%c, %c, \"%s\", \"%s\") = '%c'\n",
                    basej, basei, hb_atom1[k], hb_atom2[k], type4);
            fprintf(stdout, "\n");
        }
    }
    
    m = validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej, hb_atom1, hb_atom2);
    
    fprintf(stdout, "After validation (validate_hbonds returned %ld valid H-bonds):\n", m);
    for (k = 1; k <= num_hbonds; k++) {
        int is_n3_n2 = 0;
        if ((strncmp(hb_atom1[k], " N3 ", 4) == 0 && strncmp(hb_atom2[k], " N2 ", 4) == 0) ||
            (strncmp(hb_atom1[k], " N2 ", 4) == 0 && strncmp(hb_atom2[k], " N3 ", 4) == 0)) {
            is_n3_n2 = 1;
        }
        if (is_n3_n2) {
            fprintf(stdout, "  H-bond #%ld: %s -> %s, dist=%.6f, type='%c', lkg_type=%ld\n",
                    k, hb_atom1[k], hb_atom2[k], hb_dist[k], hb_type[k], lkg_type[k]);
            if (hb_type[k] == ' ') {
                fprintf(stdout, "    -> Type is ' ' (invalid/skipped)\n");
            } else if (hb_type[k] == '-') {
                fprintf(stdout, "    -> Type is '-' (standard H-bond) ✓\n");
            } else if (hb_type[k] == '*') {
                fprintf(stdout, "    -> Type is '*' (non-standard H-bond) ⚠\n");
            }
        }
    }
    fprintf(stdout, "\n");
    
    // Show final hb_info string
    get_hbond_ij(i, j, basei, basej, &Gvars.misc_pars, seidx, idx, AtomName, xyz, hb_info);
    fprintf(stdout, "=== Final hb_info string ===\n");
    fprintf(stdout, "%s\n\n", hb_info);
    
cleanup:
    free_cmatrix(AtomName, 1, num, 0, 4);
    free_cmatrix(ResName, 1, num, 0, 3);
    free_cvector(ChainID, 1, num);
    free_lvector(ResSeq, 1, num);
    free_dmatrix(xyz, 1, num, 1, 3);
    free_cmatrix(Miscs, 1, num, 0, NMISC);
    free_lvector(idx, 1, num);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_lvector(RY, 1, num_residue);
    free_cvector(bseq, 1, num_residue);
    if (num_hbonds > 0) {
        free_cmatrix(hb_atom1, 1, BUF512, 0, 4);
        free_cmatrix(hb_atom2, 1, BUF512, 0, 4);
        free_dvector(hb_dist, 1, BUF512);
        free_cvector(hb_type, 1, num_hbonds);
        free_lvector(lkg_type, 1, num_hbonds);
    }
    
    clear_my_globals();
    
    return 0;
}

