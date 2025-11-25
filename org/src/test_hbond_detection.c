#include "x3dna.h"
#include "json_writer.h"

// Test tool to isolate H-bond detection for debugging
// Usage: test_hbond_detection <pdb_file> <residue_i> <residue_j> [output.json]

int main(int argc, char *argv[]) {
    char *pdbfile;
    long i, j, num, num_residue;
    char **AtomName, **ResName, **Miscs;
    char *ChainID;
    long *ResSeq, *idx, *RY, **seidx;
    double **xyz;
    char *bseq;
    char hb_info[BUF512];
    FILE *json_fp = stdout;
    
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <pdb_file> <residue_i> <residue_j> [output.json]\n", argv[0]);
        fprintf(stderr, "  residue_i, residue_j: 1-based legacy residue indices\n");
        return 1;
    }
    
    pdbfile = argv[1];
    i = atol(argv[2]);
    j = atol(argv[3]);
    
    if (argc >= 5) {
        json_fp = open_file(argv[4], "w");
    }
    
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
        fprintf(stderr, "  Requested: i=%ld, j=%ld\n", i, j);
        return 1;
    }
    
    if (RY[i] < 0 || RY[j] < 0) {
        fprintf(stderr, "Error: One or both residues are not nucleotides\n");
        fprintf(stderr, "  Residue %ld: RY=%ld, base=%c\n", i, RY[i], bseq[i]);
        fprintf(stderr, "  Residue %ld: RY=%ld, base=%c\n", j, RY[j], bseq[j]);
        return 1;
    }
    
    // Calculate base frames (needed for some H-bond functions)
    char BDIR[BUF512];
    double **orien, **org;
    get_BDIR(BDIR, "Atomic_A.pdb");
    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq, Miscs,
               xyz, BDIR, orien, org);
    
    // Get base types
    char basei = bseq[i];
    char basej = bseq[j];
    
    // Call get_hbond_ij to get H-bond information
    get_hbond_ij(i, j, basei, basej, &Gvars.misc_pars, seidx, idx, AtomName, xyz, hb_info);
    
    // Output JSON with detailed information
    fprintf(json_fp, "{\n");
    fprintf(json_fp, "  \"pdb_file\": \"%s\",\n", pdbfile);
    fprintf(json_fp, "  \"residue_i\": %ld,\n", i);
    fprintf(json_fp, "  \"residue_j\": %ld,\n", j);
    fprintf(json_fp, "  \"base_i\": \"%c\",\n", basei);
    fprintf(json_fp, "  \"base_j\": \"%c\",\n", basej);
    fprintf(json_fp, "  \"hb_info\": \"%s\",\n", hb_info);
    
    // Parse hb_info to extract individual H-bonds
    // Format: [num] atom1typeatom2 dist ...
    long num_hb = 0;
    if (sscanf(hb_info, "[%ld]", &num_hb) == 1) {
        fprintf(json_fp, "  \"num_hbonds\": %ld,\n", num_hb);
        fprintf(json_fp, "  \"hbonds\": [\n");
        
        if (num_hb > 0) {
            char *pchar = strchr(hb_info, ' ');
            long k;
            for (k = 1; k <= num_hb && pchar != NULL; k++) {
                char atom1[10], atom2[10], hb_type_char;
                double dist;
                
                // Parse: " atom1typeatom2 dist"
                if (sscanf(pchar, " %4s%c%4s %lf", atom1, &hb_type_char, atom2, &dist) == 4) {
                    fprintf(json_fp, "    {\n");
                    fprintf(json_fp, "      \"index\": %ld,\n", k);
                    fprintf(json_fp, "      \"donor_atom\": \"%s\",\n", atom1);
                    fprintf(json_fp, "      \"acceptor_atom\": \"%s\",\n", atom2);
                    fprintf(json_fp, "      \"type\": \"%c\",\n", hb_type_char);
                    fprintf(json_fp, "      \"distance\": %.6f", dist);
                    
                    // Check if this is a "good" H-bond for adjust_pairQuality
                    int is_good = 0;
                    if (hb_type_char == '-' && dist >= 2.5 && dist <= 3.5) {
                        is_good = 1;
                    }
                    fprintf(json_fp, ",\n      \"is_good_for_quality\": %d", is_good);
                    
                    fprintf(json_fp, "\n    }");
                    if (k < num_hb) {
                        fprintf(json_fp, ",");
                    }
                    fprintf(json_fp, "\n");
                }
                
                // Move to next H-bond (format: " atom1typeatom2 dist" = ~15 chars)
                pchar += 15;
                while (*pchar == ' ') pchar++;
            }
        }
        
        fprintf(json_fp, "  ]\n");
    } else {
        fprintf(json_fp, "  \"num_hbonds\": 0,\n");
        fprintf(json_fp, "  \"hbonds\": []\n");
    }
    
    fprintf(json_fp, "}\n");
    
    // Cleanup
    if (json_fp != stdout) {
        close_file(json_fp);
    }
    
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
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    
    clear_my_globals();
    
    return 0;
}

