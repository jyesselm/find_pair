/**
 * @file generate_residue_ordering_json.c
 * @brief Generate JSON file with residue ordering from legacy code
 * 
 * This tool generates a JSON file containing the residue ordering information
 * using legacy's residue_idx() function, matching the format from modern code.
 * 
 * Usage: generate_residue_ordering_json <pdb_file> <output_json>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "x3dna.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <pdb_file> <output_json>\n", argv[0]);
        fprintf(stderr, "Example: %s data/pdb/3G8T.pdb data/residue_ordering_legacy/3G8T.json\n", argv[0]);
        return 1;
    }
    
    char *pdb_file = argv[1];
    char *output_json = argv[2];
    
    // Initialize
    set_my_globals(argv[0]);
    
    // Count atoms first
    long num = number_of_atoms(pdb_file, TRUE, (char*)"*");
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
    read_pdb(pdb_file, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, (char*)"*");
    
    // Get atom indices
    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);
    
    // Get residue indices (legacy's residue_idx function)
    long num_residue;
    long **seidx;
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    
    // Open output file
    FILE *out = fopen(output_json, "w");
    if (!out) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", output_json);
        return 1;
    }
    
    // Extract PDB ID from filename
    char *pdb_id = strrchr(pdb_file, '/');
    if (!pdb_id) {
        pdb_id = strrchr(pdb_file, '\\');
    }
    if (pdb_id) {
        pdb_id++;
    } else {
        pdb_id = pdb_file;
    }
    // Remove .pdb extension
    char pdb_id_clean[256];
    strncpy(pdb_id_clean, pdb_id, sizeof(pdb_id_clean) - 1);
    pdb_id_clean[sizeof(pdb_id_clean) - 1] = '\0';
    char *dot = strrchr(pdb_id_clean, '.');
    if (dot) {
        *dot = '\0';
    }
    
    // Write JSON header
    fprintf(out, "{\n");
    fprintf(out, "  \"pdb_id\": \"%s\",\n", pdb_id_clean);
    fprintf(out, "  \"total_residues\": %ld,\n", num_residue);
    fprintf(out, "  \"residues\": [\n");
    
    // Write each residue
    for (long i = 1; i <= num_residue; i++) {
        long start = seidx[i][1];
        long end = seidx[i][2];
        
        // Get residue info from first atom
        char *res_name = ResName[start];
        char chain_id = ChainID[start];
        long res_seq = ResSeq[start];
        char insertion = (Miscs == NULL) ? ' ' : Miscs[start][2];
        
        // Count atoms in this residue
        long num_atoms = end - start + 1;
        
        // Get first and last atom names
        char *first_atom = AtomName[start];
        char *last_atom = AtomName[end];
        
        // Write residue JSON
        fprintf(out, "    {\n");
        fprintf(out, "      \"legacy_index\": %ld,\n", i);
        fprintf(out, "      \"residue_name\": \"%s\",\n", res_name);
        fprintf(out, "      \"chain_id\": \"%c\",\n", chain_id);
        fprintf(out, "      \"residue_seq\": %ld,\n", res_seq);
        fprintf(out, "      \"insertion_code\": \"%c\",\n", insertion);
        fprintf(out, "      \"num_atoms\": %ld", num_atoms);
        
        if (num_atoms > 0) {
            fprintf(out, ",\n      \"first_atom\": \"%s\",\n", first_atom);
            fprintf(out, "      \"last_atom\": \"%s\"", last_atom);
        }
        
        if (i < num_residue) {
            fprintf(out, "\n    },\n");
        } else {
            fprintf(out, "\n    }\n");
        }
    }
    
    // Write JSON footer
    fprintf(out, "  ]\n");
    fprintf(out, "}\n");
    
    fclose(out);
    
    printf("Generated legacy residue ordering JSON: %s\n", output_json);
    printf("Total residues: %ld\n", num_residue);
    
    // Cleanup
    free_lvector(idx, 1, num);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cmatrix(AtomName, 1, num, 0, 4);
    free_cmatrix(ResName, 1, num, 0, 3);
    free_cvector(ChainID, 1, num);
    free_lvector(ResSeq, 1, num);
    free_dmatrix(xyz, 1, num, 1, 3);
    if (Miscs) {
        free_cmatrix(Miscs, 1, num, 0, NMISC);
    }
    
    return 0;
}

