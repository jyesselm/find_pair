#include "json_writer.h"
#include <sys/stat.h>
#include <sys/types.h>

/* Singleton state */
static long initialized = FALSE;
static FILE *json_file = NULL;
static char json_filename[BUF1K];
static char globals_filename[BUF1K];
static long first_entry = TRUE;

/* Helper: Escape string for JSON */
static void json_escape_string(const char *str, char *out, size_t out_size) {
    size_t i, j = 0;
    if (!str) {
        out[0] = '\0';
        return;
    }
    for (i = 0; str[i] != '\0' && j < out_size - 1; i++) {
        if (str[i] == '"' || str[i] == '\\') {
            if (j < out_size - 2) {
                out[j++] = '\\';
                out[j++] = str[i];
            }
        } else if (str[i] == '\n') {
            if (j < out_size - 3) {
                out[j++] = '\\';
                out[j++] = 'n';
            }
        } else {
            out[j++] = str[i];
        }
    }
    out[j] = '\0';
}

/* Helper: Write array of doubles */
static void json_write_double_array(FILE *fp, double *arr, long n) {
    long i;
    fprintf(fp, "[");
    for (i = 0; i < n; i++) {
        if (i > 0) fprintf(fp, ", ");
        if (arr[i] > EMPTY_CRITERION || arr[i] < -EMPTY_CRITERION) {
            fprintf(fp, "%.6f", arr[i]);
        } else {
            fprintf(fp, "null");
        }
    }
    fprintf(fp, "]");
}

/* Helper: Write 3x3 matrix */
static void json_write_matrix(FILE *fp, double **m) {
    long i, j;
    fprintf(fp, "[[");
    for (i = 1; i <= 3; i++) {
        if (i > 1) fprintf(fp, "], [");
        for (j = 1; j <= 3; j++) {
            if (j > 1) fprintf(fp, ", ");
            fprintf(fp, "%.6f", m[i][j]);
        }
    }
    fprintf(fp, "]]");
}

long json_writer_init(const char *pdbfile) {
    char pdb_name[BUF512];
    char dir_path[BUF1K];
    struct stat st = {0};
    
    if (initialized) {
        return TRUE; /* Already initialized */
    }
    
    /* Extract base name without extension */
    bname_noext((char *)pdbfile, pdb_name);
    
    /* Determine project root data directory */
    /* Check if ../data exists (we're in org/ directory) */
    if (stat("../data", &st) == 0) {
        sprintf(dir_path, "../data");
        sprintf(json_filename, "../data/json_legacy/%s.json", pdb_name);
        sprintf(globals_filename, "../data/json_legacy/%s_globals.json", pdb_name);
    } else if (stat("data", &st) == 0) {
        /* We're already at project root */
        sprintf(dir_path, "data");
        sprintf(json_filename, "data/json_legacy/%s.json", pdb_name);
        sprintf(globals_filename, "data/json_legacy/%s_globals.json", pdb_name);
    } else {
        /* Create data directory at current location */
        sprintf(dir_path, "data");
        sprintf(json_filename, "data/json_legacy/%s.json", pdb_name);
        sprintf(globals_filename, "data/json_legacy/%s_globals.json", pdb_name);
    }
    
    /* Create parent directory 'data' if it doesn't exist */
    if (stat(dir_path, &st) == -1) {
        if (mkdir(dir_path, 0755) != 0) {
            fprintf(stderr, "[JSON_WRITER] Warning: Could not create directory %s\n", dir_path);
            return FALSE;
        }
    }
    /* Create json_legacy directory */
    sprintf(dir_path, "%s/json_legacy", dir_path);
    if (stat(dir_path, &st) == -1) {
        if (mkdir(dir_path, 0755) != 0) {
            fprintf(stderr, "[JSON_WRITER] Warning: Could not create directory %s\n", dir_path);
            return FALSE;
        }
    }
    
    /* Create JSON filename with corrected path */
    if (stat("../data", &st) == 0) {
        sprintf(json_filename, "../data/json_legacy/%s.json", pdb_name);
        sprintf(globals_filename, "../data/json_legacy/%s_globals.json", pdb_name);
    } else {
        sprintf(json_filename, "data/json_legacy/%s.json", pdb_name);
        sprintf(globals_filename, "data/json_legacy/%s_globals.json", pdb_name);
    }
    
    /* Open file for writing */
    json_file = fopen(json_filename, "w");
    if (!json_file) {
        fprintf(stderr, "[JSON_WRITER] Error: Could not open %s for writing\n", json_filename);
        return FALSE;
    }
    
    /* Write JSON header */
    fprintf(json_file, "{\n");
    fprintf(json_file, "  \"pdb_file\": \"%s\",\n", pdbfile);
    fprintf(json_file, "  \"pdb_name\": \"%s\",\n", pdb_name);
    fprintf(json_file, "  \"calculations\": [\n");
    
    first_entry = TRUE;
    initialized = TRUE;
    
    fprintf(stderr, "[JSON_WRITER] Initialized: %s\n", json_filename);
    return TRUE;
}

void json_writer_finalize(void) {
    if (!initialized || !json_file) {
        return;
    }
    
    /* Close calculations array */
    fprintf(json_file, "\n  ],\n");
    fprintf(json_file, "  \"metadata\": {\n");
    fprintf(json_file, "    \"version\": \"%s\"\n", Gvars.X3DNA_VER);
    fprintf(json_file, "  }\n");
    fprintf(json_file, "}\n");
    
    fclose(json_file);
    json_file = NULL;
    initialized = FALSE;
    
    fprintf(stderr, "[JSON_WRITER] Finalized: %s\n", json_filename);
}

long json_writer_is_initialized(void) {
    return initialized && (json_file != NULL);
}

void json_writer_record_bpstep_params(long bp_idx1, long bp_idx2, 
                                      double *pars,
                                      double *mst_org,
                                      double **mst_orien) {
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"bpstep_params\",\n");
    fprintf(json_file, "      \"bp_idx1\": %ld,\n", bp_idx1);
    fprintf(json_file, "      \"bp_idx2\": %ld,\n", bp_idx2);
    fprintf(json_file, "      \"params\": {\n");
    fprintf(json_file, "        \"Shift\": %.6f,\n", pars[1]);
    fprintf(json_file, "        \"Slide\": %.6f,\n", pars[2]);
    fprintf(json_file, "        \"Rise\": %.6f,\n", pars[3]);
    fprintf(json_file, "        \"Tilt\": %.6f,\n", pars[4]);
    fprintf(json_file, "        \"Roll\": %.6f,\n", pars[5]);
    fprintf(json_file, "        \"Twist\": %.6f\n", pars[6]);
    fprintf(json_file, "      },\n");
    fprintf(json_file, "      \"mst_org\": ");
    json_write_double_array(json_file, &mst_org[1], 3);
    fprintf(json_file, ",\n");
    fprintf(json_file, "      \"mst_orien\": ");
    json_write_matrix(json_file, mst_orien);
    fprintf(json_file, "\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_helical_params(long bp_idx1, long bp_idx2,
                                       double *pars,
                                       double *mst_orgH,
                                       double **mst_orienH) {
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"helical_params\",\n");
    fprintf(json_file, "      \"bp_idx1\": %ld,\n", bp_idx1);
    fprintf(json_file, "      \"bp_idx2\": %ld,\n", bp_idx2);
    fprintf(json_file, "      \"params\": ");
    json_write_double_array(json_file, &pars[1], 6);
    fprintf(json_file, ",\n");
    fprintf(json_file, "      \"mst_orgH\": ");
    json_write_double_array(json_file, &mst_orgH[1], 3);
    fprintf(json_file, ",\n");
    fprintf(json_file, "      \"mst_orienH\": ");
    json_write_matrix(json_file, mst_orienH);
    fprintf(json_file, "\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_base_pair(long i, long j, char *bp_type,
                                  double *dir_xyz,
                                  double **orien_i,
                                  double **orien_j,
                                  double *org_i,
                                  double *org_j) {
    char esc_bp_type[BUF32];
    
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    json_escape_string(bp_type, esc_bp_type, sizeof(esc_bp_type));
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"base_pair\",\n");
    fprintf(json_file, "      \"base_i\": %ld,\n", i);
    fprintf(json_file, "      \"base_j\": %ld,\n", j);
    if (bp_type) {
        fprintf(json_file, "      \"bp_type\": \"%s\",\n", esc_bp_type);
    }
    if (dir_xyz) {
        fprintf(json_file, "      \"dir_xyz\": [%.6f, %.6f, %.6f],\n",
                dir_xyz[1], dir_xyz[2], dir_xyz[3]);
    }
    fprintf(json_file, "      \"orien_i\": ");
    json_write_matrix(json_file, orien_i);
    fprintf(json_file, ",\n");
    fprintf(json_file, "      \"orien_j\": ");
    json_write_matrix(json_file, orien_j);
    fprintf(json_file, ",\n");
    fprintf(json_file, "      \"org_i\": ");
    json_write_double_array(json_file, &org_i[1], 3);
    fprintf(json_file, ",\n");
    fprintf(json_file, "      \"org_j\": ");
    json_write_double_array(json_file, &org_j[1], 3);
    fprintf(json_file, "\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_ref_frame(long residue_idx, double **orien, double *org) {
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"ref_frame\",\n");
    fprintf(json_file, "      \"residue_idx\": %ld,\n", residue_idx);
    fprintf(json_file, "      \"orien\": ");
    json_write_matrix(json_file, orien);
    fprintf(json_file, ",\n");
    fprintf(json_file, "      \"org\": ");
    json_write_double_array(json_file, &org[1], 3);
    fprintf(json_file, "\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_sequence(long num_residue, char *bseq) {
    char esc_seq[BUF512];
    
    if (!json_writer_is_initialized()) return;
    if (!bseq) return;
    
    json_escape_string(bseq, esc_seq, sizeof(esc_seq));
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"sequence\",\n");
    fprintf(json_file, "      \"num_residue\": %ld,\n", num_residue);
    fprintf(json_file, "      \"sequence\": \"%s\"\n", esc_seq);
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_bp_sequence(long num_bp, char **bp_seq, long ds) {
    long i, j;
    char esc_seq[BUF32];
    
    if (!json_writer_is_initialized()) return;
    if (!bp_seq) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"bp_sequence\",\n");
    fprintf(json_file, "      \"num_bp\": %ld,\n", num_bp);
    fprintf(json_file, "      \"ds\": %ld,\n", ds);
    fprintf(json_file, "      \"pairs\": [\n");
    
    for (i = 0; i <= ds; i++) {
        if (i > 0) fprintf(json_file, ",\n");
        fprintf(json_file, "        [");
        for (j = 1; j <= num_bp; j++) {
            if (j > 1) fprintf(json_file, ", ");
            if (bp_seq[i] && bp_seq[i][j]) {
                json_escape_string(&bp_seq[i][j], esc_seq, sizeof(esc_seq));
                fprintf(json_file, "\"%s\"", esc_seq);
            } else {
                fprintf(json_file, "\"\"");
            }
        }
        fprintf(json_file, "]");
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_pdb_atoms(long num, char **AtomName, char **ResName,
                                   char *ChainID, long *ResSeq, double **xyz,
                                   char **Miscs) {
    char esc_atom[BUF32], esc_res[BUF32];
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!AtomName || !ResName || !ChainID || !ResSeq || !xyz) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"pdb_atoms\",\n");
    fprintf(json_file, "      \"num_atoms\": %ld,\n", num);
    fprintf(json_file, "      \"atoms\": [\n");
    
    for (i = 1; i <= num; i++) {
        if (i > 1) fprintf(json_file, ",\n");
        fprintf(json_file, "        {\n");
        
        /* Atom index/position in array */
        fprintf(json_file, "          \"atom_idx\": %ld,\n", i);
        
        /* Atom name */
        json_escape_string(AtomName[i], esc_atom, sizeof(esc_atom));
        fprintf(json_file, "          \"atom_name\": \"%s\",\n", esc_atom);
        
        /* Residue name */
        json_escape_string(ResName[i], esc_res, sizeof(esc_res));
        fprintf(json_file, "          \"residue_name\": \"%s\",\n", esc_res);
        
        /* Chain ID */
        fprintf(json_file, "          \"chain_id\": \"%c\",\n", ChainID[i]);
        
        /* Residue sequence number */
        fprintf(json_file, "          \"residue_seq\": %ld,\n", ResSeq[i]);
        
        /* Coordinates */
        fprintf(json_file, "          \"xyz\": [%.6f, %.6f, %.6f]", 
                xyz[i][1], xyz[i][2], xyz[i][3]);
        
        /* Miscellaneous info if available */
        if (Miscs && Miscs[i]) {
            fprintf(json_file, ",\n          \"record_type\": \"%c\"", Miscs[i][0]);
            if (Miscs[i][1] != ' ') {
                fprintf(json_file, ",\n          \"alt_loc\": \"%c\"", Miscs[i][1]);
            }
            if (Miscs[i][2] != ' ') {
                fprintf(json_file, ",\n          \"insertion\": \"%c\"", Miscs[i][2]);
            }
        }
        
        fprintf(json_file, "\n        }");
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_residue_indices(long num_residue, long **seidx) {
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!seidx) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"residue_indices\",\n");
    fprintf(json_file, "      \"num_residue\": %ld,\n", num_residue);
    fprintf(json_file, "      \"seidx\": [\n");
    
    for (i = 1; i <= num_residue; i++) {
        if (i > 1) fprintf(json_file, ",\n");
        fprintf(json_file, "        {\"residue_idx\": %ld, \"start_atom\": %ld, \"end_atom\": %ld}",
                i, seidx[i][1], seidx[i][2]);
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_base_pairs(long ds, long num_bp, long **pair_num) {
    long i, j;
    
    if (!json_writer_is_initialized()) return;
    if (!pair_num) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"base_pairs\",\n");
    fprintf(json_file, "      \"ds\": %ld,\n", ds);
    fprintf(json_file, "      \"num_bp\": %ld,\n", num_bp);
    fprintf(json_file, "      \"pair_num\": [\n");
    
    for (i = 1; i <= ds + 1; i++) {
        if (i > 1) fprintf(json_file, ",\n");
        fprintf(json_file, "        [");
        for (j = 1; j <= num_bp; j++) {
            if (j > 1) fprintf(json_file, ", ");
            fprintf(json_file, "%ld", pair_num[i][j]);
        }
        fprintf(json_file, "]");
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_all_ref_frames(long ds, long num_bp, double **orien, double **org) {
    long i, j, k, idx;
    
    if (!json_writer_is_initialized()) return;
    if (!orien || !org) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"all_ref_frames\",\n");
    fprintf(json_file, "      \"ds\": %ld,\n", ds);
    fprintf(json_file, "      \"num_bp\": %ld,\n", num_bp);
    fprintf(json_file, "      \"frames\": [\n");
    
    for (i = 1; i <= ds; i++) {
        if (i > 1) fprintf(json_file, ",\n");
        fprintf(json_file, "        {\"strand\": %ld, \"bp_frames\": [\n", i);
        
        for (j = 1; j <= num_bp; j++) {
            if (j > 1) fprintf(json_file, ",\n");
            fprintf(json_file, "          {\"bp_idx\": %ld, \"orien\": ", j);
            
            /* Write 3x3 matrix */
            fprintf(json_file, "[[");
            for (k = 1; k <= 3; k++) {
                if (k > 1) fprintf(json_file, "], [");
                idx = (j - 1) * 9 + (k - 1) * 3;
                fprintf(json_file, "%.6f, %.6f, %.6f",
                        orien[i][idx + 1], orien[i][idx + 2], orien[i][idx + 3]);
            }
            fprintf(json_file, "]], ");
            
            /* Write origin */
            idx = (j - 1) * 3;
            fprintf(json_file, "\"org\": [%.6f, %.6f, %.6f]",
                    org[i][idx + 1], org[i][idx + 2], org[i][idx + 3]);
            fprintf(json_file, "}");
        }
        
        fprintf(json_file, "\n        ]}");
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_wc_info(long num_bp, long *WC_info) {
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!WC_info) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"wc_info\",\n");
    fprintf(json_file, "      \"num_bp\": %ld,\n", num_bp);
    fprintf(json_file, "      \"wc_values\": [");
    
    for (i = 1; i <= num_bp; i++) {
        if (i > 1) fprintf(json_file, ", ");
        fprintf(json_file, "%ld", WC_info[i]);
    }
    
    fprintf(json_file, "]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_ry(long num_residue, long *RY) {
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!RY) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"ry_classification\",\n");
    fprintf(json_file, "      \"num_residue\": %ld,\n", num_residue);
    fprintf(json_file, "      \"ry_values\": [");
    
    for (i = 1; i <= num_residue; i++) {
        if (i > 1) fprintf(json_file, ", ");
        fprintf(json_file, "%ld", RY[i]);
    }
    
    fprintf(json_file, "]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_helices(long num_bp, long *bphlx) {
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!bphlx) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"helices\",\n");
    fprintf(json_file, "      \"num_bp\": %ld,\n", num_bp);
    fprintf(json_file, "      \"bphlx\": [");
    
    for (i = 1; i <= num_bp; i++) {
        if (i > 1) fprintf(json_file, ", ");
        fprintf(json_file, "%ld", bphlx[i]);
    }
    
    fprintf(json_file, "]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_twist_rise(long nbpm1, double **twist_rise) {
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!twist_rise) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"twist_rise\",\n");
    fprintf(json_file, "      \"num_steps\": %ld,\n", nbpm1);
    fprintf(json_file, "      \"steps\": [\n");
    
    for (i = 1; i <= nbpm1; i++) {
        if (i > 1) fprintf(json_file, ",\n");
        fprintf(json_file, "        {\"step_idx\": %ld, \"twist\": %.6f, \"rise\": %.6f}",
                i, twist_rise[i][1], twist_rise[i][2]);
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_input_parameters(miscPars *misc_pars, long ds, long hetatm, long ip) {
    if (!json_writer_is_initialized()) return;
    if (!misc_pars) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"input_parameters\",\n");
    fprintf(json_file, "      \"ds\": %ld,\n", ds);
    fprintf(json_file, "      \"hetatm\": %ld,\n", hetatm);
    fprintf(json_file, "      \"ip\": %ld\n", ip);
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_global_variables(void) {
    FILE *globals_file;
    long i;
    char esc_str[BUF512];
    
    if (!initialized) return;
    
    /* Open globals file for writing */
    globals_file = fopen(globals_filename, "w");
    if (!globals_file) {
        fprintf(stderr, "[JSON_WRITER] Warning: Could not open %s for writing\n", globals_filename);
        return;
    }
    
    /* Write global variables JSON */
    fprintf(globals_file, "{\n");
    fprintf(globals_file, "  \"global_variables\": {\n");
    
    /* Gvars structure values */
    fprintf(globals_file, "    \"DEBUG\": %ld,\n", Gvars.DEBUG);
    fprintf(globals_file, "    \"VERBOSE\": %ld,\n", Gvars.VERBOSE);
    fprintf(globals_file, "    \"NUM_ELE\": %ld,\n", Gvars.NUM_ELE);
    fprintf(globals_file, "    \"CHAIN_CASE\": %ld,\n", Gvars.CHAIN_CASE);
    fprintf(globals_file, "    \"ALL_MODEL\": %ld,\n", Gvars.ALL_MODEL);
    fprintf(globals_file, "    \"ATTACH_RESIDUE\": %ld,\n", Gvars.ATTACH_RESIDUE);
    fprintf(globals_file, "    \"THREE_LETTER_NTS\": %ld,\n", Gvars.THREE_LETTER_NTS);
    fprintf(globals_file, "    \"PDBV3\": %ld,\n", Gvars.PDBV3);
    fprintf(globals_file, "    \"ORIGINAL_COORDINATE\": %ld,\n", Gvars.ORIGINAL_COORDINATE);
    fprintf(globals_file, "    \"OCCUPANCY\": %ld,\n", Gvars.OCCUPANCY);
    fprintf(globals_file, "    \"HEADER\": %ld,\n", Gvars.HEADER);
    fprintf(globals_file, "    \"mmcif\": %ld,\n", Gvars.mmcif);
    fprintf(globals_file, "    \"NT_CUTOFF\": %.6f,\n", Gvars.NT_CUTOFF);
    
    json_escape_string(Gvars.X3DNA_VER, esc_str, sizeof(esc_str));
    fprintf(globals_file, "    \"X3DNA_VER\": \"%s\",\n", esc_str);
    
    json_escape_string(Gvars.X3DNA_HOMEDIR, esc_str, sizeof(esc_str));
    fprintf(globals_file, "    \"X3DNA_HOMEDIR\": \"%s\",\n", esc_str);
    
    json_escape_string(Gvars.CHAIN_MARKERS, esc_str, sizeof(esc_str));
    fprintf(globals_file, "    \"CHAIN_MARKERS\": \"%s\",\n", esc_str);
    
    json_escape_string(Gvars.REBUILD_CHAIN_IDS, esc_str, sizeof(esc_str));
    fprintf(globals_file, "    \"REBUILD_CHAIN_IDS\": \"%s\",\n", esc_str);
    
    if (Gvars.PROGNAME) {
        json_escape_string(Gvars.PROGNAME, esc_str, sizeof(esc_str));
        fprintf(globals_file, "    \"PROGNAME\": \"%s\",\n", esc_str);
    } else {
        fprintf(globals_file, "    \"PROGNAME\": \"\",\n");
    }
    
    fprintf(globals_file, "    \"NUM_SATOM\": %ld,\n", Gvars.NUM_SATOM);
    fprintf(globals_file, "    \"NUM_SBASE\": %ld,\n", Gvars.NUM_SBASE);
    fprintf(globals_file, "    \"Name0\": %ld,\n", Gvars.Name0);
    fprintf(globals_file, "    \"label_RC8_YC6\": %ld,\n", Gvars.label_RC8_YC6);
    
    /* ATOMLIST */
    fprintf(globals_file, "    \"ATOMLIST\": [\n");
    if (Gvars.ATOMLIST) {
        for (i = 1; i <= Gvars.NUM_SATOM; i++) {
            if (i > 1) fprintf(globals_file, ",\n");
            json_escape_string(Gvars.ATOMLIST[i], esc_str, sizeof(esc_str));
            fprintf(globals_file, "      \"%s\"", esc_str);
        }
    }
    fprintf(globals_file, "\n    ],\n");
    
    /* BASELIST */
    fprintf(globals_file, "    \"BASELIST\": [\n");
    if (Gvars.BASELIST) {
        for (i = 1; i <= Gvars.NUM_SBASE; i++) {
            if (i > 1) fprintf(globals_file, ",\n");
            json_escape_string(Gvars.BASELIST[i], esc_str, sizeof(esc_str));
            fprintf(globals_file, "      \"%s\"", esc_str);
        }
    }
    fprintf(globals_file, "\n    ],\n");
    
    /* misc_pars */
    fprintf(globals_file, "    \"misc_pars\": {\n");
    fprintf(globals_file, "      \"min_base_hb\": %ld,\n", Gvars.misc_pars.min_base_hb);
    fprintf(globals_file, "      \"hb_lower\": %.6f,\n", Gvars.misc_pars.hb_lower);
    fprintf(globals_file, "      \"hb_dist1\": %.6f,\n", Gvars.misc_pars.hb_dist1);
    fprintf(globals_file, "      \"hb_dist2\": %.6f,\n", Gvars.misc_pars.hb_dist2);
    fprintf(globals_file, "      \"max_dorg\": %.6f,\n", Gvars.misc_pars.max_dorg);
    fprintf(globals_file, "      \"min_dorg\": %.6f,\n", Gvars.misc_pars.min_dorg);
    fprintf(globals_file, "      \"max_dv\": %.6f,\n", Gvars.misc_pars.max_dv);
    fprintf(globals_file, "      \"min_dv\": %.6f,\n", Gvars.misc_pars.min_dv);
    fprintf(globals_file, "      \"max_plane_angle\": %.6f,\n", Gvars.misc_pars.max_plane_angle);
    fprintf(globals_file, "      \"min_plane_angle\": %.6f,\n", Gvars.misc_pars.min_plane_angle);
    fprintf(globals_file, "      \"max_dNN\": %.6f,\n", Gvars.misc_pars.max_dNN);
    fprintf(globals_file, "      \"min_dNN\": %.6f,\n", Gvars.misc_pars.min_dNN);
    fprintf(globals_file, "      \"helix_break\": %.6f,\n", Gvars.misc_pars.helix_break);
    fprintf(globals_file, "      \"std_curved\": %.6f,\n", Gvars.misc_pars.std_curved);
    fprintf(globals_file, "      \"water_dist\": %.6f,\n", Gvars.misc_pars.water_dist);
    fprintf(globals_file, "      \"water_dlow\": %.6f,\n", Gvars.misc_pars.water_dlow);
    fprintf(globals_file, "      \"o3p_dist\": %.6f,\n", Gvars.misc_pars.o3p_dist);
    
    json_escape_string(Gvars.misc_pars.alt_list, esc_str, sizeof(esc_str));
    fprintf(globals_file, "      \"alt_list\": \"%s\",\n", esc_str);
    
    json_escape_string(Gvars.misc_pars.hb_atoms, esc_str, sizeof(esc_str));
    fprintf(globals_file, "      \"hb_atoms\": \"%s\",\n", esc_str);
    
    json_escape_string(Gvars.misc_pars.water_atoms, esc_str, sizeof(esc_str));
    fprintf(globals_file, "      \"water_atoms\": \"%s\"\n", esc_str);
    fprintf(globals_file, "    }\n");
    
    fprintf(globals_file, "  },\n");
    
    /* Constants */
    fprintf(globals_file, "  \"constants\": {\n");
    fprintf(globals_file, "    \"NR_END\": %ld,\n", (long)NR_END);
    fprintf(globals_file, "    \"TRUE\": %ld,\n", (long)TRUE);
    fprintf(globals_file, "    \"FALSE\": %ld,\n", (long)FALSE);
    fprintf(globals_file, "    \"BUF32\": %d,\n", BUF32);
    fprintf(globals_file, "    \"BUF512\": %d,\n", BUF512);
    fprintf(globals_file, "    \"BUF1K\": %d,\n", BUF1K);
    fprintf(globals_file, "    \"BUF2K\": %d,\n", BUF2K);
    fprintf(globals_file, "    \"BUFBIG\": %d,\n", BUFBIG);
    fprintf(globals_file, "    \"PI\": %.15f,\n", PI);
    fprintf(globals_file, "    \"XEPS\": %.10e,\n", XEPS);
    fprintf(globals_file, "    \"XBIG\": %.10e,\n", XBIG);
    fprintf(globals_file, "    \"XBIG_CUTOFF\": %.10e,\n", XBIG_CUTOFF);
    fprintf(globals_file, "    \"MFACTOR\": %.6f,\n", MFACTOR);
    fprintf(globals_file, "    \"NMISC\": %d,\n", NMISC);
    fprintf(globals_file, "    \"DEBUG_LEVEL\": %d,\n", DEBUG_LEVEL);
    fprintf(globals_file, "    \"EMPTY_CRITERION\": %ld,\n", (long)EMPTY_CRITERION);
    fprintf(globals_file, "    \"EMPTY_NUMBER\": %.6f\n", EMPTY_NUMBER);
    fprintf(globals_file, "  }\n");
    
    fprintf(globals_file, "}\n");
    
    fclose(globals_file);
    fprintf(stderr, "[JSON_WRITER] Global variables and constants saved to: %s\n", globals_filename);
}

void json_writer_record_hbonds(long base_i, long base_j, long num_hbonds,
                                char **hb_atom1, char **hb_atom2,
                                double *hb_dist, char *hb_type,
                                long *lkg_type) {
    long i;
    char esc_atom1[BUF512], esc_atom2[BUF512];
    
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"hydrogen_bonds\",\n");
    fprintf(json_file, "      \"base_i\": %ld,\n", base_i);
    fprintf(json_file, "      \"base_j\": %ld,\n", base_j);
    fprintf(json_file, "      \"num_hbonds\": %ld,\n", num_hbonds);
    fprintf(json_file, "      \"hbonds\": [\n");
    
    for (i = 1; i <= num_hbonds; i++) {
        if (i > 1) fprintf(json_file, ",\n");
        json_escape_string(hb_atom1[i], esc_atom1, sizeof(esc_atom1));
        json_escape_string(hb_atom2[i], esc_atom2, sizeof(esc_atom2));
        fprintf(json_file, "        {\n");
        fprintf(json_file, "          \"hbond_idx\": %ld,\n", i);
        fprintf(json_file, "          \"donor_atom\": \"%s\",\n", esc_atom1);
        fprintf(json_file, "          \"acceptor_atom\": \"%s\",\n", esc_atom2);
        fprintf(json_file, "          \"distance\": %.6f,\n", hb_dist ? fabs(hb_dist[i]) : 0.0);
        fprintf(json_file, "          \"type\": \"%c\",\n", hb_type && hb_type[i] ? hb_type[i] : ' ');
        fprintf(json_file, "          \"linkage_type\": %ld\n", lkg_type ? lkg_type[i] : 0);
        fprintf(json_file, "        }");
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_base_frame_calc(long residue_idx, char base_type,
                                         const char *standard_template,
                                         double rms_fit, long num_matched,
                                         char **matched_atoms, long num_atoms) {
    long i;
    char esc_atom[BUF512], esc_template[BUF512];
    
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"base_frame_calc\",\n");
    fprintf(json_file, "      \"residue_idx\": %ld,\n", residue_idx);
    fprintf(json_file, "      \"base_type\": \"%c\",\n", base_type);
    json_escape_string(standard_template, esc_template, sizeof(esc_template));
    fprintf(json_file, "      \"standard_template\": \"%s\",\n", esc_template);
    fprintf(json_file, "      \"rms_fit\": %.6f,\n", rms_fit);
    fprintf(json_file, "      \"num_matched_atoms\": %ld,\n", num_matched);
    fprintf(json_file, "      \"matched_atoms\": [");
    
    for (i = 1; i <= num_matched && i <= num_atoms; i++) {
        if (i > 1) fprintf(json_file, ", ");
        if (matched_atoms && matched_atoms[i]) {
            json_escape_string(matched_atoms[i], esc_atom, sizeof(esc_atom));
            fprintf(json_file, "\"%s\"", esc_atom);
        } else {
            fprintf(json_file, "\"\"");
        }
    }
    
    fprintf(json_file, "]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_pair_validation(long base_i, long base_j,
                                         long is_valid, long bp_type_id,
                                         double dir_x, double dir_y, double dir_z,
                                         double *rtn_val, miscPars *misc_pars) {
    if (!json_writer_is_initialized()) return;
    if (!rtn_val || !misc_pars) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"pair_validation\",\n");
    fprintf(json_file, "      \"base_i\": %ld,\n", base_i);
    fprintf(json_file, "      \"base_j\": %ld,\n", base_j);
    fprintf(json_file, "      \"is_valid\": %ld,\n", is_valid);
    fprintf(json_file, "      \"bp_type_id\": %ld,\n", bp_type_id);
    fprintf(json_file, "      \"direction_vectors\": {\n");
    fprintf(json_file, "        \"dir_x\": %.6f,\n", dir_x);
    fprintf(json_file, "        \"dir_y\": %.6f,\n", dir_y);
    fprintf(json_file, "        \"dir_z\": %.6f\n", dir_z);
    fprintf(json_file, "      },\n");
    fprintf(json_file, "      \"calculated_values\": {\n");
    fprintf(json_file, "        \"dorg\": %.6f,\n", rtn_val[1]);
    fprintf(json_file, "        \"d_v\": %.6f,\n", rtn_val[2]);
    fprintf(json_file, "        \"plane_angle\": %.6f,\n", rtn_val[3]);
    fprintf(json_file, "        \"dNN\": %.6f,\n", rtn_val[4]);
    fprintf(json_file, "        \"quality_score\": %.6f\n", rtn_val[5]);
    fprintf(json_file, "      },\n");
    fprintf(json_file, "      \"validation_checks\": {\n");
    fprintf(json_file, "        \"distance_check\": %s,\n",
            (dval_in_range(rtn_val[1], misc_pars->min_dorg, misc_pars->max_dorg)) ? "true" : "false");
    fprintf(json_file, "        \"d_v_check\": %s,\n",
            (dval_in_range(rtn_val[2], misc_pars->min_dv, misc_pars->max_dv)) ? "true" : "false");
    fprintf(json_file, "        \"plane_angle_check\": %s,\n",
            (dval_in_range(rtn_val[3], misc_pars->min_plane_angle, misc_pars->max_plane_angle)) ? "true" : "false");
    fprintf(json_file, "        \"dNN_check\": %s\n",
            (dval_in_range(rtn_val[4], misc_pars->min_dNN, misc_pars->max_dNN)) ? "true" : "false");
    fprintf(json_file, "      },\n");
    fprintf(json_file, "      \"thresholds\": {\n");
    fprintf(json_file, "        \"min_dorg\": %.6f,\n", misc_pars->min_dorg);
    fprintf(json_file, "        \"max_dorg\": %.6f,\n", misc_pars->max_dorg);
    fprintf(json_file, "        \"min_dv\": %.6f,\n", misc_pars->min_dv);
    fprintf(json_file, "        \"max_dv\": %.6f,\n", misc_pars->max_dv);
    fprintf(json_file, "        \"min_plane_angle\": %.6f,\n", misc_pars->min_plane_angle);
    fprintf(json_file, "        \"max_plane_angle\": %.6f,\n", misc_pars->max_plane_angle);
    fprintf(json_file, "        \"min_dNN\": %.6f,\n", misc_pars->min_dNN);
    fprintf(json_file, "        \"max_dNN\": %.6f\n", misc_pars->max_dNN);
    fprintf(json_file, "      }\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_hbond_list(long base_i, long base_j, long num_hbonds,
                                    char **hb_atom1, char **hb_atom2,
                                    double *hb_dist, char *hb_type,
                                    const char *hb_info_string) {
    long i;
    char esc_atom1[BUF512], esc_atom2[BUF512], esc_info[BUF1K];
    
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"hbond_list\",\n");
    fprintf(json_file, "      \"base_i\": %ld,\n", base_i);
    fprintf(json_file, "      \"base_j\": %ld,\n", base_j);
    fprintf(json_file, "      \"num_hbonds\": %ld,\n", num_hbonds);
    json_escape_string(hb_info_string, esc_info, sizeof(esc_info));
    fprintf(json_file, "      \"hb_info_string\": \"%s\",\n", esc_info);
    fprintf(json_file, "      \"hbonds\": [\n");
    
    for (i = 1; i <= num_hbonds && hb_atom1 && hb_atom2; i++) {
        if (i > 1) fprintf(json_file, ",\n");
        json_escape_string(hb_atom1[i], esc_atom1, sizeof(esc_atom1));
        json_escape_string(hb_atom2[i], esc_atom2, sizeof(esc_atom2));
        fprintf(json_file, "        {\n");
        fprintf(json_file, "          \"hbond_idx\": %ld,\n", i);
        fprintf(json_file, "          \"donor_atom\": \"%s\",\n", esc_atom1);
        fprintf(json_file, "          \"acceptor_atom\": \"%s\",\n", esc_atom2);
        if (hb_dist) {
            fprintf(json_file, "          \"distance\": %.6f", fabs(hb_dist[i]));
        } else {
            fprintf(json_file, "          \"distance\": null");
        }
        if (hb_type && hb_type[i]) {
            fprintf(json_file, ",\n          \"type\": \"%c\"", hb_type[i]);
        }
        fprintf(json_file, "\n        }");
    }
    
    fprintf(json_file, "\n      ]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_frame_calc(long residue_idx, char base_type,
                                    const char *template_file, double rms_fit,
                                    long num_matched_atoms,
                                    double **matched_std_xyz,
                                    double **matched_exp_xyz) {
    long i;
    char esc_template[BUF512];
    
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"frame_calc\",\n");
    fprintf(json_file, "      \"residue_idx\": %ld,\n", residue_idx);
    fprintf(json_file, "      \"base_type\": \"%c\",\n", base_type);
    json_escape_string(template_file, esc_template, sizeof(esc_template));
    fprintf(json_file, "      \"template_file\": \"%s\",\n", esc_template);
    fprintf(json_file, "      \"rms_fit\": %.6f,\n", rms_fit);
    fprintf(json_file, "      \"num_matched_atoms\": %ld,\n", num_matched_atoms);
    
    if (matched_std_xyz && matched_exp_xyz) {
        fprintf(json_file, "      \"matched_coordinates\": [\n");
        for (i = 1; i <= num_matched_atoms; i++) {
            if (i > 1) fprintf(json_file, ",\n");
            fprintf(json_file, "        {\n");
            fprintf(json_file, "          \"atom_idx\": %ld,\n", i);
            fprintf(json_file, "          \"std_xyz\": [%.6f, %.6f, %.6f],\n",
                    matched_std_xyz[i][1], matched_std_xyz[i][2], matched_std_xyz[i][3]);
            fprintf(json_file, "          \"exp_xyz\": [%.6f, %.6f, %.6f]\n",
                    matched_exp_xyz[i][1], matched_exp_xyz[i][2], matched_exp_xyz[i][3]);
            fprintf(json_file, "        }");
        }
        fprintf(json_file, "\n      ]\n");
    }
    
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_ring_atoms(long residue_idx, long *ring_atom_indices,
                                    long num_ring_atoms) {
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!ring_atom_indices) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"ring_atoms\",\n");
    fprintf(json_file, "      \"residue_idx\": %ld,\n", residue_idx);
    fprintf(json_file, "      \"num_ring_atoms\": %ld,\n", num_ring_atoms);
    fprintf(json_file, "      \"ring_atom_indices\": [");
    
    for (i = 1; i <= num_ring_atoms; i++) {
        if (i > 1) fprintf(json_file, ", ");
        fprintf(json_file, "%ld", ring_atom_indices[i]);
    }
    
    fprintf(json_file, "]\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_distance_checks(long base_i, long base_j,
                                         double dorg, double dNN,
                                         double plane_angle, double d_v,
                                         double overlap_area) {
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"distance_checks\",\n");
    fprintf(json_file, "      \"base_i\": %ld,\n", base_i);
    fprintf(json_file, "      \"base_j\": %ld,\n", base_j);
    fprintf(json_file, "      \"values\": {\n");
    fprintf(json_file, "        \"dorg\": %.6f,\n", dorg);
    fprintf(json_file, "        \"dNN\": %.6f,\n", dNN);
    fprintf(json_file, "        \"plane_angle\": %.6f,\n", plane_angle);
    fprintf(json_file, "        \"d_v\": %.6f,\n", d_v);
    fprintf(json_file, "        \"overlap_area\": %.6f\n", overlap_area);
    fprintf(json_file, "      }\n");
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_ls_fitting(long residue_idx, long num_points,
                                    double rms_fit, double **rotation_matrix,
                                    double *translation) {
    long i, j;
    
    if (!json_writer_is_initialized()) return;
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"ls_fitting\",\n");
    fprintf(json_file, "      \"residue_idx\": %ld,\n", residue_idx);
    fprintf(json_file, "      \"num_points\": %ld,\n", num_points);
    fprintf(json_file, "      \"rms_fit\": %.6f,\n", rms_fit);
    
    if (rotation_matrix) {
        fprintf(json_file, "      \"rotation_matrix\": ");
        json_write_matrix(json_file, rotation_matrix);
        fprintf(json_file, ",\n");
    }
    
    if (translation) {
        fprintf(json_file, "      \"translation\": ");
        json_write_double_array(json_file, &translation[1], 3);
        fprintf(json_file, "\n");
    }
    
    fprintf(json_file, "    }");
    
    fflush(json_file);
}
