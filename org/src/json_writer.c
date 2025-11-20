#include "json_writer.h"
#include <sys/stat.h>
#include <sys/types.h>

/* Singleton state */
static long initialized = FALSE;
static FILE *json_file = NULL;
static char json_filename[BUF1K];
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
    char data_dir[BUF1K];
    struct stat st = {0};
    
    if (initialized) {
        return TRUE; /* Already initialized */
    }
    
    /* Extract base name without extension */
    bname_noext((char *)pdbfile, pdb_name);
    
    /* Determine project root data directory */
    /* Check if ../data exists (we're in org/ directory) */
    if (stat("../data", &st) == 0) {
        sprintf(data_dir, "../data");
        sprintf(dir_path, "../data/json_legacy");
    } else if (stat("data", &st) == 0) {
        /* We're already at project root */
        sprintf(data_dir, "data");
        sprintf(dir_path, "data/json_legacy");
    } else {
        /* Create data directory at current location */
        sprintf(data_dir, "data");
        sprintf(dir_path, "data/json_legacy");
    }
    
    /* Create parent directory 'data' if it doesn't exist */
    if (stat(data_dir, &st) == -1) {
        if (mkdir(data_dir, 0755) != 0) {
            fprintf(stderr, "[JSON_WRITER] Warning: Could not create directory %s\n", data_dir);
            return FALSE;
        }
    }
    /* Create json_legacy directory */
    if (stat(dir_path, &st) == -1) {
        if (mkdir(dir_path, 0755) != 0) {
            fprintf(stderr, "[JSON_WRITER] Warning: Could not create directory %s\n", dir_path);
            return FALSE;
        }
    }
    
    /* Create JSON filename */
    sprintf(json_filename, "%s/%s.json", dir_path, pdb_name);
    
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

