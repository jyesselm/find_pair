#include "json_writer.h"
#include <sys/stat.h>
#include <sys/types.h>

/* Structure to track open type files */
typedef struct {
    char calc_type[BUF32];
    FILE *file;
    long entry_count;
} TypeFileHandle;

/* PDB line cache structure */
typedef struct {
    char **lines;        /* Array of line strings */
    long max_lines;      /* Maximum number of lines allocated */
    long num_lines;      /* Actual number of lines loaded */
    char file_path[BUF1K]; /* Path to cached file */
} PdbLineCache;

/* Singleton state */
static long initialized = FALSE;
static long json_disabled = FALSE; /* Flag to disable JSON output */
static FILE *json_file = NULL;
static TypeFileHandle type_files[32] = {{0}}; /* File handles for each calculation type */
static long num_type_files = 0;
static char json_filename[BUF1K];
static char json_base_name[BUF512]; /* Base name without extension */
static char json_dir_path[BUF1K]; /* Directory path for split files */
static char globals_filename[BUF1K];
static char pdb_file_path[BUF1K] = "";
static long first_entry = TRUE;
static PdbLineCache pdb_line_cache = {NULL, 0, 0, ""}; /* Cache for PDB lines */
static long use_split_files = TRUE; /* Write to separate files per type */

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

void json_writer_disable(void) {
    json_disabled = TRUE;
}

long json_writer_init(const char *pdbfile) {
    char pdb_name[BUF512];
    char dir_path[BUF1K];
    struct stat st = {0};
    
    if (json_disabled) {
        return FALSE; /* JSON output is disabled */
    }
    
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
    
    /* Store base name and directory for split files */
    strncpy(json_base_name, pdb_name, sizeof(json_base_name) - 1);
    json_base_name[sizeof(json_base_name) - 1] = '\0';
    strncpy(json_dir_path, dir_path, sizeof(json_dir_path) - 1);
    json_dir_path[sizeof(json_dir_path) - 1] = '\0';
    
    /* Create JSON filename with corrected path (for metadata file) */
    if (stat("../data", &st) == 0) {
        sprintf(json_filename, "../data/json_legacy/%s.json", pdb_name);
        sprintf(globals_filename, "../data/json_legacy/%s_globals.json", pdb_name);
    } else {
        sprintf(json_filename, "data/json_legacy/%s.json", pdb_name);
        sprintf(globals_filename, "data/json_legacy/%s_globals.json", pdb_name);
    }
    
    /* Store PDB file path for later use */
    strncpy(pdb_file_path, pdbfile, sizeof(pdb_file_path) - 1);
    pdb_file_path[sizeof(pdb_file_path) - 1] = '\0';
    
    /* No longer create main JSON file - only split files are written */
    json_file = NULL;
    first_entry = TRUE;
    initialized = TRUE;
    
    fprintf(stderr, "[JSON_WRITER] Initialized for split files in %s/json_legacy/%s/*.json\n", dir_path, json_base_name);
    return TRUE;
}

/* Helper: Get or create file handle for a calculation type */
static FILE* get_type_file_handle(const char *calc_type, long *is_first_entry) {
    char type_filename[BUF1K];
    FILE *fp;
    struct stat st;
    long i;
    
    fprintf(stderr, "[DEBUG] get_type_file_handle: ENTRY calc_type=%s\n", calc_type ? calc_type : "NULL");
    
    if (!json_writer_is_initialized() || !calc_type) {
        fprintf(stderr, "[DEBUG] get_type_file_handle: Not initialized or calc_type is NULL, returning NULL\n");
        return NULL;
    }
    
    /* Check if file is already open */
    for (i = 0; i < num_type_files; i++) {
        if (strcmp(type_files[i].calc_type, calc_type) == 0 && type_files[i].file != NULL) {
            fprintf(stderr, "[DEBUG] get_type_file_handle: Found existing file handle for %s\n", calc_type);
            *is_first_entry = (type_files[i].entry_count == 0);
            type_files[i].entry_count++;
            return type_files[i].file;
        }
    }
    
    fprintf(stderr, "[DEBUG] get_type_file_handle: Creating new file for %s\n", calc_type);
    
    /* Determine file path - new structure: <record_type>/<PDB_ID>.json */
    char type_dir[BUF1K];
    if (stat("../data", &st) == 0) {
        sprintf(type_dir, "../data/json_legacy/%s", calc_type);
        sprintf(type_filename, "../data/json_legacy/%s/%s.json", calc_type, json_base_name);
    } else {
        sprintf(type_dir, "data/json_legacy/%s", calc_type);
        sprintf(type_filename, "data/json_legacy/%s/%s.json", calc_type, json_base_name);
    }
    
    fprintf(stderr, "[DEBUG] get_type_file_handle: type_dir=%s type_filename=%s\n", type_dir, type_filename);
    
    /* Create record-type directory if it doesn't exist */
    if (stat(type_dir, &st) == -1) {
        fprintf(stderr, "[DEBUG] get_type_file_handle: Creating directory %s\n", type_dir);
        if (mkdir(type_dir, 0755) != 0) {
            fprintf(stderr, "[JSON_WRITER] Warning: Could not create directory %s\n", type_dir);
            return NULL;
        }
    }
    
    /* Create new file */
    fprintf(stderr, "[DEBUG] get_type_file_handle: Opening file %s\n", type_filename);
    fp = fopen(type_filename, "w");
    if (!fp) {
        fprintf(stderr, "[JSON_WRITER] Warning: Could not create %s\n", type_filename);
        return NULL;
    }
    fprintf(stderr, "[DEBUG] get_type_file_handle: File opened successfully, writing array start\n");
    fprintf(fp, "[\n");
    
    /* Store in type_files array */
    if (num_type_files < 32) {
        strncpy(type_files[num_type_files].calc_type, calc_type, sizeof(type_files[num_type_files].calc_type) - 1);
        type_files[num_type_files].calc_type[sizeof(type_files[num_type_files].calc_type) - 1] = '\0';
        type_files[num_type_files].file = fp;
        type_files[num_type_files].entry_count = 1;
        num_type_files++;
        fprintf(stderr, "[DEBUG] get_type_file_handle: Stored in type_files array, num_type_files=%ld\n", num_type_files);
    } else {
        fprintf(stderr, "[DEBUG] get_type_file_handle: WARNING num_type_files >= 32, cannot store more files\n");
    }
    
    *is_first_entry = TRUE;
    fprintf(stderr, "[DEBUG] get_type_file_handle: EXIT successful, returning file handle\n");
    return fp;
}


void json_writer_finalize(void) {
    long i;
    
    fprintf(stderr, "[DEBUG] json_writer_finalize: ENTRY\n");
    
    if (!initialized) {
        fprintf(stderr, "[DEBUG] json_writer_finalize: Not initialized, returning\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_finalize: num_type_files=%ld\n", num_type_files);
    
    /* Close all open type files */
    for (i = 0; i < num_type_files; i++) {
        fprintf(stderr, "[DEBUG] json_writer_finalize: Processing type file %ld: calc_type=%s file=%p\n", 
                i, type_files[i].calc_type, (void*)type_files[i].file);
        if (type_files[i].file != NULL) {
            fprintf(stderr, "[DEBUG] json_writer_finalize: Closing file for %s\n", type_files[i].calc_type);
            fprintf(type_files[i].file, "\n]");
            fprintf(stderr, "[DEBUG] json_writer_finalize: Flushing file for %s\n", type_files[i].calc_type);
            fflush(type_files[i].file);
            fprintf(stderr, "[DEBUG] json_writer_finalize: Closing file handle for %s\n", type_files[i].calc_type);
            fclose(type_files[i].file);
            type_files[i].file = NULL;
            fprintf(stderr, "[DEBUG] json_writer_finalize: Closed file for %s\n", type_files[i].calc_type);
        }
    }
    
    /* No main JSON file to close - only split files are used */
    if (json_file) {
        fprintf(stderr, "[DEBUG] json_writer_finalize: Closing main json_file\n");
        fclose(json_file);
        json_file = NULL;
    }
    initialized = FALSE;
    num_type_files = 0;
    
    fprintf(stderr, "[DEBUG] json_writer_finalize: About to print final message\n");
    fprintf(stderr, "[JSON_WRITER] Finalized: split files written to %s/json_legacy/%s/*.json\n", json_dir_path, json_base_name);
    fprintf(stderr, "[DEBUG] json_writer_finalize: EXIT successful\n");
}

long json_writer_is_initialized(void) {
    return initialized && !json_disabled; /* No longer require json_file to be open */
}

void json_writer_record_bpstep_params(long bp_idx1, long bp_idx2, 
                                      double *pars,
                                      double *mst_org,
                                      double **mst_orien) {
    if (!json_writer_is_initialized()) return;
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    FILE *type_file;
    long is_first;
    char esc_bp_type[BUF32];
    
    if (!json_writer_is_initialized()) return;
    
    /* Write to separate file for easier comparison (original identification from find_pair) */
    type_file = get_type_file_handle("base_pair", &is_first);
    if (!type_file) return;
    
    if (!is_first) fprintf(type_file, ",\n");
    
    json_escape_string(bp_type, esc_bp_type, sizeof(esc_bp_type));
    
    fprintf(type_file, "    {\n");
    fprintf(type_file, "      \"type\": \"base_pair\",\n");
    fprintf(type_file, "      \"base_i\": %ld,\n", i);
    fprintf(type_file, "      \"base_j\": %ld,\n", j);
    if (bp_type) {
        fprintf(type_file, "      \"bp_type\": \"%s\",\n", esc_bp_type);
    }
    if (dir_xyz) {
        fprintf(type_file, "      \"dir_xyz\": [%.6f, %.6f, %.6f],\n",
                dir_xyz[1], dir_xyz[2], dir_xyz[3]);
    }
    fprintf(type_file, "      \"orien_i\": ");
    json_write_matrix(type_file, orien_i);
    fprintf(type_file, ",\n");
    fprintf(type_file, "      \"orien_j\": ");
    json_write_matrix(type_file, orien_j);
    fprintf(type_file, ",\n");
    fprintf(type_file, "      \"org_i\": ");
    json_write_double_array(type_file, &org_i[1], 3);
    fprintf(type_file, ",\n");
    fprintf(type_file, "      \"org_j\": ");
    json_write_double_array(type_file, &org_j[1], 3);
    fprintf(type_file, "\n");
    fprintf(type_file, "    }");
    
    fflush(type_file);
}

void json_writer_record_ref_frame(long residue_idx, double **orien, double *org) {
    if (!json_writer_is_initialized()) return;
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    
    fprintf(stderr, "[DEBUG] json_writer_record_sequence: ENTRY num_residue=%ld\n", num_residue);
    fflush(stderr);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_sequence: Not initialized, returning\n");
        return;
    }
    if (!bseq) {
        fprintf(stderr, "[DEBUG] json_writer_record_sequence: bseq is NULL, returning\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_sequence: json_file=%p\n", (void*)json_file);
    fflush(stderr);
    
    if (!json_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_sequence: json_file is NULL, returning (using split files)\n");
        return;
    }
    
    json_escape_string(bseq, esc_seq, sizeof(esc_seq));
    
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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

/* Helper: Load PDB lines into cache */
static void load_pdb_lines_cache(const char *pdbfile) {
    FILE *fp;
    char *line = NULL;
    long i;
    
    /* Check if already cached for this file */
    if (pdb_line_cache.lines != NULL && strcmp(pdb_line_cache.file_path, pdbfile) == 0) {
        return; /* Already cached */
    }
    
    /* Free old cache if different file */
    if (pdb_line_cache.lines != NULL) {
        for (i = 1; i <= pdb_line_cache.num_lines; i++) {
            if (pdb_line_cache.lines[i]) {
                free(pdb_line_cache.lines[i]);
            }
        }
        free(pdb_line_cache.lines);
        pdb_line_cache.lines = NULL;
    }
    
    /* Initialize cache */
    pdb_line_cache.max_lines = 100000; /* Initial size */
    pdb_line_cache.num_lines = 0;
    pdb_line_cache.lines = (char **)calloc(pdb_line_cache.max_lines + 1, sizeof(char *));
    if (!pdb_line_cache.lines) {
        fprintf(stderr, "[JSON_WRITER] Warning: Could not allocate PDB line cache\n");
        return;
    }
    strncpy(pdb_line_cache.file_path, pdbfile, sizeof(pdb_line_cache.file_path) - 1);
    pdb_line_cache.file_path[sizeof(pdb_line_cache.file_path) - 1] = '\0';
    
    /* Load all lines from file */
    fp = fopen(pdbfile, "r");
    if (!fp) {
        fprintf(stderr, "[JSON_WRITER] Warning: Could not open PDB file for caching: %s\n", pdbfile);
        return;
    }
    
    while ((line = my_getline(fp)) != NULL) {
        pdb_line_cache.num_lines++;
        if (pdb_line_cache.num_lines > pdb_line_cache.max_lines) {
            /* Reallocate if needed */
            long new_max = pdb_line_cache.max_lines * 2;
            char **new_lines = (char **)realloc(pdb_line_cache.lines, (new_max + 1) * sizeof(char *));
            if (new_lines) {
                pdb_line_cache.lines = new_lines;
                pdb_line_cache.max_lines = new_max;
                /* Initialize new entries */
                for (i = pdb_line_cache.num_lines + 1; i <= new_max; i++) {
                    pdb_line_cache.lines[i] = NULL;
                }
            }
        }
        pdb_line_cache.lines[pdb_line_cache.num_lines] = line;
    }
    
    fclose(fp);
    fprintf(stderr, "[JSON_WRITER] Cached %ld PDB lines from %s\n", pdb_line_cache.num_lines, pdbfile);
}

/* Helper: Get PDB line from cache by line number */
static char* get_pdb_line_by_number(const char *pdbfile, long line_num, char *buffer, size_t buf_size) {
    size_t len;
    
    if (!pdbfile || line_num <= 0 || !buffer || buf_size == 0) {
        buffer[0] = '\0';
        return buffer;
    }
    
    /* Load cache if not already loaded or different file */
    if (pdb_line_cache.lines == NULL || strcmp(pdb_line_cache.file_path, pdbfile) != 0) {
        load_pdb_lines_cache(pdbfile);
    }
    
    /* Check if line number is valid */
    if (line_num > pdb_line_cache.num_lines || !pdb_line_cache.lines[line_num]) {
        buffer[0] = '\0';
        return buffer;
    }
    
    /* Copy line from cache */
    strncpy(buffer, pdb_line_cache.lines[line_num], buf_size - 1);
    buffer[buf_size - 1] = '\0';
    
    /* Remove trailing newline if present */
    len = strlen(buffer);
    if (len > 0 && buffer[len - 1] == '\n') {
        buffer[len - 1] = '\0';
    }
    
    return buffer;
}

void json_writer_record_pdb_atoms(long num, char **AtomName, char **ResName,
                                   char *ChainID, long *ResSeq, double **xyz,
                                   char **Miscs, long *line_numbers) {
    char esc_atom[BUF32], esc_res[BUF32], esc_pdb_line[BUF512];
    char pdb_line_buffer[BUF512];
    FILE *type_file;
    long i;
    long is_first = TRUE;
    
    if (!json_writer_is_initialized()) return;
    if (!AtomName || !ResName || !ChainID || !ResSeq || !xyz) return;
    
    /* Get file handle for pdb_atoms type */
    type_file = get_type_file_handle("pdb_atoms", &is_first);
    if (!type_file) return;
    
    /* Write entry to type-specific file */
    if (!is_first) {
        fprintf(type_file, ",\n");
    }
    fprintf(type_file, "  {\n");
    fprintf(type_file, "    \"num_atoms\": %ld,\n", num);
    fprintf(type_file, "    \"atoms\": [\n");
    
    for (i = 1; i <= num; i++) {
        if (i > 1) fprintf(type_file, ",\n");
        fprintf(type_file, "      {\n");
        
        /* Atom index/position in array */
        fprintf(type_file, "        \"atom_idx\": %ld,\n", i);
        
        /* Line number in PDB file */
        if (line_numbers && line_numbers[i] > 0) {
            fprintf(type_file, "        \"line_number\": %ld,\n", line_numbers[i]);
            
            /* Get and include the original PDB line for debugging */
            if (pdb_file_path[0] != '\0') {
                get_pdb_line_by_number(pdb_file_path, line_numbers[i], pdb_line_buffer, sizeof(pdb_line_buffer));
                if (pdb_line_buffer[0] != '\0') {
                    json_escape_string(pdb_line_buffer, esc_pdb_line, sizeof(esc_pdb_line));
                    fprintf(type_file, "        \"pdb_line\": \"%s\",\n", esc_pdb_line);
                }
            }
        }
        
        /* Atom name */
        json_escape_string(AtomName[i], esc_atom, sizeof(esc_atom));
        fprintf(type_file, "        \"atom_name\": \"%s\",\n", esc_atom);
        
        /* Residue name */
        json_escape_string(ResName[i], esc_res, sizeof(esc_res));
        fprintf(type_file, "        \"residue_name\": \"%s\",\n", esc_res);
        
        /* Chain ID */
        fprintf(type_file, "        \"chain_id\": \"%c\",\n", ChainID[i]);
        
        /* Residue sequence number */
        fprintf(type_file, "        \"residue_seq\": %ld,\n", ResSeq[i]);
        
        /* Coordinates */
        fprintf(type_file, "        \"xyz\": [%.6f, %.6f, %.6f]", 
                xyz[i][1], xyz[i][2], xyz[i][3]);
        
        /* Miscellaneous info if available */
        if (Miscs && Miscs[i]) {
            fprintf(type_file, ",\n        \"record_type\": \"%c\"", Miscs[i][0]);
            if (Miscs[i][1] != ' ') {
                fprintf(type_file, ",\n        \"alt_loc\": \"%c\"", Miscs[i][1]);
            }
            if (Miscs[i][2] != ' ') {
                fprintf(type_file, ",\n        \"insertion\": \"%c\"", Miscs[i][2]);
            }
        }
        
        fprintf(type_file, "\n      }");
    }
    
    fprintf(type_file, "\n    ]\n");
    fprintf(type_file, "  }");
    
    /* Don't close file - keep it open for more entries of this type */
    fflush(type_file);
}

void json_writer_record_residue_indices(long num_residue, long **seidx) {
    long i;
    FILE *type_file;
    long is_first;
    
    fprintf(stderr, "[DEBUG] json_writer_record_residue_indices: ENTRY num_residue=%ld\n", num_residue);
    fflush(stderr);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_residue_indices: Not initialized, returning\n");
        return;
    }
    if (!seidx) {
        fprintf(stderr, "[DEBUG] json_writer_record_residue_indices: seidx is NULL, returning\n");
        return;
    }
    
    /* Write to separate file for easier comparison (matches other record functions) */
    type_file = get_type_file_handle("residue_indices", &is_first);
    fprintf(stderr, "[DEBUG] json_writer_record_residue_indices: get_type_file_handle returned type_file=%p is_first=%ld\n", 
            (void*)type_file, is_first);
    if (!type_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_residue_indices: type_file is NULL, returning\n");
        return;
    }
    
    if (!is_first) fprintf(type_file, ",\n");
    
    fprintf(type_file, "  {\n");
    fprintf(type_file, "    \"type\": \"residue_indices\",\n");
    fprintf(type_file, "    \"num_residue\": %ld,\n", num_residue);
    fprintf(type_file, "    \"seidx\": [\n");
    
    for (i = 1; i <= num_residue; i++) {
        if (i > 1) fprintf(type_file, ",\n");
        fprintf(type_file, "      {\"residue_idx\": %ld, \"start_atom\": %ld, \"end_atom\": %ld}",
                i, seidx[i][1], seidx[i][2]);
    }
    
    fprintf(type_file, "\n    ]\n");
    fprintf(type_file, "  }");
    
    fflush(type_file);
    fprintf(stderr, "[DEBUG] json_writer_record_residue_indices: EXIT successful\n");
}

void json_writer_record_base_pairs(long ds, long num_bp, long **pair_num) {
    long i, j;
    
    fprintf(stderr, "[DEBUG] json_writer_record_base_pairs: ENTRY ds=%ld num_bp=%ld\n", ds, num_bp);
    fflush(stderr);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_base_pairs: Not initialized, returning\n");
        return;
    }
    if (!pair_num) {
        fprintf(stderr, "[DEBUG] json_writer_record_base_pairs: pair_num is NULL, returning\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_base_pairs: json_file=%p\n", (void*)json_file);
    fflush(stderr);
    
    if (!json_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_base_pairs: json_file is NULL, returning (using split files)\n");
        return;
    }
    
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
    
    fprintf(stderr, "[DEBUG] json_writer_record_all_ref_frames: ENTRY ds=%ld num_bp=%ld\n", ds, num_bp);
    fflush(stderr);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_all_ref_frames: Not initialized, returning\n");
        return;
    }
    if (!orien || !org) {
        fprintf(stderr, "[DEBUG] json_writer_record_all_ref_frames: NULL pointers, returning\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_all_ref_frames: json_file=%p\n", (void*)json_file);
    fflush(stderr);
    
    if (!json_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_all_ref_frames: json_file is NULL, returning (using split files)\n");
        return;
    }
    
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
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    
    fprintf(stderr, "[DEBUG] json_writer_record_helices: ENTRY num_bp=%ld\n", num_bp);
    fflush(stderr);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_helices: Not initialized, returning\n");
        return;
    }
    if (!bphlx) {
        fprintf(stderr, "[DEBUG] json_writer_record_helices: bphlx is NULL, returning\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_helices: json_file=%p\n", (void*)json_file);
    fflush(stderr);
    
    if (!json_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_helices: json_file is NULL, returning (using split files)\n");
        return;
    }
    
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
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    if (!json_file) return; /* Using split files - json_file is NULL */
    
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
    /* NOTE: This function is deprecated - we now use json_writer_record_hbond_list
     * which writes to split files. This function is kept for backward compatibility
     * but disabled when using split files (json_file is NULL).
     */
    if (!json_writer_is_initialized()) return;
    
    /* Skip if using split files (json_file is NULL) - hbond_list already recorded */
    /* This prevents crashes from trying to write to NULL json_file */
    return;
}

void json_writer_record_base_frame_calc(long residue_idx, char base_type,
                                         const char *standard_template,
                                         double rms_fit, long num_matched,
                                         char **matched_atoms, long num_atoms,
                                         const char *residue_name, char chain_id, long residue_seq, char insertion_code) {
    long i;
    char esc_atom[BUF512], esc_template[BUF512], esc_resname[BUF512];
    FILE *type_file;
    long is_first = TRUE;
    
    if (!json_writer_is_initialized()) return;
    
    type_file = get_type_file_handle("base_frame_calc", &is_first);
    if (!type_file) return;
    
    if (!is_first) fprintf(type_file, ",\n");
    fprintf(type_file, "  {\n");
    fprintf(type_file, "    \"residue_idx\": %ld,\n", residue_idx);
    fprintf(type_file, "    \"base_type\": \"%c\",\n", base_type);
    
    /* Add residue identification information */
    if (residue_name) {
        json_escape_string(residue_name, esc_resname, sizeof(esc_resname));
        fprintf(type_file, "    \"residue_name\": \"%s\",\n", esc_resname);
    }
    fprintf(type_file, "    \"chain_id\": \"%c\",\n", chain_id);
    fprintf(type_file, "    \"residue_seq\": %ld,\n", residue_seq);
    if (insertion_code != ' ') {
        fprintf(type_file, "    \"insertion\": \"%c\",\n", insertion_code);
    }
    
    json_escape_string(standard_template, esc_template, sizeof(esc_template));
    fprintf(type_file, "    \"standard_template\": \"%s\",\n", esc_template);
    fprintf(type_file, "    \"rms_fit\": %.6f,\n", rms_fit);
    fprintf(type_file, "    \"num_matched_atoms\": %ld,\n", num_matched);
    fprintf(type_file, "    \"matched_atoms\": [");
    
    for (i = 1; i <= num_matched && i <= num_atoms; i++) {
        if (i > 1) fprintf(type_file, ", ");
        if (matched_atoms && matched_atoms[i]) {
            json_escape_string(matched_atoms[i], esc_atom, sizeof(esc_atom));
            fprintf(type_file, "\"%s\"", esc_atom);
        } else {
            fprintf(type_file, "\"\"");
        }
    }
    
    fprintf(type_file, "]\n");
    fprintf(type_file, "  }");
    
    fflush(type_file);
}

void json_writer_record_pair_validation(long base_i, long base_j,
                                         long is_valid, long bp_type_id,
                                         double dir_x, double dir_y, double dir_z,
                                         double *rtn_val, miscPars *misc_pars) {
    static long call_count = 0;
    call_count++;
    FILE *type_file;
    long is_first;
    
    if (call_count % 50000 == 0) {
        fprintf(stderr, "[DEBUG] json_writer_record_pair_validation: Call #%ld for pair (%ld, %ld)\n", call_count, base_i, base_j);
        fflush(stderr);
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_pair_validation: ENTRY pair (%ld, %ld) call_count=%ld\n", base_i, base_j, call_count);
    fflush(stderr);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_pair_validation: Not initialized\n");
        return;
    }
    if (!rtn_val || !misc_pars) {
        fprintf(stderr, "[DEBUG] json_writer_record_pair_validation: NULL pointers\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_pair_validation: About to get file handle\n");
    fflush(stderr);
    /* Write to separate file for easier comparison */
    type_file = get_type_file_handle("pair_validation", &is_first);
    fprintf(stderr, "[DEBUG] json_writer_record_pair_validation: Got file handle=%p is_first=%ld\n", (void*)type_file, is_first);
    fflush(stderr);
    if (!type_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_pair_validation: File handle is NULL, returning\n");
        return;
    }
    
    if (!is_first) fprintf(type_file, ",\n");
    
    fprintf(type_file, "    {\n");
    fprintf(type_file, "      \"type\": \"pair_validation\",\n");
    fprintf(type_file, "      \"base_i\": %ld,\n", base_i);
    fprintf(type_file, "      \"base_j\": %ld,\n", base_j);
    fprintf(type_file, "      \"is_valid\": %ld,\n", is_valid);
    fprintf(type_file, "      \"bp_type_id\": %ld,\n", bp_type_id);
    fprintf(type_file, "      \"direction_vectors\": {\n");
    fprintf(type_file, "        \"dir_x\": %.6f,\n", dir_x);
    fprintf(type_file, "        \"dir_y\": %.6f,\n", dir_y);
    fprintf(type_file, "        \"dir_z\": %.6f\n", dir_z);
    fprintf(type_file, "      },\n");
    fprintf(type_file, "      \"calculated_values\": {\n");
    fprintf(type_file, "        \"dorg\": %.6f,\n", rtn_val[1]);
    fprintf(type_file, "        \"d_v\": %.6f,\n", rtn_val[2]);
    fprintf(type_file, "        \"plane_angle\": %.6f,\n", rtn_val[3]);
    fprintf(type_file, "        \"dNN\": %.6f,\n", rtn_val[4]);
    fprintf(type_file, "        \"quality_score\": %.6f\n", rtn_val[5]);
    fprintf(type_file, "      },\n");
    fprintf(type_file, "      \"validation_checks\": {\n");
    fprintf(type_file, "        \"distance_check\": %s,\n",
            (dval_in_range(rtn_val[1], misc_pars->min_dorg, misc_pars->max_dorg)) ? "true" : "false");
    fprintf(type_file, "        \"d_v_check\": %s,\n",
            (dval_in_range(rtn_val[2], misc_pars->min_dv, misc_pars->max_dv)) ? "true" : "false");
    fprintf(type_file, "        \"plane_angle_check\": %s,\n",
            (dval_in_range(rtn_val[3], misc_pars->min_plane_angle, misc_pars->max_plane_angle)) ? "true" : "false");
    fprintf(type_file, "        \"dNN_check\": %s\n",
            (dval_in_range(rtn_val[4], misc_pars->min_dNN, misc_pars->max_dNN)) ? "true" : "false");
    fprintf(type_file, "      },\n");
    fprintf(type_file, "      \"thresholds\": {\n");
    fprintf(type_file, "        \"min_dorg\": %.6f,\n", misc_pars->min_dorg);
    fprintf(type_file, "        \"max_dorg\": %.6f,\n", misc_pars->max_dorg);
    fprintf(type_file, "        \"min_dv\": %.6f,\n", misc_pars->min_dv);
    fprintf(type_file, "        \"max_dv\": %.6f,\n", misc_pars->max_dv);
    fprintf(type_file, "        \"min_plane_angle\": %.6f,\n", misc_pars->min_plane_angle);
    fprintf(type_file, "        \"max_plane_angle\": %.6f,\n", misc_pars->max_plane_angle);
    fprintf(type_file, "        \"min_dNN\": %.6f,\n", misc_pars->min_dNN);
    fprintf(type_file, "        \"max_dNN\": %.6f\n", misc_pars->max_dNN);
    fprintf(type_file, "      }\n");
    fprintf(type_file, "    }");
    
    fflush(type_file);
}

void json_writer_record_hbond_list(long base_i, long base_j, long num_hbonds,
                                    char **hb_atom1, char **hb_atom2,
                                    double *hb_dist, char *hb_type,
                                    const char *hb_info_string) {
    FILE *type_file;
    long is_first;
    long i;
    char esc_atom1[BUF512], esc_atom2[BUF512], esc_info[BUF1K];
    
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: ENTRY base_i=%ld base_j=%ld num_hbonds=%ld\n", base_i, base_j, num_hbonds);
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: hb_atom1=%p hb_atom2=%p hb_dist=%p hb_type=%p hb_info_string=%p\n", 
            (void*)hb_atom1, (void*)hb_atom2, (void*)hb_dist, (void*)hb_type, (void*)hb_info_string);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: NOT INITIALIZED, returning\n");
        return;
    }
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: Initialized check passed\n");
    
    /* Write to separate file for easier comparison (matches other record functions) */
    type_file = get_type_file_handle("hbond_list", &is_first);
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: get_type_file_handle returned type_file=%p is_first=%ld\n", 
            (void*)type_file, is_first);
    if (!type_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: type_file is NULL, returning\n");
        return;
    }
    
    if (!is_first) fprintf(type_file, ",\n");
    
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: About to write record header\n");
    fprintf(type_file, "    {\n");
    fprintf(type_file, "      \"type\": \"hbond_list\",\n");
    fprintf(type_file, "      \"base_i\": %ld,\n", base_i);
    fprintf(type_file, "      \"base_j\": %ld,\n", base_j);
    fprintf(type_file, "      \"num_hbonds\": %ld,\n", num_hbonds);
    
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: About to escape hb_info_string\n");
    if (hb_info_string) {
        json_escape_string(hb_info_string, esc_info, sizeof(esc_info));
        fprintf(type_file, "      \"hb_info_string\": \"%s\",\n", esc_info);
    } else {
        fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: WARNING hb_info_string is NULL\n");
        fprintf(type_file, "      \"hb_info_string\": \"\",\n");
    }
    
    fprintf(type_file, "      \"hbonds\": [\n");
    
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: About to loop through %ld hbonds\n", num_hbonds);
    for (i = 1; i <= num_hbonds && hb_atom1 && hb_atom2; i++) {
        fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: Processing hbond %ld of %ld\n", i, num_hbonds);
        if (i > 1) fprintf(type_file, ",\n");
        
        if (hb_atom1[i]) {
            json_escape_string(hb_atom1[i], esc_atom1, sizeof(esc_atom1));
        } else {
            fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: WARNING hb_atom1[%ld] is NULL\n", i);
            strcpy(esc_atom1, "");
        }
        
        if (hb_atom2[i]) {
            json_escape_string(hb_atom2[i], esc_atom2, sizeof(esc_atom2));
        } else {
            fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: WARNING hb_atom2[%ld] is NULL\n", i);
            strcpy(esc_atom2, "");
        }
        
        fprintf(type_file, "        {\n");
        fprintf(type_file, "          \"hbond_idx\": %ld,\n", i);
        fprintf(type_file, "          \"donor_atom\": \"%s\",\n", esc_atom1);
        fprintf(type_file, "          \"acceptor_atom\": \"%s\",\n", esc_atom2);
        if (hb_dist) {
            fprintf(type_file, "          \"distance\": %.6f", fabs(hb_dist[i]));
        } else {
            fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: WARNING hb_dist is NULL\n");
            fprintf(type_file, "          \"distance\": null");
        }
        if (hb_type && hb_type[i]) {
            fprintf(type_file, ",\n          \"type\": \"%c\"", hb_type[i]);
        }
        fprintf(type_file, "\n        }");
        fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: Finished hbond %ld\n", i);
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: Finished loop, closing record\n");
    fprintf(type_file, "\n      ]\n");
    fprintf(type_file, "    }");
    
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: About to fflush\n");
    fflush(type_file);
    fprintf(stderr, "[DEBUG] json_writer_record_hbond_list: EXIT successful\n");
}

void json_writer_record_frame_calc(long residue_idx, char base_type,
                                    const char *template_file, double rms_fit,
                                    long num_matched_atoms,
                                    double **matched_std_xyz,
                                    double **matched_exp_xyz,
                                    const char *residue_name, char chain_id, long residue_seq, char insertion_code) {
    static long call_count = 0;
    call_count++;
    long i;
    char esc_template[BUF512], esc_resname[BUF512];
    FILE *type_file;
    long is_first = TRUE;
    
    if (call_count % 1000 == 0) {
        fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: Call #%ld residue_idx=%ld\n", call_count, residue_idx);
        fflush(stderr);
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: ENTRY residue_idx=%ld call_count=%ld\n", residue_idx, call_count);
    fflush(stderr);
    
    if (!json_writer_is_initialized()) {
        fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: Not initialized, returning\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: About to get file handle\n");
    fflush(stderr);
    type_file = get_type_file_handle("frame_calc", &is_first);
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: Got file handle=%p is_first=%ld\n", (void*)type_file, is_first);
    fflush(stderr);
    if (!type_file) {
        fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: File handle is NULL, returning\n");
        return;
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: About to write record\n");
    fflush(stderr);
    
    if (!is_first) fprintf(type_file, ",\n");
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: Wrote comma if needed\n");
    fflush(stderr);
    
    fprintf(type_file, "  {\n");
    fprintf(type_file, "    \"residue_idx\": %ld,\n", residue_idx);
    fprintf(type_file, "    \"base_type\": \"%c\",\n", base_type);
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: Wrote header, about to write residue info\n");
    fflush(stderr);
    
    /* Add residue identification information */
    if (residue_name) {
        json_escape_string(residue_name, esc_resname, sizeof(esc_resname));
        fprintf(type_file, "    \"residue_name\": \"%s\",\n", esc_resname);
    }
    fprintf(type_file, "    \"chain_id\": \"%c\",\n", chain_id);
    fprintf(type_file, "    \"residue_seq\": %ld,\n", residue_seq);
    if (insertion_code != ' ') {
        fprintf(type_file, "    \"insertion\": \"%c\",\n", insertion_code);
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: About to escape template file\n");
    fflush(stderr);
    json_escape_string(template_file, esc_template, sizeof(esc_template));
    fprintf(type_file, "    \"template_file\": \"%s\",\n", esc_template);
    fprintf(type_file, "    \"rms_fit\": %.6f,\n", rms_fit);
    fprintf(type_file, "    \"num_matched_atoms\": %ld,\n", num_matched_atoms);
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: About to write matched coordinates\n");
    fflush(stderr);
    if (matched_std_xyz && matched_exp_xyz) {
        fprintf(type_file, "    \"matched_coordinates\": [\n");
        for (i = 1; i <= num_matched_atoms; i++) {
            if (i > 1) fprintf(type_file, ",\n");
            fprintf(type_file, "      {\n");
            fprintf(type_file, "        \"atom_idx\": %ld,\n", i);
            fprintf(type_file, "        \"std_xyz\": [%.6f, %.6f, %.6f],\n",
                    matched_std_xyz[i][1], matched_std_xyz[i][2], matched_std_xyz[i][3]);
            fprintf(type_file, "        \"exp_xyz\": [%.6f, %.6f, %.6f]\n",
                    matched_exp_xyz[i][1], matched_exp_xyz[i][2], matched_exp_xyz[i][3]);
            fprintf(type_file, "      }");
        }
        fprintf(type_file, "\n    ]\n");
    }
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: About to close record\n");
    fflush(stderr);
    fprintf(type_file, "  }");
    
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: About to fflush\n");
    fflush(stderr);
    fflush(type_file);
    fprintf(stderr, "[DEBUG] json_writer_record_frame_calc: EXIT successful\n");
    fflush(stderr);
}

void json_writer_record_ring_atoms(long residue_idx, long *ring_atom_indices,
                                    long num_ring_atoms) {
    long i;
    FILE *type_file;
    long is_first;
    
    if (!json_writer_is_initialized()) return;
    if (!ring_atom_indices) return;
    
    /* Write to separate file for easier comparison */
    type_file = get_type_file_handle("ring_atoms", &is_first);
    if (!type_file) return;
    
    if (!is_first) fprintf(type_file, ",\n");
    
    fprintf(type_file, "    {\n");
    fprintf(type_file, "      \"type\": \"ring_atoms\",\n");
    fprintf(type_file, "      \"residue_idx\": %ld,\n", residue_idx);
    fprintf(type_file, "      \"num_ring_atoms\": %ld,\n", num_ring_atoms);
    fprintf(type_file, "      \"ring_atom_indices\": [");
    
    for (i = 1; i <= num_ring_atoms; i++) {
        if (i > 1) fprintf(type_file, ", ");
        fprintf(type_file, "%ld", ring_atom_indices[i]);
    }
    
    fprintf(type_file, "]\n");
    fprintf(type_file, "    }");
    
    fflush(type_file);
}

void json_writer_record_distance_checks(long base_i, long base_j,
                                         double dorg, double dNN,
                                         double plane_angle, double d_v,
                                         double overlap_area) {
    FILE *type_file;
    long is_first;
    
    if (!json_writer_is_initialized()) return;
    
    /* Write to separate file for easier comparison */
    type_file = get_type_file_handle("distance_checks", &is_first);
    if (!type_file) return;
    
    if (!is_first) fprintf(type_file, ",\n");
    
    fprintf(type_file, "    {\n");
    fprintf(type_file, "      \"type\": \"distance_checks\",\n");
    fprintf(type_file, "      \"base_i\": %ld,\n", base_i);
    fprintf(type_file, "      \"base_j\": %ld,\n", base_j);
    fprintf(type_file, "      \"values\": {\n");
    fprintf(type_file, "        \"dorg\": %.6f,\n", dorg);
    fprintf(type_file, "        \"dNN\": %.6f,\n", dNN);
    fprintf(type_file, "        \"plane_angle\": %.6f,\n", plane_angle);
    fprintf(type_file, "        \"d_v\": %.6f,\n", d_v);
    fprintf(type_file, "        \"overlap_area\": %.6f\n", overlap_area);
    fprintf(type_file, "      }\n");
    fprintf(type_file, "    }");
    
    fflush(type_file);
}

void json_writer_record_ls_fitting(long residue_idx, long num_points,
                                    double rms_fit, double **rotation_matrix,
                                    double *translation,
                                    const char *residue_name, char chain_id, long residue_seq, char insertion_code) {
    long i, j;
    char esc_resname[BUF512];
    FILE *type_file;
    long is_first = TRUE;
    
    if (!json_writer_is_initialized()) return;
    
    type_file = get_type_file_handle("ls_fitting", &is_first);
    if (!type_file) return;
    
    if (!is_first) fprintf(type_file, ",\n");
    fprintf(type_file, "  {\n");
    fprintf(type_file, "    \"residue_idx\": %ld,\n", residue_idx);
    
    /* Add residue identification information */
    if (residue_name) {
        json_escape_string(residue_name, esc_resname, sizeof(esc_resname));
        fprintf(type_file, "    \"residue_name\": \"%s\",\n", esc_resname);
    }
    fprintf(type_file, "    \"chain_id\": \"%c\",\n", chain_id);
    fprintf(type_file, "    \"residue_seq\": %ld,\n", residue_seq);
    if (insertion_code != ' ') {
        fprintf(type_file, "    \"insertion\": \"%c\",\n", insertion_code);
    }
    
    fprintf(type_file, "    \"num_points\": %ld,\n", num_points);
    fprintf(type_file, "    \"rms_fit\": %.6f,\n", rms_fit);
    
    if (rotation_matrix) {
        fprintf(type_file, "    \"rotation_matrix\": ");
        json_write_matrix(type_file, rotation_matrix);
        fprintf(type_file, ",\n");
    }
    
    if (translation) {
        fprintf(type_file, "    \"translation\": ");
        json_write_double_array(type_file, &translation[1], 3);
        fprintf(type_file, "\n");
    }
    
    fprintf(type_file, "  }");
    
    fflush(type_file);
}

void json_writer_record_removed_atom(const char* pdb_line, const char* reason,
                                     long atom_serial, const char* atom_name,
                                     const char* residue_name, char chain_id,
                                     long residue_seq, double* xyz, long model_num) {
    char esc_line[BUF1K], esc_reason[BUF512], esc_atom[BUF32], esc_res[BUF32];
    long has_fields = FALSE;
    
    if (!json_writer_is_initialized()) return;
    if (!json_file) return; /* Using split files - json_file is NULL */
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    json_escape_string(pdb_line ? pdb_line : "", esc_line, sizeof(esc_line));
    json_escape_string(reason ? reason : "unknown", esc_reason, sizeof(esc_reason));
    json_escape_string(atom_name ? atom_name : "", esc_atom, sizeof(esc_atom));
    json_escape_string(residue_name ? residue_name : "", esc_res, sizeof(esc_res));
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"removed_atom\",\n");
    fprintf(json_file, "      \"reason\": \"%s\"", esc_reason);
    has_fields = TRUE;
    
    if (pdb_line && strlen(pdb_line) > 0) {
        fprintf(json_file, ",\n      \"pdb_line\": \"%s\"", esc_line);
        has_fields = TRUE;
    }
    
    if (atom_serial > 0) {
        fprintf(json_file, ",\n      \"atom_serial\": %ld", atom_serial);
    }
    
    if (atom_name && strlen(atom_name) > 0) {
        fprintf(json_file, ",\n      \"atom_name\": \"%s\"", esc_atom);
    }
    
    if (residue_name && strlen(residue_name) > 0) {
        fprintf(json_file, ",\n      \"residue_name\": \"%s\"", esc_res);
    }
    
    if (chain_id != ' ') {
        fprintf(json_file, ",\n      \"chain_id\": \"%c\"", chain_id);
    }
    
    if (residue_seq > 0) {
        fprintf(json_file, ",\n      \"residue_seq\": %ld", residue_seq);
    }
    
    if (xyz) {
        fprintf(json_file, ",\n      \"xyz\": [%.6f, %.6f, %.6f]",
                xyz[0], xyz[1], xyz[2]);
    }
    
    if (model_num >= 0) {
        fprintf(json_file, ",\n      \"model_num\": %ld", model_num);
    }
    
    fprintf(json_file, "\n    }");
    
    fflush(json_file);
}

void json_writer_record_removed_atoms_summary(long num_removed) {
    if (!json_writer_is_initialized()) return;
    if (!json_file) return; /* Using split files - json_file is NULL */
    
    if (!first_entry) fprintf(json_file, ",\n");
    first_entry = FALSE;
    
    fprintf(json_file, "    {\n");
    fprintf(json_file, "      \"type\": \"removed_atoms_summary\",\n");
    fprintf(json_file, "      \"num_removed\": %ld\n", num_removed);
    fprintf(json_file, "    }");
    
    fflush(json_file);
}

void json_writer_record_find_bestpair_selection(long num_bp, long **base_pairs) {
    FILE *type_file;
    long is_first;
    long i;
    
    if (!json_writer_is_initialized()) return;
    if (!base_pairs) return;
    
    /* Write to separate file for easier comparison (original selection from find_bestpair) */
    type_file = get_type_file_handle("find_bestpair_selection", &is_first);
    if (!type_file) return;
    
    if (!is_first) fprintf(type_file, ",\n");
    
    fprintf(type_file, "    {\n");
    fprintf(type_file, "      \"type\": \"find_bestpair_selection\",\n");
    fprintf(type_file, "      \"num_bp\": %ld,\n", num_bp);
    fprintf(type_file, "      \"pairs\": [\n");
    
    for (i = 1; i <= num_bp; i++) {
        if (i > 1) fprintf(type_file, ",\n");
        fprintf(type_file, "        [%ld, %ld]", base_pairs[i][1], base_pairs[i][2]);
    }
    
    fprintf(type_file, "\n      ]\n");
    fprintf(type_file, "    }");
    
    fflush(type_file);
}
