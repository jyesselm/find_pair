#include "x3dna.h"
#include "json_writer.h"

// Define struct_args - it's defined differently in find_pair.c and analyze.c
// We'll use a union or define both versions
// For find_pair.c:
typedef struct {
    char pdbfile[BUF512];
    char outfile[BUF512];
    char map[BUF512];
    long ds;
    long curves;
    long curves_plus;
    long divide;
    long hetatm;
    long pairs;
    long detailed;
    long waters;
    long hjb;
} struct_args_fp;

// For analyze.c:
typedef struct {
    char torsion[BUF512];
    long istart;
    long istep;
    long icnt;
    long waters;
    long bz;
    long ring;
    long simple_pars;
    long abi;
    long circular;
} struct_args_ana;

// Forward declarations - these functions are now non-static
// They're defined in find_pair.c and analyze.c which are included in the build
extern void handle_str(struct_args_fp *args);  // from find_pair.c
extern void process_str(char *inpfile, struct_args_ana *args);  // from analyze.c
extern void fp_cmdline(int argc, char *argv[], struct_args_fp *args);  // from find_pair.c
extern void analyze_cmdline(int argc, char *argv[], struct_args_ana *args);  // from analyze.c

static void combined_usage(void)
{
    fprintf(stderr, "Usage: find_pair_analyze [options] <pdb_file> [outfile]\n\n");
    fprintf(stderr, "This program runs both find_pair and analyze in sequence.\n\n");
    fprintf(stderr, "Options from both find_pair and analyze are supported.\n");
    fprintf(stderr, "See 'find_pair -h' and 'analyze -h' for detailed options.\n\n");
    fprintf(stderr, "Example: find_pair_analyze 1EHZ.pdb\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    struct_args_fp fp_args;
    struct_args_ana ana_args;
    time_t time0;
    char inpfile[BUF512], parfile[BUF512];
    char **fp_argv;
    int fp_argc;
    int i, j;
    
    time(&time0);
    set_my_globals(argv[0]);
    
    if (argc < 2)
        combined_usage();
    
    // Determine PDB file name (last non-option argument)
    for (i = argc - 1; i >= 1; i--) {
        if (*argv[i] != '-')
            break;
    }
    if (i < 1)
        combined_usage();
    
    // Build argument array for find_pair
    fp_argc = i + 1;
    fp_argv = (char**)malloc((fp_argc + 2) * sizeof(char*));
    fp_argv[0] = "find_pair";
    for (j = 1; j <= i; j++)
        fp_argv[j] = argv[j];
    if (i + 1 < argc)
        fp_argv[fp_argc++] = argv[i + 1];  // optional outfile
    fp_argv[fp_argc] = NULL;
    
    // Parse find_pair arguments
    memset(&fp_args, 0, sizeof(fp_args));
    fp_cmdline(fp_argc, fp_argv, &fp_args);
    
    // Determine input file name for analyze
    del_extension(fp_args.pdbfile, parfile);
    if (is_equal_string(fp_args.outfile, "stdout")) {
        sprintf(inpfile, "%s.inp", parfile);
        strcpy(fp_args.outfile, inpfile);
    } else {
        strcpy(inpfile, fp_args.outfile);
    }
    
    // Initialize analyze args with defaults
    memset(&ana_args, 0, sizeof(ana_args));
    strcpy(ana_args.torsion, "");
    ana_args.istart = 1;
    ana_args.istep = 1;
    ana_args.icnt = FALSE;
    ana_args.waters = fp_args.waters;
    ana_args.bz = TRUE;
    ana_args.ring = FALSE;
    ana_args.simple_pars = TRUE;
    ana_args.abi = FALSE;
    ana_args.circular = FALSE;
    
    // Parse analyze-specific options from command line
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            continue;
        if (check_global_options(argv[i]))
            continue;
        if (lux_ncmatch(argv[i], "^--?t")) {
            get_strvalue(argv[i], ana_args.torsion, FALSE);
            continue;
        }
        if (lux_ncmatch(argv[i], "^--?bz")) {
            ana_args.bz = set_switch_default_true(argv[i]);
            continue;
        }
        if (lux_ncmatch(argv[i], "^--?ri")) {
            ana_args.ring = set_switch_default_true(argv[i]);
            continue;
        }
        if (lux_ncmatch(argv[i], "^--?si")) {
            long k_val = FALSE;
            if (!lux_ncmatch(argv[i], "no|false|off")) {
                k_val = TRUE;
                if (lux_ncmatch(argv[i], "n1|n9"))
                    k_val |= 2;  // SIMPLE_BP_LONG_AXIS_RN9_YN1
                if (lux_ncmatch(argv[i], "heli"))
                    k_val |= 4;  // SIMPLE_STEP_HELICAL_PARS
            }
            ana_args.simple_pars = k_val;
            continue;
        }
        if (lux_ncmatch(argv[i], "^--?abi?")) {
            ana_args.abi = set_switch_default_true(argv[i]);
            continue;
        }
        if (lux_ncmatch(argv[i], "^--?circ")) {
            ana_args.circular = set_switch_default_true(argv[i]);
            continue;
        }
        if (strchr(argv[i], 'C') && !strchr(argv[i], '=') && !strstr(argv[i], "curves")) {
            ana_args.icnt = TRUE;
        }
        if (strchr(argv[i], 'W')) {
            ana_args.waters = TRUE;
        }
        if (strstr(argv[i], "-S=")) {
            long istart, istep;
            if (sscanf(argv[i], "-S=%ld,%ld", &istep, &istart) == 2) {
                if (istart < 0)
                    istart = -istart;
                ana_args.istart = istart;
                ana_args.istep = istep;
                fprintf(stderr, "***start at %ld, with step size: %ld***\n", istart, istep);
            }
        }
    }
    
    fprintf(stderr, "\n=== Step 1: Running find_pair on <%s> ===\n", fp_args.pdbfile);
    
    // Initialize JSON writer BEFORE running find_pair to suppress file creation
    json_writer_init(fp_args.pdbfile);
    json_writer_record_global_variables();
    
    // Step 1: Run find_pair
    handle_str(&fp_args);
    
    // Check if .inp file was created
    if (!exist_file(inpfile)) {
        fprintf(stderr, "\nError: find_pair did not create input file <%s>\n", inpfile);
        fprintf(stderr, "Cannot proceed with analyze step.\n");
        free(fp_argv);
        clear_my_globals();
        return 1;
    }
    
    fprintf(stderr, "\n=== Step 2: Running analyze on <%s> ===\n", inpfile);
    
    // Step 2: Run analyze
    // Clean up analyze output files first (as analyze main does)
    remove_file(AUX_FILE);
    remove_file(BPSTEP_FILE);
    remove_file(HLXSTEP_FILE);
    remove_file(STACK_FILE);
    remove_file(HSTACK_FILE);
    remove_file(REF_FILE);
    remove_file(POC_FILE);
    remove_file(SEVEN_FILE);
    
    process_str(inpfile, &ana_args);
    
    fprintf(stderr, "\n=== Combined analysis complete ===\n");
    
    // Finalize JSON writer
    json_writer_finalize();
    
    free(fp_argv);
    clear_my_globals();
    print_used_time(time0);
    return 0;
}
