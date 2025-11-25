#include "x3dna.h"
#include "json_writer.h"

// Test tool to isolate donor/acceptor type determination
// Usage: test_donor_acceptor <base1> <base2> <atom1> <atom2> [output.json]
// Example: test_donor_acceptor C G " N3 " " N2 " > test_n3_n2.json

int main(int argc, char *argv[]) {
    char base1, base2;
    char *atom1, *atom2;
    char hb_type;
    FILE *json_fp = stdout;
    
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <base1> <base2> <atom1> <atom2> [output.json]\n", argv[0]);
        fprintf(stderr, "  base1, base2: One-letter base codes (A, C, G, T, U)\n");
        fprintf(stderr, "  atom1, atom2: Atom names with spaces (e.g., \" N3 \", \" N2 \")\n");
        fprintf(stderr, "Example: %s C G \" N3 \" \" N2 \"\n", argv[0]);
        return 1;
    }
    
    base1 = argv[1][0];
    base2 = argv[2][0];
    atom1 = argv[3];
    atom2 = argv[4];
    
    if (argc >= 6) {
        json_fp = open_file(argv[5], "w");
    }
    
    // Initialize globals (needed for some static data)
    set_my_globals(argv[0]);
    
    // Call donor_acceptor function
    hb_type = donor_acceptor(base1, base2, atom1, atom2);
    
    // Output JSON with detailed information
    fprintf(json_fp, "{\n");
    fprintf(json_fp, "  \"base1\": \"%c\",\n", base1);
    fprintf(json_fp, "  \"base2\": \"%c\",\n", base2);
    fprintf(json_fp, "  \"atom1\": \"%s\",\n", atom1);
    fprintf(json_fp, "  \"atom2\": \"%s\",\n", atom2);
    fprintf(json_fp, "  \"hbond_type\": \"%c\",\n", hb_type);
    fprintf(json_fp, "  \"type_description\": \"%s\"\n", 
            (hb_type == '-') ? "standard" : (hb_type == '*') ? "non-standard" : "invalid");
    fprintf(json_fp, "}\n");
    
    // Also print to stderr for quick viewing
    fprintf(stderr, "Base pair: %c-%c\n", base1, base2);
    fprintf(stderr, "Atoms: %s -> %s\n", atom1, atom2);
    fprintf(stderr, "H-bond type: %c (%s)\n", hb_type,
            (hb_type == '-') ? "standard" : (hb_type == '*') ? "non-standard" : "invalid");
    
    // Cleanup
    if (json_fp != stdout) {
        close_file(json_fp);
    }
    
    clear_my_globals();
    
    return 0;
}

