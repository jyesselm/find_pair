#ifndef JSON_WRITER_H
#define JSON_WRITER_H

#include "x3dna.h"

/* JSON Writer Singleton for logging key calculations */

/* Initialize the JSON writer with PDB filename */
/* Extracts PDB name from full path and creates data/json_legacy/{pdb_name}.json
 */
long json_writer_init(const char* pdbfile);

/* Finalize and write JSON file */
void json_writer_finalize(void);

/* Record base pair step parameters */
void json_writer_record_bpstep_params(
    long bp_idx1, long bp_idx2, double* pars, /* 6 params: Shift, Slide, Rise, Tilt, Roll, Twist */
    double* mst_org,                          /* 3 coords */
    double** mst_orien /* 3x3 matrix */);

/* Record helical parameters */
void json_writer_record_helical_params(long bp_idx1, long bp_idx2, double* pars, /* 6 params */
                                       double* mst_orgH,                         /* 3 coords */
                                       double** mst_orienH /* 3x3 matrix */);

/* Record base pair information */
void json_writer_record_base_pair(long i, long j, char* bp_type,
                                  double* dir_xyz,  /* dir_x, dir_y, dir_z */
                                  double** orien_i, /* 3x3 matrix for base i */
                                  double** orien_j, /* 3x3 matrix for base j */
                                  double* org_i,    /* 3 coords for base i */
                                  double* org_j);   /* 3 coords for base j */

/* Record reference frame calculation */
void json_writer_record_ref_frame(long residue_idx, double** orien, double* org);

/* Record sequence information */
void json_writer_record_sequence(long num_residue, char* bseq);

/* Record base pair sequence */
void json_writer_record_bp_sequence(long num_bp, char** bp_seq, long ds);

/* Record PDB atom data */
void json_writer_record_pdb_atoms(long num, char** AtomName, char** ResName, char* ChainID,
                                  long* ResSeq, double** xyz, char** Miscs, long* line_numbers);

/* Record residue index mapping (seidx) */
void json_writer_record_residue_indices(long num_residue, long** seidx);

/* Record base pair information */
void json_writer_record_base_pairs(long ds, long num_bp, long** pair_num);

/* Record reference frames for all residues */
void json_writer_record_all_ref_frames(long ds, long num_bp, double** orien, double** org);

/* Record WC_info (Watson-Crick information) */
void json_writer_record_wc_info(long num_bp, long* WC_info);

/* Record RY (Purine/Pyrimidine classification) */
void json_writer_record_ry(long num_residue, long* RY);

/* Record helix information (bphlx) */
void json_writer_record_helices(long num_bp, long* bphlx);

/* Record twist and rise for each step */
void json_writer_record_twist_rise(long nbpm1, double** twist_rise);

/* Record input parameters and thresholds */
void json_writer_record_input_parameters(miscPars* misc_pars, long ds, long hetatm, long ip);

/* Record all global variables (writes to separate file for constants) */
void json_writer_record_global_variables(void);

/* Record hydrogen bond details for a base pair */
void json_writer_record_hbonds(long base_i, long base_j, long num_hbonds, char** hb_atom1,
                               char** hb_atom2, double* hb_dist, char* hb_type, long* lkg_type);

/* Record base frame calculation details (from base_frame or ref_frames) */
void json_writer_record_base_frame_calc(long residue_idx, char base_type,
                                        const char* standard_template, double rms_fit,
                                        long num_matched, char** matched_atoms, long num_atoms,
                                        const char* residue_name, char chain_id, long residue_seq,
                                        char insertion_code);

/* Record pair validation results (from check_pair) */
void json_writer_record_pair_validation(long base_i, long base_j, long is_valid, long bp_type_id,
                                        double dir_x, double dir_y, double dir_z,
                                        double* rtn_val, /* rtn_val[1..5] */
                                        miscPars* misc_pars);

/* Record hydrogen bond list (detailed from get_hbond_ij) */
void json_writer_record_hbond_list(long base_i, long base_j, long num_hbonds, char** hb_atom1,
                                   char** hb_atom2, double* hb_dist, char* hb_type,
                                   const char* hb_info_string);

/* Record base frame calculation with full details (alternative signature) */
void json_writer_record_frame_calc(long residue_idx, char base_type, const char* template_file,
                                   double rms_fit, long num_matched_atoms, double** matched_std_xyz,
                                   double** matched_exp_xyz, const char* residue_name,
                                   char chain_id, long residue_seq, char insertion_code);

/* Record ring atom indices for a residue */
void json_writer_record_ring_atoms(long residue_idx, long* ring_atom_indices, long num_ring_atoms);

/* Record distance and angle checks */
void json_writer_record_distance_checks(long base_i, long base_j, double dorg, double dNN,
                                        double plane_angle, double d_v, double overlap_area);

/* Record least squares fitting details */
void json_writer_record_ls_fitting(long residue_idx, long num_points, double rms_fit,
                                   double** rotation_matrix, double* translation,
                                   const char* residue_name, char chain_id, long residue_seq,
                                   char insertion_code);

/* Record removed atoms during PDB parsing */
void json_writer_record_removed_atom(const char* pdb_line, const char* reason, long atom_serial,
                                     const char* atom_name, const char* residue_name, char chain_id,
                                     long residue_seq, double* xyz, long model_num);

/* Record all removed atoms (call after parsing) */
void json_writer_record_removed_atoms_summary(long num_removed);

/* Record original base pair selection from find_bestpair (before reordering) */
void json_writer_record_find_bestpair_selection(long num_bp, long **base_pairs);

/* Check if JSON writer is initialized */
long json_writer_is_initialized(void);

#endif /* JSON_WRITER_H */
