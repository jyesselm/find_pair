#ifndef JSON_WRITER_H
#define JSON_WRITER_H

#include "x3dna.h"

/* JSON Writer Singleton for logging key calculations */

/* Initialize the JSON writer with PDB filename */
/* Extracts PDB name from full path and creates data/json_legacy/{pdb_name}.json
 */
long json_writer_init(const char *pdbfile);

/* Finalize and write JSON file */
void json_writer_finalize(void);

/* Record base pair step parameters */
void json_writer_record_bpstep_params(
    long bp_idx1, long bp_idx2,
    double *pars,    /* 6 params: Shift, Slide, Rise, Tilt, Roll, Twist */
    double *mst_org, /* 3 coords */
    double **mst_orien /* 3x3 matrix */);

/* Record helical parameters */
void json_writer_record_helical_params(long bp_idx1, long bp_idx2,
                                       double *pars,     /* 6 params */
                                       double *mst_orgH, /* 3 coords */
                                       double **mst_orienH /* 3x3 matrix */);

/* Record base pair information */
void json_writer_record_base_pair(long i, long j, char *bp_type,
                                  double *dir_xyz,  /* dir_x, dir_y, dir_z */
                                  double **orien_i, /* 3x3 matrix for base i */
                                  double **orien_j, /* 3x3 matrix for base j */
                                  double *org_i,    /* 3 coords for base i */
                                  double *org_j);   /* 3 coords for base j */

/* Record reference frame calculation */
void json_writer_record_ref_frame(long residue_idx, double **orien,
                                  double *org);

/* Record sequence information */
void json_writer_record_sequence(long num_residue, char *bseq);

/* Record base pair sequence */
void json_writer_record_bp_sequence(long num_bp, char **bp_seq, long ds);

/* Check if JSON writer is initialized */
long json_writer_is_initialized(void);

#endif /* JSON_WRITER_H */
