# Legacy Code Workflows

**Date**: 2025-01-XX  
**Purpose**: Complete workflow diagrams showing execution order, function calls, and data flow  
**Status**: Comprehensive workflow reference

---

## Table of Contents

1. [Complete find_pair Workflow](#complete-find_pair-workflow)
2. [Complete analyze Workflow](#complete-analyze-workflow)
3. [Base Frame Calculation Workflow](#base-frame-calculation-workflow)
4. [Base Pair Finding Workflow](#base-pair-finding-workflow)
5. [H-Bond Detection Workflow](#h-bond-detection-workflow)
6. [Function Dependency Map](#function-dependency-map)

---

## Complete find_pair Workflow

### High-Level Flow

```
find_pair_main(argc, argv)
├─ 1. Initialization
│  ├─ set_my_globals("find_pair")
│  ├─ fp_cmdline(argc, argv, &args)  # Parse command line
│  ├─ get_3dna_pars(&misc_pars)       # Load parameters
│  └─ json_writer_init(pdbfile)      # Initialize JSON (if enabled)
│
├─ 2. Main Processing
│  └─ handle_str(&args)
│     ├─ 2.1. PDB File Reading
│     │   ├─ number_of_atoms(pdbfile, hetatm, alt_list)
│     │   │   └─ Returns: num_atoms
│     │   ├─ Allocate arrays:
│     │   │   ├─ AtomName = cmatrix(1, num_atoms, 0, 4)
│     │   │   ├─ ResName = cmatrix(1, num_atoms, 0, 3)
│     │   │   ├─ ChainID = cvector(1, num_atoms)
│     │   │   ├─ ResSeq = lvector(1, num_atoms)
│     │   │   ├─ xyz = dmatrix(1, num_atoms, 1, 3)
│     │   │   └─ Miscs = cmatrix(1, num_atoms, 0, NMISC)
│     │   ├─ read_pdb(pdbfile, ..., AtomName, ResName, ..., xyz, Miscs, ...)
│     │   │   └─ Parses ATOM/HETATM records, fills arrays
│     │   └─ json_writer_record_pdb_atoms(...)  # Record atom data
│     │
│     ├─ 2.2. Residue Indexing
│     │   ├─ atom_idx(num_atoms, AtomName, NULL, idx)
│     │   │   └─ Classifies atoms: N=1, O=2, etc.
│     │   ├─ residue_idx(num_atoms, ResSeq, Miscs, ChainID, ResName, &num_residue)
│     │   │   ├─ Creates unique residue identifiers
│     │   │   ├─ Finds residue boundaries
│     │   │   └─ Returns: seidx[1..num_residue][1..2]
│     │   └─ json_writer_record_residue_indices(...)  # Record residue mapping
│     │
│     ├─ 2.3. Sequence Extraction
│     │   ├─ get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY)
│     │   │   ├─ For each residue i:
│     │   │   │   ├─ residue_ident(...) → RY[i] (1=purine, 0=pyrimidine, -1=not nucleotide)
│     │   │   │   └─ base_ident(...) → bseq[i] ('A', 'C', 'G', 'T', etc.)
│     │   │   └─ Returns: bseq[1..num_residue], RY[1..num_residue]
│     │   └─ json_writer_record_ry(...)  # Record RY classification
│     │
│     ├─ 2.4. Base Frame Calculation
│     │   ├─ base_info(num_residue, bseq, seidx, RY, ..., xyz, orien, org, NC1xyz, o3_p)
│     │   │   ├─ base_frame(num_residue, bseq, seidx, RY, ..., xyz, BDIR, orien, org)
│     │   │   │   ├─ For each residue i in [1, num_residue]:
│     │   │   │   │   ├─ if RY[i] < 0: continue  # Skip non-nucleotides
│     │   │   │   │   ├─ Load template: Atomic_{bseq[i]}.pdb
│     │   │   │   │   ├─ Match ring atoms (RA_LIST)
│     │   │   │   │   ├─ if nmatch < 3: skip
│     │   │   │   │   ├─ ls_fitting(...) → R, orgi, RMS
│     │   │   │   │   ├─ mst2orien(orien[i], 0, R)  # Store rotation
│     │   │   │   │   └─ org[i] = orgi  # Store origin
│     │   │   │   └─ Returns: orien[][], org[][]
│     │   │   ├─ Extract N and C1' coordinates → NC1xyz[][]
│     │   │   └─ Extract O3' and P coordinates → o3_p[][]
│     │   └─ json_writer_record_base_frame_calc(...)  # Record frame calculations
│     │
│     ├─ 2.5. Ring Atom Indexing
│     │   ├─ ring_oidx(num_atoms, num_residue, RY, seidx, AtomName, xyz, idx, ring_atom)
│     │   │   ├─ get_bonds(...) → connect[][]
│     │   │   ├─ all_bring_atoms(...) → ring_atom[][]
│     │   │   └─ get_cntatom(...) → connected atoms
│     │   └─ json_writer_record_ring_atoms(...)  # Record ring atoms
│     │
│     ├─ 2.6. Base Pair Finding
│     │   ├─ find_bestpair(nout, base_pairs, num_residue, bseq, seidx, RY, AtomName, xyz, idx, orien, org, NC1xyz, ring_atom, misc_pars)
│     │   │   ├─ Initialize: matched_idx[] = zeros, num_bp = 0
│     │   │   ├─ while matched_count increases:
│     │   │   │   ├─ For each residue i in [1, num_residue]:
│     │   │   │   │   ├─ if RY[i] < 0 or matched_idx[i]: continue
│     │   │   │   │   ├─ best_pair(i, ..., pair_istat)
│     │   │   │   │   │   ├─ For each j != i:
│     │   │   │   │   │   │   ├─ if RY[j] < 0 or matched_idx[j]: continue
│     │   │   │   │   │   │   ├─ check_pair(i, j, ..., rtn_val, &bpid, ..., 0)
│     │   │   │   │   │   │   │   ├─ get_bp_zoave(i, j, orien, org, oave, zave)
│     │   │   │   │   │   │   │   ├─ Calculate: dorg, dv, plane_angle, dNN
│     │   │   │   │   │   │   │   ├─ Check constraints (cdns)
│     │   │   │   │   │   │   │   ├─ Check overlap: get_oarea(...)
│     │   │   │   │   │   │   │   ├─ Count H-bonds: get_hbond_ij(...)
│     │   │   │   │   │   │   │   ├─ if valid: calculate_more_bppars(...) → bpid
│     │   │   │   │   │   │   │   └─ Calculate quality: rtn_val[5]
│     │   │   │   │   │   │   └─ if bpid != 0 and rtn_val[5] < best: update best
│     │   │   │   │   │   └─ Return: j, bpid, quality scores
│     │   │   │   │   ├─ if pair_istat[1] != 0:  # Found partner j
│     │   │   │   │   │   ├─ best_pair(j, ..., pair_jstat)  # Check mutual
│     │   │   │   │   │   └─ if i == pair_jstat[1]:  # Mutual match!
│     │   │   │   │   │     ├─ matched_idx[i] = 1
│     │   │   │   │   │     ├─ matched_idx[j] = 1
│     │   │   │   │   │     └─ base_pairs[++num_bp] = (i, j, bpid, ...)
│     │   │   │   └─ Check if matched_count increased
│     │   │   └─ Returns: num_bp, base_pairs[][]
│     │   └─ json_writer_record_find_bestpair_selection(...)  # Record selected pairs
│     │
│     ├─ 2.7. Base Pair Reordering
│     │   ├─ re_ordering(num_bp, base_pairs, bp_idx, helix_marker, helix_idx, misc_pars, &num_helix, o3_p, ...)
│     │   │   ├─ bp_context(...)  # Analyze neighbor relationships
│     │   │   ├─ locate_helix(...)  # Group pairs into helical regions
│     │   │   ├─ five2three(...)  # Reorder to 5'→3' direction
│     │   │   └─ check_zdna(...)  # Special Z-DNA handling
│     │   └─ Returns: bp_idx[], helix_idx[][], helix_marker[], num_helix
│     │
│     └─ 2.8. Output Generation
│        ├─ x3dna_input(outfile, ...)  # Write .inp file
│        └─ write_bestpairs(...)      # Write pair information
│
└─ 3. Cleanup
   ├─ Free all allocated arrays
   ├─ clear_my_globals()
   ├─ json_writer_finalize()  # Close JSON files
   └─ print_used_time(time0)
```

---

## Complete analyze Workflow

### High-Level Flow

```
analyze_main(argc, argv)
├─ 1. Initialization
│  ├─ set_my_globals("analyze")
│  ├─ analyze_cmdline(argc, argv, &args)
│  └─ json_writer_init(pdbfile)  # If not already initialized
│
├─ 2. Main Processing
│  └─ process_str(inpfile, &args)
│     ├─ 2.1. Read Input File
│     │   ├─ read_input(inpfile, pdbfile, outfile, &ds, &num_bp, &ip, &hetatm)
│     │   │   ├─ Parse .inp file from find_pair
│     │   │   ├─ Extract: PDB file name, duplex status, pair numbers
│     │   │   └─ Load pair_num[1..ds][1..num_bp]: Residue indices for each base pair
│     │   └─ json_writer_record_base_pairs(...)  # Record pair numbers
│     │
│     ├─ 2.2. Read PDB File (again)
│     │   ├─ number_of_atoms(pdbfile, hetatm, alt_list)
│     │   ├─ Allocate atom arrays (same as find_pair)
│     │   ├─ read_pdb(pdbfile, ..., AtomName, ..., xyz, Miscs, ...)
│     │   └─ json_writer_record_pdb_atoms(...)  # Record atom data
│     │
│     ├─ 2.3. Residue Indexing (again)
│     │   ├─ atom_idx(num_atoms, AtomName, NULL, idx)
│     │   ├─ residue_idx(num_atoms, ResSeq, Miscs, ChainID, ResName, &num_residue)
│     │   └─ json_writer_record_residue_indices(...)
│     │
│     ├─ 2.4. Base Pair Validation
│     │   ├─ pair_checking(ip, ds, num_residue, pdbfile, &num_bp, pair_num)
│     │   │   └─ Validate pair indices are in valid range
│     │   └─ Result: Validated pair_num[][]
│     │
│     ├─ 2.5. Extract Base Sequences
│     │   ├─ get_bpseq(ds, num_bp, pair_num, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bp_seq, RY)
│     │   │   ├─ For each base pair k:
│     │   │   │   ├─ i = pair_num[1][k], j = pair_num[2][k]
│     │   │   │   ├─ residue_ident(...) → RY[i], RY[j]
│     │   │   │   └─ base_ident(...) → bp_seq[1][k], bp_seq[2][k]
│     │   │   └─ Returns: bp_seq[1..ds][1..num_bp], RY[]
│     │   ├─ get_seq(num_residue, seidx, ..., bseq, RY)  # Full sequence
│     │   └─ json_writer_record_ry(...)
│     │
│     ├─ 2.6. Recalculate Base Frames
│     │   ├─ ref_frames(ds, num_bp, pair_num, bp_seq, seidx, RY, ..., xyz, fp, orien, org, WC_info, ...)
│     │   │   ├─ For each base pair k:
│     │   │   │   ├─ i = pair_num[1][k], j = pair_num[2][k]
│     │   │   │   ├─ base_frame(...) → Calculate frame for base i
│     │   │   │   ├─ base_frame(...) → Calculate frame for base j
│     │   │   │   ├─ Check orientation: dir_z = dot(z_i, z_j)
│     │   │   │   ├─ If needed: reverse one frame
│     │   │   │   └─ Store: orien[1..ds][offset], org[1..ds][offset]
│     │   │   ├─ check_Watson_Crick(num_bp, bp_seq, WC_info)
│     │   │   │   └─ Classify pairs: WC, wobble, other
│     │   │   └─ Returns: orien[][], org[][], WC_info[]
│     │   └─ json_writer_record_all_ref_frames(...)  # Record frames
│     │
│     └─ 2.7. Calculate Step Parameters
│        ├─ get_parameters(ds, num_bp, bp_seq, orien, org, WC_info, fp, ...)
│        │   ├─ For each consecutive pair (k, k+1):
│        │   │   ├─ refs_i_j(k, k+1, orien, org, r1, o1, r2, o2)
│        │   │   │   └─ Extract frames from orien[][], org[][]
│        │   │   ├─ bpstep_par(r1, o1, r2, o2, pars, mst_orien, mst_org)
│        │   │   │   ├─ Calculate: Shift, Slide, Rise, Tilt, Roll, Twist
│        │   │   │   └─ Returns: pars[1..6], mst_orien[], mst_org[]
│        │   │   ├─ helical_par(r1, o1, r2, o2, pars, ...)
│        │   │   │   └─ Calculate helical parameters
│        │   │   └─ Store: bp_step_par[][], bp_heli_par[][]
│        │   └─ Returns: Step parameters for all pairs
│        └─ json_writer_record_bpstep_params(...)  # Record step parameters
│
└─ 3. Cleanup
   ├─ Free all arrays
   ├─ clear_my_globals()
   └─ print_used_time(time0)
```

---

## Base Frame Calculation Workflow

### Detailed Flow

```
base_frame(num_residue, bseq, seidx, res_type, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, BDIR, orien, org)

For each residue i in [1, num_residue]:
│
├─ Step 1: Check if nucleotide
│  └─ if res_type[i] < 0: continue  # Skip non-nucleotides
│
├─ Step 2: Get atom range
│  ├─ ib = seidx[i][1]  # Start atom index
│  └─ ie = seidx[i][2]  # End atom index
│
├─ Step 3: Load standard template
│  ├─ set_std_base_pdb(BDIR, FALSE, bseq[i], spdb)
│  │   └─ Construct: spdb = "{BDIR}/Atomic_{bseq[i]}.pdb"
│  ├─ read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, sMiscs, 1, "*")
│  │   └─ Load template PDB file
│  └─ snum = number of atoms in template
│
├─ Step 4: Match ring atoms
│  ├─ RingAtom[] = RA_LIST = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "}
│  ├─ RingAtom_num = (res_type[i] == 1) ? 9 : 6  # Purine=9, Pyrimidine=6
│  ├─ nmatch = 0
│  ├─ For j = 0 to RingAtom_num-1:
│  │   ├─ exp_katom = find_1st_atom(RingAtom[j], AtomName, ib, ie, idmsg)
│  │   ├─ std_katom = find_1st_atom(RingAtom[j], sAtomName, 1, snum, sidmsg)
│  │   └─ if exp_katom && std_katom:
│  │      ├─ cpxyz(xyz[exp_katom], eRing_xyz[++nmatch])
│  │      └─ cpxyz(sxyz[std_katom], sRing_xyz[nmatch])
│  └─ if nmatch < 3: continue  # Need minimum 3 atoms
│
├─ Step 5: Least-squares fitting
│  ├─ ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi)
│  │   ├─ Calculate covariance matrix U
│  │   ├─ Build quaternion matrix N
│  │   ├─ Eigenvalue decomposition: jacobi(N, 4, D, V)
│  │   ├─ Extract quaternion: q = V[1..4][4]
│  │   ├─ Convert to rotation matrix R
│  │   ├─ Calculate translation: orgi = ave_exyz - R × ave_sxyz
│  │   └─ Calculate RMS
│  └─ Returns: RMS, R[3][3], orgi[3]
│
└─ Step 6: Store frame
   ├─ mst2orien(orien[i], 0, R)
   │   └─ Convert R[1..3][1..3] to orien[i][1..9]
   └─ cpxyz(orgi, org[i])
      └─ Store orgi[1..3] in org[i][1..3]
```

**JSON Recording**: `json_writer_record_base_frame_calc(i, bseq[i], template_file, rms, orien[i], org[i])`

---

## Base Pair Finding Workflow

### Detailed Flow

```
find_bestpair(nout, base_pairs, num_residue, bseq, seidx, RY, AtomName, xyz, idx, orien, org, NC1xyz, ring_atom, misc_pars)

Initialize:
├─ matched_idx[1..num_residue] = all zeros
├─ num_bp = 0
└─ num1 = 0, num2 = 1

while (num1 < num2):  # Continue until no new pairs found
│
├─ num1 = num2  # Save current count
│
├─ For each residue i in [1, num_residue]:
│  │
│  ├─ Skip if:
│  │  ├─ RY[i] < 0  # Not a nucleotide
│  │  └─ matched_idx[i] == 1  # Already matched
│  │
│  ├─ Find best partner for i:
│  │  └─ best_pair(i, num_residue, RY, seidx, xyz, idx, NC1xyz, matched_idx, orien, org, ring_atom, AtomName, bseq, misc_pars, pair_istat)
│  │     │
│  │     ├─ Initialize: best_score = XBIG, pair_istat[1] = 0
│  │     │
│  │     ├─ For each residue j in [1, num_residue]:
│  │     │  ├─ Skip if: j == i, RY[j] < 0, matched_idx[j] == 1
│  │     │  │
│  │     │  ├─ Validate pair:
│  │     │  │  └─ check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx, AtomName, ring_atom, misc_pars, rtn_val, &bpid, dir_x, dir_y, dir_z, 0)
│  │     │  │     │
│  │     │  │     ├─ Calculate geometry:
│  │     │  │     │  ├─ get_bp_zoave(i, j, orien, org, oave, zave)
│  │     │  │     │  ├─ dorg = org[j] - org[i]
│  │     │  │     │  ├─ dNN_vec = NC1xyz[j] - NC1xyz[i]
│  │     │  │     │  └─ Calculate: dir_x, dir_y, dir_z
│  │     │  │     │
│  │     │  │     ├─ Calculate metrics:
│  │     │  │     │  ├─ rtn_val[1] = veclen(dorg)  # dorg
│  │     │  │     │  ├─ rtn_val[2] = fabs(dot(dorg, zave))  # dv
│  │     │  │     │  ├─ rtn_val[3] = vec_ang(z_i, z_j, NULL)  # plane_angle
│  │     │  │     │  └─ rtn_val[4] = veclen(dNN_vec)  # dNN
│  │     │  │     │
│  │     │  │     ├─ Check constraints (cdns):
│  │     │  │     │  ├─ cdns = (rtn_val[1] in [min_dorg, max_dorg]) &&
│  │     │  │     │  │        (rtn_val[2] in [min_dv, max_dv]) &&
│  │     │  │     │  │        (rtn_val[3] <= max_plane_angle) &&
│  │     │  │     │  │        (rtn_val[4] in [min_dNN, max_dNN])
│  │     │  │     │  └─ if !cdns: return (bpid = 0)
│  │     │  │     │
│  │     │  │     ├─ Check overlap:
│  │     │  │     │  ├─ overlap = get_oarea(i, j, ring_atom, oave, zave, xyz, 0)
│  │     │  │     │  └─ if overlap >= OVERLAP: return (bpid = 0)
│  │     │  │     │
│  │     │  │     ├─ Count H-bonds:
│  │     │  │     │  ├─ num_base_hb = 0
│  │     │  │     │  ├─ For m in [seidx[i][1], seidx[i][2]]:
│  │     │  │     │  │   ├─ For n in [seidx[j][1], seidx[j][2]]:
│  │     │  │     │  │   │   ├─ if good_hbatoms(...) && within_limits(...):
│  │     │  │     │  │   │   │   └─ num_base_hb++
│  │     │  │     │  └─ if num_base_hb < min_base_hb: return (bpid = 0)
│  │     │  │     │
│  │     │  │     ├─ Calculate pair type:
│  │     │  │     │  └─ calculate_more_bppars(...) → bpid
│  │     │  │     │
│  │     │  │     └─ Calculate quality:
│  │     │  │        ├─ rtn_val[5] = dorg + 2*dv + plane_angle/20
│  │     │  │        └─ rtn_val[5] += adjust_pairQuality(...)
│  │     │  │
│  │     │  └─ if bpid != 0 && rtn_val[5] < best_score:
│  │     │     ├─ best_score = rtn_val[5]
│  │     │     ├─ pair_istat[1] = j  # Best partner
│  │     │     ├─ pair_istat[2] = bpid
│  │     │     └─ pair_istat[3..] = rtn_val[...]
│  │     │
│  │     └─ Returns: pair_istat[] (best partner, bpid, scores)
│  │
│  ├─ Check mutual match:
│  │  ├─ if pair_istat[1] != 0:  # Found a partner j
│  │  │   ├─ best_pair(pair_istat[1], ..., pair_jstat)  # Check j's best
│  │  │   └─ if i == pair_jstat[1]:  # Mutual match!
│  │  │     ├─ matched_idx[i] = 1
│  │  │     ├─ matched_idx[pair_istat[1]] = 1
│  │  │     ├─ base_pairs[++num_bp][1] = i
│  │  │     └─ For j = 1 to nout:
│  │  │       └─ base_pairs[num_bp][j+1] = pair_istat[j]
│  │
│  └─ (Continue to next residue)
│
├─ Count matched residues:
│  └─ num2 = count(matched_idx == 1)
│
└─ if num2 == num1: break  # No new pairs found

Return: num_bp, base_pairs[][]
```

---

## H-Bond Detection Workflow

### Three-Phase Process

```
get_hbond_ij(i, j, basei, basej, misc_pars, seidx, idx, AtomName, xyz, hb_info)

├─ Phase 1: Initial Detection
│  ├─ num_hbonds = 0
│  ├─ For m in [seidx[i][1], seidx[i][2]]:
│  │   ├─ For n in [seidx[j][1], seidx[j][2]]:
│  │   │   ├─ if good_hbatoms(misc_pars, AtomName[m], AtomName[n], idx[m], idx[n]):
│  │   │   │   └─ if within_limits(xyz[n], xyz[m], hb_lower, hb_dist1):
│  │   │   │     ├─ strcpy(hb_atom1[++num_hbonds], AtomName[m])
│  │   │   │     ├─ strcpy(hb_atom2[num_hbonds], AtomName[n])
│  │   │   │     └─ hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m])
│  │   │   └─ (Continue to next atom pair)
│  │   └─ (Continue to next atom in residue i)
│  └─ Result: num_hbonds initial H-bonds found
│
├─ Phase 2: Conflict Resolution
│  └─ hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, lkg_type, misc_pars)
│     │
│     ├─ Initialize: matched_idx[] = zeros, num_iter = 0
│     │
│     ├─ while not all H-bonds processed:
│     │   ├─ Find best H-bond for each atom (by distance)
│     │   ├─ Check for conflicts (same H-bond best for both atoms)
│     │   ├─ if conflict found:
│     │   │   ├─ hb_dist[k] = -hb_dist[k]  # Mark conflict (negate)
│     │   │   └─ Mark related H-bonds as processed
│     │   └─ num_iter++
│     │
│     ├─ Calculate linkage types:
│     │   ├─ For conflicted H-bonds: idx2[k][1] = 9, idx2[k][2] = 9
│     │   ├─ For non-conflicted sharing atoms: mark appropriately
│     │   └─ lkg_type[k] = idx2[k][1] + idx2[k][2]
│     │
│     └─ Returns: Modified hb_dist[] (negative = conflict), lkg_type[]
│
├─ Phase 3: Validation
│  └─ validate_hbonds(num_hbonds, hb_dist, lkg_type, hb_type, basei, basej, hb_atom1, hb_atom2)
│     │
│     ├─ Initialize: hb_type[k] = ' ' for all k
│     │
│     ├─ For each H-bond k:
│     │   ├─ if hb_dist[k] <= 0.0:  # Conflicted
│     │   │   ├─ hb_dist[k] = fabs(hb_dist[k])  # Make positive
│     │   │   └─ hb_type[k] = donor_acceptor(basei, basej, hb_atom1[k], hb_atom2[k])
│     │   └─ else:
│     │     └─ hb_type[k] = ' '  # Normal
│     │
│     ├─ Filter invalid H-bonds:
│     │   ├─ if hb_dist[k] > 3.6: hb_type[k] = ' '
│     │   └─ if hb_type[k] == '*' && lkg_type != 18 && distance not in [2.6, 3.2]:
│     │     └─ hb_type[k] = ' '
│     │
│     └─ Returns: count(hb_type != ' ')
│
└─ Format output:
   └─ sprintf(hb_info, "[%ld] %s%s %g ...", count, atom1, atom2, dist, ...)
```

**JSON Recording**: 
- `json_writer_record_hbond_list(...)` after Phase 1
- `json_writer_record_hbonds(...)` after Phase 3

---

## Function Dependency Map

### Level 1: Core Mathematical Functions

```
ls_fitting()
├─ cov_matrix()          # Calculate covariance matrix
├─ jacobi()              # Eigenvalue decomposition
├─ dot()                 # Vector dot product
├─ ave_dmatrix()         # Calculate average coordinates
└─ ddxyz()               # Vector difference

get_oarea()
├─ ratom_xyz()           # Get ring atom coordinates
├─ get_zoave()           # Get average z-axis
├─ align2zaxis()        # Align to z-axis
└─ pia_inter()           # Polygon intersection
```

### Level 2: Geometric Calculation Functions

```
base_frame()
├─ ls_fitting()          # Least-squares fitting
├─ read_pdb()            # Read template PDB
├─ find_1st_atom()       # Find atom by name
├─ set_std_base_pdb()    # Construct template filename
└─ mst2orien()           # Convert matrix to array

check_pair()
├─ get_bp_zoave()        # Average z-axis
├─ get_oarea()           # Overlap calculation
├─ good_hbatoms()        # Atom pair validation
├─ within_limits()       # Distance check
├─ get_hbond_ij()        # H-bond detection
├─ calculate_more_bppars()  # Pair type calculation
└─ adjust_pairQuality()  # Quality adjustment
```

### Level 3: Base Pair Finding Functions

```
find_bestpair()
├─ best_pair()
│  └─ check_pair()
│     ├─ get_hbond_ij()
│     │  ├─ hb_atompair()
│     │  │  ├─ update_hb_idx()  # Find best H-bond for atom
│     │  │  └─ Array operations
│     │  └─ validate_hbonds()
│     │     └─ donor_acceptor()  # Determine H-bond type
│     └─ (other geometric checks)
└─ Array management
```

### Level 4: High-Level Workflow Functions

```
handle_str() (find_pair)
├─ read_pdb()            # Parse PDB file
├─ residue_idx()         # Map atoms to residues
├─ get_seq()             # Extract sequence
├─ base_info()
│  └─ base_frame()       # Calculate frames
├─ ring_oidx()           # Index ring atoms
├─ find_bestpair()
│  └─ (dependency chain above)
├─ re_ordering()         # Reorder pairs
└─ x3dna_input()         # Write output

process_str() (analyze)
├─ read_input()          # Read .inp file
├─ read_pdb()            # Parse PDB file
├─ ref_frames()
│  └─ base_frame()       # Recalculate frames
└─ get_parameters()
   ├─ refs_i_j()         # Extract frames
   ├─ bpstep_par()       # Calculate step parameters
   └─ helical_par()      # Calculate helical parameters
```

---

## Data Flow Diagrams

### Data Flow Through find_pair

```
Input: PDB file
  ↓
[PDB Parsing]
  ├─ AtomName[][], ResName[][], xyz[][], etc.
  ↓
[Residue Mapping]
  ├─ seidx[][] (atom ranges per residue)
  ↓
[Sequence Extraction]
  ├─ bseq[] (base sequence)
  ├─ RY[] (purine/pyrimidine classification)
  ↓
[Frame Calculation]
  ├─ orien[][] (rotation matrices)
  ├─ org[][] (origins)
  ├─ NC1xyz[][] (N and C1' coordinates)
  └─ o3_p[][] (O3' and P coordinates)
  ↓
[Ring Atom Indexing]
  ├─ ring_atom[][] (ring atom indices)
  ↓
[Pair Finding]
  ├─ base_pairs[][] (pair data)
  ↓
[Reordering]
  ├─ bp_idx[] (reordered indices)
  ├─ helix_idx[][] (helix regions)
  └─ helix_marker[] (helix boundaries)
  ↓
Output: .inp file
```

### Data Flow Through analyze

```
Input: .inp file + PDB file
  ↓
[Input Parsing]
  ├─ pair_num[][] (residue indices)
  ↓
[PDB Parsing] (same as find_pair)
  ├─ Atom arrays
  ↓
[Frame Recalculation]
  ├─ orien[][] (frames for pairs)
  └─ org[][] (origins for pairs)
  ↓
[Parameter Calculation]
  ├─ bp_step_par[][] (step parameters)
  └─ bp_heli_par[][] (helical parameters)
  ↓
Output: Parameter files
```

---

**Next**: [Helper Functions](05_HELPER_FUNCTIONS.md) for utility function reference

