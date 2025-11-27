# Root Cause Investigation: Why Legacy .par Files Are Not Created

## Problem
When running legacy `analyze_original` with `--no-json`, the `bp_step.par` and `bp_helical.par` files were not created in the working directory.

## Root Causes Found and Fixed

### Issue 1: File writes redirected to `/dev/null` ✅ FIXED

**Location**: `org/src/cmn_fncs.c` (lines 406-425)

The `open_file()` function was redirecting writes to `/dev/null` when JSON writer is initialized:

```c
if (json_writer_is_initialized() && write_mode) {
    if (is .inp file) {
        /* Allow .inp files */
    } else {
        /* Discard writes to other output files */
        fp = fopen("/dev/null", "w");
    }
}
```

**Fix Applied**: Added `.par` files to the exception list (similar to `.inp` files) so they're written to disk even when JSON writer is active.

**Code Change**:
```c
} else if (len >= 4 && strcmp(filename + len - 4, ".par") == 0) {
    /* Allow .par files to be written - needed for parameter comparison */
    fp = fopen(filename, filemode);
}
```

### Issue 2: `--no-json` being treated as input file ✅ FIXED

**Location**: `org/src/analyze.c` (line 602)

When `--no-json` appeared AFTER the input file(s) in the argument list, it was being treated as a second input file: "Processing structure #2: <--no-json>".

**Fix Applied**: Skip arguments that look like flags (start with `-`) in the file processing loop.

**Code Change**:
```c
for (j = i; j < argc; j++) {
    /* Skip arguments that look like flags (start with -) */
    if (argv[j][0] == '-' && strcmp(argv[j], "tmpfile") != 0)
        continue;
    /* Process file */
}
```

## Verification

✅ **Fixed**: Added `.par` file exception to `open_file()` in `cmn_fncs.c`  
✅ **Fixed**: Argument parsing now skips flag-like arguments  
✅ **Verified**: `bp_step.par` and `bp_helical.par` are now created correctly  

## Test Results

After fixes, running:
```bash
org/build/bin/analyze_original 6V9Q_legacy.inp --no-json
```

**Result**: 
- ✅ `bp_step.par` created (1.0 KB)
- ✅ `bp_helical.par` created (1.0 KB)
- ✅ Files contain expected parameter data
- ✅ No errors about `--no-json` being treated as input file

## Next Steps

Now that `.par` files are created, we can:
1. Parse legacy `.par` files in comparison scripts
2. Compare parameters between legacy and modern code
3. Verify parameter matching for perfect pair matches
