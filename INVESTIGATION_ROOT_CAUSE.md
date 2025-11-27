# Root Cause Investigation: Why Legacy .par Files Are Not Created

## Problem
When running legacy `analyze_original` with `--no-json`, the `bp_step.par` and `bp_helical.par` files are not created in the working directory.

## Root Causes Found

### Issue 1: File writes redirected to `/dev/null`

**Location**: `org/src/cmn_fncs.c` (lines 406-425)

The `open_file()` function redirects writes to `/dev/null` when JSON writer is initialized:

```c
if (json_writer_is_initialized() && 
    (filemode[0] == 'w' || filemode[0] == 'a' || strstr(filemode, "w") != NULL)) {
    if (len >= 4 && strcmp(filename + len - 4, ".inp") == 0) {
        /* Allow .inp files */
    } else {
        /* Discard writes to other output files */
        fp = fopen("/dev/null", "w");
    }
}
```

**Fix Applied**: Added `.par` files to the exception list (similar to `.inp` files) so they're written to disk.

### Issue 2: `--no-json` being treated as input file

**Location**: `org/src/analyze.c`

When `--no-json` appears AFTER the input file(s) in the argument list, it's being treated as a second input file and processed as "Processing structure #2: <--no-json>".

**Fix Needed**: Ensure `--no-json` is processed before file iteration begins, OR filter it out from the file list.

## Fix Status

✅ **Fixed**: Added `.par` file exception to `open_file()` in `cmn_fncs.c`  
⚠️ **Pending**: Rebuild legacy binaries to apply the fix  
⚠️ **Pending**: Fix argument parsing so `--no-json` isn't treated as input file

## Next Steps

1. Rebuild legacy analyze binary
2. Test that `.par` files are created
3. Fix argument parsing order issue (if it persists)
