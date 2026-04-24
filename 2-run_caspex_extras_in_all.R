# =============================================================================
# 2-run_caspex_extras_in_all.R — helper to regenerate the supplementary plots
# for an existing run_caspex() result in the current R session without
# re-running the primary pipeline.
#
# Expects `result` to be in the workspace (e.g. produced by 1-try.R or
# 1-try-hocomoco.R). Writes A1-A6, B1-B3, C1-C3, D1-D3 into result$out_dir/extras
# (the NNLS-only B4/C4 steps auto-skip under coverage-aware runs).
# =============================================================================

source("caspex_extras.R")

if (!exists("result")) {
  stop("`result` not found in the workspace. Run 1-try.R (or 1-try-hocomoco.R) ",
       "first so that `result` is defined, then source this script.")
}

out_dir <- file.path(result$out_dir %||% "caspex_output", "extras")
run_caspex_extras(result, out_dir = out_dir)
