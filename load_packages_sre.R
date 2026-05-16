# =============================================================================
# load_packages_sre.R
# =============================================================================
# Run this script OR paste its contents at the top of rsf_homelessness.R.
#
# It tells R where to find the locally installed packages.
# No internet connection or package installation is required.
#
# Instructions
# ------------
#   Option A (recommended):
#     Paste the .libPaths() line at the very top of rsf_homelessness.R,
#     before any library() calls.
#
#   Option B:
#     source("load_packages_sre.R") at the top of rsf_homelessness.R.
# =============================================================================

# Path to the folder containing the stripped package files.
# Update this path if you copied the folder to a different location.
SRE_PKG_DIR <- "installed"   # relative path -- update to absolute if needed
                               # e.g. "C:/research/Homeless/R_packages/installed"

# Add the package folder to R's search path
.libPaths(c(SRE_PKG_DIR, .libPaths()))

# Verify all required packages are available
required <- c("ranger", "readstata13", "dplyr", "ggplot2",
               "survival", "tibble", "tidyr", "scales")

cat("Checking packages in:", SRE_PKG_DIR, "
")
all_ok <- TRUE
for (pkg in required) {
  ok <- requireNamespace(pkg, lib.loc = SRE_PKG_DIR, quietly = TRUE)
  cat(sprintf("  %-20s %s
", pkg, if (ok) "OK" else "MISSING"))
  if (!ok) all_ok <- FALSE
}

if (all_ok) {
  cat("
All packages found. Ready to run rsf_homelessness.R
")
} else {
  cat("
WARNING: Some packages are missing.
")
  cat("Check that the installed/ folder was copied correctly to the SRE.
")
}

