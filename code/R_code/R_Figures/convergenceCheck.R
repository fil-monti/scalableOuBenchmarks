## Convergence check for simulation_GP_exact run
## Usage: source this file or run via Rscript convergenceCheck.R
## Output: prints ESS table, saves trace and ACF PDFs

log_path <- "/Users/filippomonti/Desktop/npCTMCBenchmarks/simulationGP/results/simulationGP_exact/simulation_GP_exact.log"
output_dir <- "/Users/filippomonti/Desktop/npCTMCBenchmarks/simulationGP/results/simulationGP_exact"

burn_in <- 0.1   # discard first 10%

# --- source infrastructure ---
stopifnot(requireNamespace("rprojroot", quietly = TRUE))
stopifnot(requireNamespace("fs",        quietly = TRUE))
root <- rprojroot::find_root(rprojroot::has_dir("code"),
                             path = "/Users/filippomonti/Desktop/suchard_lab/NonParametricModelingOfCTMCs")
R_path <- fs::path(root, "code", "R_code")
sapply(file.path(R_path, c("R_libraries.R", "R_paths.R")), source)
r_files <- c(
  list.files(R_functions_path, pattern = "\\.R$", full.names = TRUE, recursive = TRUE),
  list.files(R_classes_path,   pattern = "\\.R$", full.names = TRUE)
)
for (f in r_files) source(f)

# --- load log ---
cat("Loading log:", log_path, "\n")
PJL <- AnalyserFromLog$new(
  project            = "NPRates",
  path               = log_path,
  jobsAttributes     = list(),
  additionalAttributes = list(),
  burnIn             = burn_in,
  thinning           = 1,
  burnOut            = 1,
  actions            = list()
)

mcmcEval <- PJL$mcmc(1)
n_post   <- nrow(mcmcEval$data)
cat(sprintf("Post-burnin samples: %d\n", n_post))

# --- ESS ---
cat("\n=== ESS (sorted ascending) ===\n")
ess_df <- mcmcEval$ESSSorted()
print(head(ess_df, 30))

# summary of low-ESS parameters
low_ess <- ess_df[ess_df$ESS > 0 & ess_df$ESS < 200, ]
if (nrow(low_ess) > 0) {
  cat("\nParameters with ESS < 200:\n")
  print(low_ess)
} else {
  cat("\nAll parameters have ESS >= 200.\n")
}

# --- trace plots: key scalars ---
key_cols <- c("joint", "prior", "likelihood", "hyperScaleExp1", "hyperLengthExp1", "clock.rate")
key_cols <- intersect(key_cols, names(mcmcEval$data))
traces_scalar <- mcmcEval$plottingTraces(ncol = 3, cols = match(key_cols, names(mcmcEval$data)))
out_trace_scalar <- file.path(output_dir, "traces_scalars.pdf")
ggsave(out_trace_scalar, traces_scalar, width = 14, height = 8)
cat("\nSaved:", out_trace_scalar, "\n")

# --- trace plots: 12 representative host.model entries ---
hm_cols <- grep("^host\\.model\\.", names(mcmcEval$data), value = TRUE)
set.seed(42)
hm_sample <- sort(sample(hm_cols, min(12, length(hm_cols))))
traces_hm <- mcmcEval$plottingTraces(ncol = 4, cols = match(hm_sample, names(mcmcEval$data)))
out_trace_hm <- file.path(output_dir, "traces_hostModel.pdf")
ggsave(out_trace_hm, traces_hm, width = 16, height = 10)
cat("Saved:", out_trace_hm, "\n")

# --- ACF plots: key scalars ---
acf_scalar <- mcmcEval$plottingAcf(n_cols = length(key_cols), n_lags = 50)
out_acf_scalar <- file.path(output_dir, "acf_scalars.pdf")
ggsave(out_acf_scalar, acf_scalar, width = 14, height = 8)
cat("Saved:", out_acf_scalar, "\n")
