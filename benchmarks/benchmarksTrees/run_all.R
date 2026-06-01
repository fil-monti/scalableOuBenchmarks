#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else file.path(getwd(), "benchmarks/benchmarksTrees/run_all.R")
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

scripts <- c(
  "bench_pcmbase.R",
  "bench_mvmorph.R",
  "bench_mvslouch.R",
  "bench_pcmfit.R"
)

for (script in scripts) {
  path <- file.path(script_dir, script)
  cat("\n== ", script, " ==\n", sep = "")
  status <- system2(file.path(R.home("bin"), "Rscript"), path)
  if (!identical(status, 0L)) {
    warning(sprintf("%s exited with status %s", script, status), call. = FALSE)
  }
}
