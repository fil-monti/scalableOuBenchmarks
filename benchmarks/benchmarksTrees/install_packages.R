#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- sub("^--file=", "", file_arg[1L] %||% file.path(getwd(), "benchmarks/benchmarksTrees/install_packages.R"))
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
lib <- file.path(script_dir, "r-lib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(normalizePath(lib), .libPaths())))

repos <- c(
  "https://cran.r-universe.dev",
  "https://cloud.r-project.org",
  "https://venelin.r-universe.dev"
)

cran_packages <- c("PCMBase", "PCMBaseCpp", "mvSLOUCH", "mvMORPH", "ape", "phytools", "microbenchmark")
missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, lib = lib, repos = repos)
}

if (tolower(Sys.getenv("INSTALL_PCMFIT", "false")) %in% c("1", "true", "yes")) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = lib, repos = repos)
  }
  if (!requireNamespace("PCMFit", quietly = TRUE)) {
    remotes::install_github("venelin/PCMFit", lib = lib, upgrade = "never")
  }
}

cat("R library path for tree benchmarks:\n", normalizePath(lib), "\n", sep = "")
