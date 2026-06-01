#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else file.path(getwd(), "benchmarks/benchmarksTrees/bench_mvmorph.R")
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "common.R"))
bt_setup_libpaths()
bt_require("mvMORPH")
bt_require("ape")

suppressPackageStartupMessages({
  library(mvMORPH)
  library(ape)
})

prepare_mvmorph_case <- function(K, N, seed) {
  case <- bt_make_case(K, N, seed)
  decomp <- bt_env_chr("MVMORPH_DECOMP", "schur")
  method <- bt_env_chr("MVMORPH_METHOD", "rpf")
  fit <- mvOU(
    tree = case$tree,
    data = case$data,
    model = "OU1",
    param = list(vcv = "fixedRoot", decomp = decomp),
    method = method,
    optimization = "fixed",
    diagnostic = FALSE,
    echo = FALSE
  )
  par <- c(fit$param$alpha, fit$param$sigma)
  list(fit = fit, par = par, decomp = decomp, method = method)
}

run_mvmorph_one <- function(K, N, seed, warmup, reps) {
  prepared <- prepare_mvmorph_case(K, N, seed)
  timing <- bt_time_expr(function() {
    prepared$fit$llik(prepared$par)
  }, warmup = warmup, reps = reps)
  logLik <- as.numeric(prepared$fit$llik(prepared$par))
  bt_status_row(
    "mvMORPH", K, N,
    sprintf("mvOU_fixed_%s_%s", prepared$decomp, prepared$method),
    "ok",
    metrics = c(as.list(timing[1L, ]), list(logLik = logLik, npar = length(prepared$par)))
  )
}

main <- function() {
  grid <- bt_default_grid()
  out <- bt_env_chr("OUT", file.path(bt_script_dir(), "results_mvmorph.csv"))
  bt_remove_if_exists(out)
  for (K in grid$K_values) {
    for (N in grid$N_values) {
      seed <- grid$seed + 2000L * K + N
      row <- tryCatch(
        run_mvmorph_one(K, N, seed, grid$warmup, grid$reps),
        error = function(e) bt_status_row("mvMORPH", K, N, "mvOU_fixed", "error", conditionMessage(e))
      )
      bt_append_result(out, row)
      print(row)
    }
  }
  cat("wrote ", out, "\n", sep = "")
}

if (bt_is_main_script(script_file)) {
  main()
}
