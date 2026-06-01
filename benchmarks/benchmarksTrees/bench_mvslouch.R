#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else file.path(getwd(), "benchmarks/benchmarksTrees/bench_mvslouch.R")
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "common.R"))
bt_setup_libpaths()
bt_require("mvSLOUCH")
bt_require("ape")

suppressPackageStartupMessages({
  library(mvSLOUCH)
  library(ape)
})

bt_parse_pair <- function(value, default) {
  parsed <- bt_parse_ints(value, default)
  if (length(parsed) < 2L || any(is.na(parsed))) {
    return(default)
  }
  parsed[seq_len(2L)]
}

prepare_mvslouch_case <- function(K, N, seed) {
  case <- bt_make_case(K, N, seed)
  list(
    tree = case$tree,
    data = case$data,
    Atype = bt_env_chr("MVSLOUCH_ATYPE", "Any"),
    Syytype = bt_env_chr("MVSLOUCH_SYYTYPE", "Any"),
    diagA = bt_env_chr("MVSLOUCH_DIAGA", "Positive"),
    maxiter = bt_parse_pair(Sys.getenv("MVSLOUCH_MAXITER", "2,2"), c(2L, 2L))
  )
}

run_mvslouch_fit <- function(prepared) {
  ouchModel(
    phyltree = prepared$tree,
    mData = prepared$data,
    Atype = prepared$Atype,
    Syytype = prepared$Syytype,
    diagA = prepared$diagA,
    estimate.root.state = TRUE,
    maxiter = prepared$maxiter
  )
}

mvslouch_loglik <- function(fit) {
  candidates <- list(
    bt_nested(fit, c("FinalFound", "LogLik")),
    bt_nested(fit, c("FinalFound", "ParamSummary", "LogLik")),
    bt_nested(fit, c("MaxLikFound", "LogLik")),
    bt_nested(fit, c("MaxLikFound", "ParamSummary", "LogLik"))
  )
  for (candidate in candidates) {
    if (!is.null(candidate) && length(candidate) && is.finite(candidate[[1L]])) {
      return(as.numeric(candidate[[1L]]))
    }
  }
  NA_real_
}

run_mvslouch_one <- function(K, N, seed, warmup, reps) {
  prepared <- prepare_mvslouch_case(K, N, seed)
  fit <- NULL
  timing <- bt_time_expr(function() {
    fit <<- run_mvslouch_fit(prepared)
  }, warmup = warmup, reps = reps)
  if (is.null(fit)) {
    fit <- run_mvslouch_fit(prepared)
  }
  npar <- length(fit$FinalFound$HeuristicSearchPointFinalFind %||% numeric())
  bt_status_row(
    "mvSLOUCH", K, N,
    sprintf("ouchModel_fit_%s_%s", prepared$Atype, prepared$Syytype),
    "ok",
    "fit_time_proxy_no_public_likelihood_closure",
    metrics = c(as.list(timing[1L, ]),
                list(logLik = mvslouch_loglik(fit),
                     npar = npar,
                     maxiter_outer = prepared$maxiter[[1L]],
                     maxiter_inner = prepared$maxiter[[2L]]))
  )
}

main <- function() {
  grid <- bt_default_grid()
  out <- bt_env_chr("OUT", file.path(bt_script_dir(), "results_mvslouch.csv"))
  bt_remove_if_exists(out)
  for (K in grid$K_values) {
    for (N in grid$N_values) {
      seed <- grid$seed + 3000L * K + N
      row <- tryCatch(
        run_mvslouch_one(K, N, seed, grid$warmup, grid$reps),
        error = function(e) bt_status_row("mvSLOUCH", K, N, "ouchModel_fit", "error", conditionMessage(e))
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
