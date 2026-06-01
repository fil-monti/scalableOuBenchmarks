#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else file.path(getwd(), "benchmarks/benchmarksTrees/bench_pcmfit.R")
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "common.R"))
bt_setup_libpaths()
bt_require("PCMBase")
bt_require("ape")

suppressPackageStartupMessages({
  library(PCMBase)
  library(ape)
})

make_pcmfit_model <- function(case) {
  K <- ncol(case$data)
  model <- PCM("OU", k = K, regimes = 1)
  model$X0[] <- case$params$x0
  model$H[, , 1] <- case$params$A
  model$Theta[, 1] <- case$params$theta
  model$Sigma_x[, , 1] <- t(chol(case$params$diffusion))
  model$Sigmae_x[, , 1] <- 0.0
  model
}

prepare_pcmfit_case <- function(K, N, seed) {
  case <- bt_make_case(K, N, seed)
  tree <- PCMTree(case$tree)
  PCMTreeSetLabels(tree)
  X <- t(case$data)
  colnames(X) <- tree$tip.label
  model <- make_pcmfit_model(case)
  vec <- double(PCMParamCount(model))
  PCMParamLoadOrStore(model, vec, offset = 0L, load = FALSE)
  list(X = X, tree = tree, model = model, vec = vec)
}

pcmfit_loglik <- function(fit) {
  candidates <- list(
    bt_nested(fit, c("logLik")),
    bt_nested(fit, c("LogLik")),
    bt_nested(fit, c("lik")),
    bt_nested(fit, c("value")),
    bt_nested(fit, c("logLikOptim")),
    bt_nested(fit, c("Optim", "value")),
    bt_nested(fit, c("optim", "value")),
    bt_nested(fit, c("BestModel", "logLik")),
    bt_nested(fit, c("BestModel", "LogLik"))
  )
  for (candidate in candidates) {
    if (!is.null(candidate) && length(candidate) && is.finite(candidate[[1L]])) {
      return(as.numeric(candidate[[1L]]))
    }
  }
  if (inherits(try(stats::logLik(fit), silent = TRUE), "logLik")) {
    return(as.numeric(stats::logLik(fit)))
  }
  NA_real_
}

run_pcmfit_one <- function(K, N, seed, warmup, reps) {
  if (!requireNamespace("PCMFit", quietly = TRUE)) {
    return(bt_status_row(
      "PCMFit", K, N, "PCMFit_fit", "skipped",
      "PCMFit is not installed; run INSTALL_PCMFIT=true Rscript install_packages.R"
    ))
  }

  prepared <- prepare_pcmfit_case(K, N, seed)
  fit <- NULL
  maxit <- bt_env_int("PCMFIT_MAXIT", 1L)
  metaI <- if (requireNamespace("PCMBaseCpp", quietly = TRUE)) PCMBaseCpp::PCMInfoCpp else PCMInfo
  matParInit <- matrix(prepared$vec, nrow = 1L)

  timing <- bt_time_expr(function() {
    fit <<- PCMFit::PCMFit(
      prepared$X,
      prepared$tree,
      prepared$model,
      metaI = metaI,
      matParInit = matParInit,
      numRunifInitVecParams = 0L,
      numGuessInitVecParams = 0L,
      numCallsOptim = 1L,
      control = list(maxit = maxit),
      doParallel = FALSE,
      verbose = FALSE
    )
  }, warmup = warmup, reps = reps)
  if (is.null(fit)) {
    fit <- PCMFit::PCMFit(
      prepared$X,
      prepared$tree,
      prepared$model,
      metaI = metaI,
      matParInit = matParInit,
      numRunifInitVecParams = 0L,
      numGuessInitVecParams = 0L,
      numCallsOptim = 1L,
      control = list(maxit = maxit),
      doParallel = FALSE,
      verbose = FALSE
    )
  }

  bt_status_row(
    "PCMFit", K, N, "PCMFit_fit", "ok",
    "fit_time_proxy_PCMBase_wrapper",
    metrics = c(as.list(timing[1L, ]),
                list(logLik = pcmfit_loglik(fit),
                     npar = length(prepared$vec),
                     maxit = maxit))
  )
}

main <- function() {
  grid <- bt_default_grid()
  out <- bt_env_chr("OUT", file.path(bt_script_dir(), "results_pcmfit.csv"))
  bt_remove_if_exists(out)
  for (K in grid$K_values) {
    for (N in grid$N_values) {
      seed <- grid$seed + 4000L * K + N
      row <- tryCatch(
        run_pcmfit_one(K, N, seed, grid$warmup, grid$reps),
        error = function(e) bt_status_row("PCMFit", K, N, "PCMFit_fit", "error", conditionMessage(e))
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
