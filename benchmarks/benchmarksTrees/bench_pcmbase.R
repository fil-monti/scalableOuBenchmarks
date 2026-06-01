#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else file.path(getwd(), "benchmarks/benchmarksTrees/bench_pcmbase.R")
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "common.R"))
bt_setup_libpaths()
bt_require("PCMBase")
bt_require("ape")

suppressPackageStartupMessages({
  library(PCMBase)
  library(ape)
})

make_pcmbase_model <- function(case, model_kind = "dense") {
  K <- ncol(case$data)
  model <- switch(
    model_kind,
    dense = PCM("OU", k = K, regimes = 1),
    schur = PCM(PCMDefaultModelTypes()[["F"]], k = K, regimes = 1),
    stop(sprintf("Unknown PCMBase model_kind '%s'", model_kind), call. = FALSE)
  )
  model$X0[] <- case$params$x0
  model$H[, , 1] <- case$params$A
  model$Theta[, 1] <- case$params$theta
  model$Sigma_x[, , 1] <- t(chol(case$params$diffusion))
  if ("Sigmae_x" %in% names(model)) {
    model$Sigmae_x[, , 1] <- 0.0
  }
  model
}

prepare_pcmbase_case <- function(K, N, seed, use_cpp = FALSE, model_kind = "dense") {
  case <- bt_make_case(K, N, seed)
  tree <- PCMTree(case$tree)
  PCMTreeSetLabels(tree)
  X <- t(case$data)
  colnames(X) <- tree$tip.label
  model <- make_pcmbase_model(case, model_kind = model_kind)
  meta <- if (use_cpp) {
    bt_require("PCMBaseCpp")
    PCMBaseCpp::PCMInfoCpp(X, tree, model)
  } else {
    PCMInfo(X, tree, model)
  }
  vec <- double(PCMParamCount(model))
  PCMParamLoadOrStore(model, vec, offset = 0L, load = FALSE)
  list(X = X, tree = tree, model = model, meta = meta, vec = vec, model_kind = model_kind)
}

run_pcmbase_one <- function(K, N, seed, warmup, reps, use_cpp = FALSE, model_kind = "dense") {
  prepared <- prepare_pcmbase_case(K, N, seed, use_cpp = use_cpp, model_kind = model_kind)
  likelihood <- PCMCreateLikelihood(prepared$X, prepared$tree, prepared$model, metaI = prepared$meta)

  direct <- bt_time_expr(function() {
    PCMLik(prepared$X, prepared$tree, prepared$model, metaI = prepared$meta)
  }, warmup = warmup, reps = reps)

  closure <- bt_time_expr(function() {
    likelihood(prepared$vec)
  }, warmup = warmup, reps = reps)

  mode_prefix <- sprintf("%s_", prepared$model_kind)
  rbind(
    bt_status_row(if (use_cpp) "PCMBaseCpp" else "PCMBase", K, N, paste0(mode_prefix, "PCMLik"), "ok",
                  metrics = c(as.list(direct[1L, ]),
                              list(logLik = as.numeric(PCMLik(prepared$X, prepared$tree, prepared$model, metaI = prepared$meta))))),
    bt_status_row(if (use_cpp) "PCMBaseCpp" else "PCMBase", K, N, paste0(mode_prefix, "PCMCreateLikelihood"), "ok",
                  metrics = c(as.list(closure[1L, ]),
                              list(logLik = as.numeric(likelihood(prepared$vec)))))
  )
}

main <- function() {
  grid <- bt_default_grid()
  out <- bt_env_chr("OUT", file.path(bt_script_dir(), "results_pcmbase.csv"))
  bt_remove_if_exists(out)
  for (K in grid$K_values) {
    for (N in grid$N_values) {
      seed <- grid$seed + 1000L * K + N
      for (use_cpp in c(FALSE, TRUE)) {
        for (model_kind in c("dense", "schur")) {
          row <- tryCatch(
            run_pcmbase_one(K, N, seed, grid$warmup, grid$reps, use_cpp = use_cpp, model_kind = model_kind),
            error = function(e) bt_status_row(if (use_cpp) "PCMBaseCpp" else "PCMBase", K, N, paste0(model_kind, "_likelihood"), "error", conditionMessage(e))
          )
          bt_append_result(out, row)
          print(row)
        }
      }
    }
  }
  cat("wrote ", out, "\n", sep = "")
}

if (bt_is_main_script(script_file)) {
  main()
}
