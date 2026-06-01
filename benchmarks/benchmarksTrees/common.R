`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

bt_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  frames <- sys.frames()
  files <- vapply(frames, function(env) {
    value <- env$ofile
    if (is.null(value)) "" else value
  }, character(1))
  files <- files[nzchar(files)]
  if (length(files) == 0L) {
    return(getwd())
  }
  dirname(normalizePath(files[[length(files)]], mustWork = FALSE))
}

bt_setup_libpaths <- function() {
  lib <- file.path(bt_script_dir(), "r-lib")
  if (dir.exists(lib)) {
    .libPaths(unique(c(normalizePath(lib), .libPaths())))
  }
  invisible(.libPaths())
}

bt_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Run benchmarks/benchmarksTrees/install_packages.R.", pkg),
         call. = FALSE)
  }
  invisible(TRUE)
}

bt_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(trimws(value))) default else as.integer(value)
}

bt_env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(trimws(value))) default else value
}

bt_remove_if_exists <- function(path) {
  if (file.exists(path)) {
    unlink(path)
  }
  invisible(TRUE)
}

bt_nested <- function(x, path, default = NULL) {
  value <- x
  for (name in path) {
    if (!is.list(value) || is.null(value[[name]])) {
      return(default)
    }
    value <- value[[name]]
  }
  value
}

bt_is_main_script <- function(script_file) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  length(file_arg) &&
    identical(
      normalizePath(script_file, mustWork = FALSE),
      normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
    )
}

bt_parse_ints <- function(value, default) {
  if (!nzchar(trimws(value))) {
    return(default)
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
}

bt_default_grid <- function() {
  list(
    K_values = bt_parse_ints(Sys.getenv("K_VALUES", "2,4"), c(2L, 4L)),
    N_values = bt_parse_ints(Sys.getenv("N_VALUES", "50,100"), c(50L, 100L)),
    reps = bt_env_int("REPS", 10L),
    warmup = bt_env_int("WARMUP", 3L),
    seed = bt_env_int("SEED", 1L)
  )
}

bt_random_orthogonal <- function(K) {
  qr.Q(qr(matrix(stats::rnorm(K * K), K, K)))
}

bt_make_dense_ou <- function(K) {
  stopifnot(K >= 1L)
  R <- bt_random_orthogonal(K)
  D <- matrix(0.0, K, K)
  i <- 1L
  while (i <= K) {
    if (i < K) {
      decay <- stats::runif(1L, 0.35, 1.4)
      u <- stats::runif(1L, 0.25, 0.9)
      v <- -stats::runif(1L, 0.25, 0.9)
      D[i, i] <- decay
      D[i + 1L, i + 1L] <- decay
      D[i, i + 1L] <- u
      D[i + 1L, i] <- v
      i <- i + 2L
    } else {
      D[i, i] <- stats::runif(1L, 0.35, 1.4)
      i <- i + 1L
    }
  }
  A <- R %*% D %*% t(R)
  G <- matrix(stats::rnorm(K * K), K, K)
  V <- crossprod(G) / K + diag(0.25, K)
  list(A = A, diffusion = V, theta = rep(0.0, K), x0 = rep(0.0, K))
}

bt_stationary_variance <- function(A, V) {
  K <- nrow(A)
  lhs <- kronecker(diag(K), A) + kronecker(A, diag(K))
  matrix(solve(lhs, as.vector(V)), K, K)
}

bt_exp_stable <- function(A, dt) {
  ev <- eigen(A)
  F <- ev$vectors %*% diag(exp(-dt * ev$values), nrow(A), nrow(A)) %*% solve(ev$vectors)
  Re(F)
}

bt_simulate_dense_ou_tips <- function(tree, params) {
  bt_require("ape")
  tree <- ape::reorder.phylo(tree, "cladewise")
  K <- length(params$x0)
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  states <- matrix(0.0, n_nodes, K)
  root <- n_tips + 1L
  states[root, ] <- params$x0
  Sigma <- bt_stationary_variance(params$A, params$diffusion)
  Sigma <- (Sigma + t(Sigma)) / 2.0

  for (edge_index in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[edge_index, 1L]
    child <- tree$edge[edge_index, 2L]
    dt <- tree$edge.length[edge_index]
    F <- bt_exp_stable(params$A, dt)
    Q <- Sigma - F %*% Sigma %*% t(F)
    Q <- (Q + t(Q)) / 2.0
    chol_Q <- chol(Q + diag(1.0e-10, K))
    mean_child <- params$theta + as.vector(F %*% (states[parent, ] - params$theta))
    states[child, ] <- mean_child + as.vector(t(chol_Q) %*% stats::rnorm(K))
  }

  tips <- states[seq_len(n_tips), , drop = FALSE]
  rownames(tips) <- tree$tip.label
  colnames(tips) <- paste0("trait", seq_len(K))
  tips
}

bt_make_case <- function(K, N, seed) {
  bt_require("ape")
  set.seed(seed)
  tree <- ape::rtree(N)
  params <- bt_make_dense_ou(K)
  data <- bt_simulate_dense_ou_tips(tree, params)
  list(tree = tree, data = data, params = params)
}

bt_time_expr <- function(expr, warmup = 3L, reps = 10L) {
  force(expr)
  for (i in seq_len(warmup)) {
    expr()
  }
  gc()
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    expr()
    times[[i]] <- proc.time()[["elapsed"]] - t0
  }
  data.frame(
    median_sec = stats::median(times),
    mean_sec = mean(times),
    min_sec = min(times),
    max_sec = max(times)
  )
}

bt_append_result <- function(path, row) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(row, file = path, sep = ",", row.names = FALSE, col.names = !file.exists(path),
              append = file.exists(path), quote = TRUE)
}

bt_status_row <- function(package, K, N, mode, status, message = "", metrics = list()) {
  base <- data.frame(
    package = package,
    K = K,
    N = N,
    mode = mode,
    status = status,
    message = message,
    stringsAsFactors = FALSE
  )
  if (length(metrics)) {
    cbind(base, as.data.frame(metrics, stringsAsFactors = FALSE))
  } else {
    base
  }
}
