#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else file.path(getwd(), "benchmarks/benchmarksTrees/generate_beast_xml.R")
source(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), "common.R"))
bt_setup_libpaths()
bt_require("ape")

suppressPackageStartupMessages({
  library(ape)
})

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

fmt_num <- function(x) {
  formatC(as.numeric(x), digits = 10L, format = "fg", flag = "#")
}

collapse_num <- function(x) {
  paste(fmt_num(x), collapse = " ")
}

matrix_parameter_xml <- function(id, M, indent = "    ") {
  lines <- c(sprintf("%s<matrixParameter id=\"%s\">", indent, xml_escape(id)))
  for (j in seq_len(ncol(M))) {
    lines <- c(lines, sprintf("%s    <parameter id=\"%s.col%d\" value=\"%s\"/>",
                              indent, xml_escape(id), j, collapse_num(M[, j])))
  }
  c(lines, sprintf("%s</matrixParameter>", indent))
}

parameter_xml <- function(id, values, indent = "    ", lower = NULL, upper = NULL) {
  attrs <- sprintf("id=\"%s\" value=\"%s\"", xml_escape(id), collapse_num(values))
  if (!is.null(lower)) {
    attrs <- paste(attrs, sprintf("lower=\"%s\"", collapse_num(lower)))
  }
  if (!is.null(upper)) {
    attrs <- paste(attrs, sprintf("upper=\"%s\"", collapse_num(upper)))
  }
  sprintf("%s<parameter %s/>", indent, attrs)
}

givens_matrix <- function(K, angles) {
  R <- diag(K)
  pos <- 1L
  for (i in seq_len(K - 1L)) {
    for (j in seq.int(i + 1L, K)) {
      theta <- angles[[pos]]
      cth <- cos(theta)
      sth <- sin(theta)
      old_i <- R[, i]
      old_j <- R[, j]
      R[, i] <- old_i * cth + old_j * sth
      R[, j] <- -old_i * sth + old_j * cth
      pos <- pos + 1L
    }
  }
  R
}

make_block_case <- function(K, N, seed) {
  if (K %% 2L != 0L) {
    stop("The orthogonal-block BEAST XML case requires even K.", call. = FALSE)
  }
  set.seed(seed)
  tree <- ape::rtree(N)
  n_blocks <- K %/% 2L
  n_angles <- (K * (K - 1L)) %/% 2L
  rho <- seq(0.55, 1.25, length.out = n_blocks)
  block_theta <- seq(-0.18, 0.22, length.out = n_blocks)
  block_t <- seq(-0.30, 0.30, length.out = n_blocks)
  angles <- 0.08 * sin(seq_len(n_angles))

  R <- givens_matrix(K, angles)
  D <- matrix(0.0, K, K)
  for (b in seq_len(n_blocks)) {
    i <- 2L * b - 1L
    cth <- cos(block_theta[[b]])
    sth <- sin(block_theta[[b]])
    D[i, i] <- rho[[b]] * cth
    D[i + 1L, i + 1L] <- rho[[b]] * cth
    D[i, i + 1L] <- rho[[b]] * sth - block_t[[b]]
    D[i + 1L, i] <- rho[[b]] * sth + block_t[[b]]
  }

  G <- matrix(stats::rnorm(K * K), K, K)
  diffusion <- crossprod(G) / K + diag(0.35, K)
  params <- list(
    A = R %*% D %*% t(R),
    diffusion = diffusion,
    theta = rep(0.0, K),
    x0 = rep(0.0, K),
    block_rho = rho,
    block_theta = block_theta,
    block_t = block_t,
    rotation_angles = angles
  )
  data <- bt_simulate_dense_ou_tips(tree, params)
  list(tree = tree, data = data, params = params)
}

taxa_xml <- function(case, trait_name) {
  lines <- c("    <taxa id=\"taxa\">")
  for (taxon in case$tree$tip.label) {
    values <- case$data[taxon, , drop = TRUE]
    lines <- c(lines,
               sprintf("        <taxon id=\"%s\"><attr name=\"%s\">%s</attr></taxon>",
                       xml_escape(taxon), xml_escape(trait_name), collapse_num(values)))
  }
  c(lines, "    </taxa>")
}

tree_model_xml <- function(trait_name, K) {
  c(
    "    <treeModel id=\"treeModel\">",
    "        <newick idref=\"tree\"/>",
    "        <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>",
    "        <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>",
    "        <nodeHeights internalNodes=\"true\" rootNode=\"true\"><parameter id=\"treeModel.allInternalNodeHeights\"/></nodeHeights>",
    sprintf("        <nodeTraits name=\"%s\" rootNode=\"false\" internalNodes=\"false\" leafNodes=\"true\" traitDimension=\"%d\" asMatrix=\"true\">",
            xml_escape(paste0(trait_name, ".tipTraits")), K),
    sprintf("            <parameter id=\"%s.leafTraits\"/>", xml_escape(trait_name)),
    "        </nodeTraits>",
    "    </treeModel>"
  )
}

diffusion_xml <- function(case, trait_name) {
  precision <- solve(case$params$diffusion)
  c(
    sprintf("    <multivariateDiffusionModel id=\"%s.diffusionModel\">", xml_escape(trait_name)),
    "        <precisionMatrix>",
    matrix_parameter_xml(paste0(trait_name, ".precision.matrix"), precision, indent = "            "),
    "        </precisionMatrix>",
    "    </multivariateDiffusionModel>"
  )
}

mean_parameter_xml <- function(trait_name, K) {
  opt_ids <- paste0(trait_name, ".opt.", seq_len(K))
  lines <- unlist(Map(function(id) parameter_xml(id, 0.0, indent = "    "), opt_ids), use.names = FALSE)
  lines <- c(lines, sprintf("    <compoundParameter id=\"%s.meanParameter\">", xml_escape(trait_name)))
  for (id in opt_ids) {
    lines <- c(lines, sprintf("        <parameter idref=\"%s\"/>", xml_escape(id)))
  }
  c(lines, "    </compoundParameter>")
}

optimal_traits_xml <- function(trait_name, K) {
  lines <- sprintf("        <optimalTraits id=\"%s.opt\">", xml_escape(trait_name))
  for (i in seq_len(K)) {
    lines <- c(lines,
               "            <strictClockBranchRates>",
               sprintf("                <rate><parameter idref=\"%s.opt.%d\"/></rate>", xml_escape(trait_name), i),
               "            </strictClockBranchRates>")
  }
  c(lines, "        </optimalTraits>")
}

dense_selection_xml <- function(case, trait_name) {
  c(
    "        <strengthOfSelectionMatrix>",
    matrix_parameter_xml(paste0(trait_name, ".selection.matrix"), case$params$A, indent = "            "),
    "        </strengthOfSelectionMatrix>"
  )
}

block_selection_definitions_xml <- function(case, trait_name, K) {
  n_blocks <- K %/% 2L
  angle_count <- (K * (K - 1L)) %/% 2L
  c(
    parameter_xml(paste0(trait_name, ".blockRho"), case$params$block_rho, indent = "    ", lower = rep(0.0, n_blocks)),
    parameter_xml(paste0(trait_name, ".blockTheta"), case$params$block_theta, indent = "    "),
    parameter_xml(paste0(trait_name, ".blockT"), case$params$block_t, indent = "    "),
    parameter_xml(paste0(trait_name, ".rotation.angle"), case$params$rotation_angles, indent = "    ",
                  lower = rep(-pi, angle_count), upper = rep(pi, angle_count)),
    sprintf("    <givensRotationMatrixParameter id=\"%s.rotation\" dimension=\"%d\">", xml_escape(trait_name), K),
    sprintf("        <angles><parameter idref=\"%s.rotation.angle\"/></angles>", xml_escape(trait_name)),
    "    </givensRotationMatrixParameter>",
    sprintf("    <orthogonalBlockDiagonalPolarStableMatrixParameter id=\"%s.orthBlock\">", xml_escape(trait_name)),
    sprintf("        <orthogonalRotationMatrix><givensRotationMatrixParameter idref=\"%s.rotation\"/></orthogonalRotationMatrix>", xml_escape(trait_name)),
    sprintf("        <blockRho><parameter idref=\"%s.blockRho\"/></blockRho>", xml_escape(trait_name)),
    sprintf("        <blockTheta><parameter idref=\"%s.blockTheta\"/></blockTheta>", xml_escape(trait_name)),
    sprintf("        <blockT><parameter idref=\"%s.blockT\"/></blockT>", xml_escape(trait_name)),
    "    </orthogonalBlockDiagonalPolarStableMatrixParameter>"
  )
}

block_selection_xml <- function(trait_name) {
  c(
    "        <strengthOfSelectionMatrix>",
    sprintf("            <orthogonalBlockDiagonalPolarStableMatrixParameter idref=\"%s.orthBlock\"/>", xml_escape(trait_name)),
    "        </strengthOfSelectionMatrix>"
  )
}

trait_likelihood_xml <- function(case, trait_name, K, selection_mode) {
  selection <- switch(
    selection_mode,
    dense = dense_selection_xml(case, trait_name),
    smbp = block_selection_xml(trait_name),
    stop(sprintf("Unknown BEAST selection mode '%s'", selection_mode), call. = FALSE)
  )

  c(
    sprintf("    <traitDataLikelihood id=\"%s.traitLikelihood\" traitName=\"%s\" selectionChart=\"%s\" implementation=\"canonical\" reconstructTraits=\"false\"",
            xml_escape(trait_name), xml_escape(trait_name),
            if (identical(selection_mode, "dense")) "dense" else "orthogonalBlock"),
    "        useTreeLength=\"false\" scaleByTime=\"true\"",
    "        forceFullPrecision=\"true\" allowSingular=\"true\"",
    "        reportAsMultivariate=\"true\" integrateInternalTraits=\"true\">",
    sprintf("        <multivariateDiffusionModel idref=\"%s.diffusionModel\"/>", xml_escape(trait_name)),
    "        <treeModel idref=\"treeModel\"/>",
    sprintf("        <traitParameter><parameter idref=\"%s.leafTraits\"/></traitParameter>", xml_escape(trait_name)),
    "        <conjugateRootPrior>",
    sprintf("            <meanParameter><compoundParameter idref=\"%s.meanParameter\"/></meanParameter>", xml_escape(trait_name)),
    "            <priorSampleSize><parameter value=\"1.0\"/></priorSampleSize>",
    "        </conjugateRootPrior>",
    optimal_traits_xml(trait_name, K),
    selection,
    "    </traitDataLikelihood>"
  )
}

selection_definitions_xml <- function(case, trait_name, K, selection_mode) {
  switch(
    selection_mode,
    dense = character(),
    smbp = block_selection_definitions_xml(case, trait_name, K),
    stop(sprintf("Unknown BEAST selection mode '%s'", selection_mode), call. = FALSE)
  )
}

benchmark_xml <- function(trait_name, iteration_count) {
  c(
    sprintf("    <benchmarker id=\"%s.benchmark\" iterationCount=\"%d\">", xml_escape(trait_name), iteration_count),
    sprintf("        <traitDataLikelihood idref=\"%s.traitLikelihood\"/>", xml_escape(trait_name)),
    "    </benchmarker>"
  )
}

write_beast_xml <- function(case, K, N, selection_mode, out_dir, iteration_count, seed) {
  trait_name <- sprintf("trait%d", K)
  file <- file.path(out_dir, sprintf("beast_tree_%s_k%d_n%d_seed%d.xml", selection_mode, K, N, seed))
  newick <- ape::write.tree(case$tree)
  lines <- c(
    "<?xml version=\"1.0\" standalone=\"yes\"?>",
    sprintf("<!-- Generated by generate_beast_xml.R: selection=%s, K=%d, N=%d, seed=%d -->",
            selection_mode, K, N, seed),
    "<beast version=\"1.10.4\">",
    "",
    taxa_xml(case, trait_name),
    "",
    sprintf("    <newick id=\"tree\" fixTree=\"true\">%s</newick>", xml_escape(newick)),
    "",
    tree_model_xml(trait_name, K),
    "",
    diffusion_xml(case, trait_name),
    "",
    mean_parameter_xml(trait_name, K),
    "",
    selection_definitions_xml(case, trait_name, K, selection_mode),
    "",
    trait_likelihood_xml(case, trait_name, K, selection_mode),
    "",
    benchmark_xml(trait_name, iteration_count),
    "",
    "</beast>"
  )
  writeLines(lines, file)
  normalizePath(file, mustWork = FALSE)
}

main <- function() {
  grid <- bt_default_grid()
  modes <- strsplit(bt_env_chr("BEAST_XML_MODES", "smbp,dense"), ",", fixed = TRUE)[[1L]]
  modes <- trimws(modes[nzchar(trimws(modes))])
  out_dir <- bt_env_chr("OUT_DIR", file.path(bt_script_dir(), "beast_xml"))
  iteration_count <- bt_env_int("BEAST_ITERATIONS", 1000L)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  written <- character()
  for (K in grid$K_values) {
    for (N in grid$N_values) {
      for (mode in modes) {
        seed <- grid$seed + 5000L * K + N + if (identical(mode, "smbp")) 17L else 0L
        case <- switch(
          mode,
          dense = bt_make_case(K, N, seed),
          smbp = make_block_case(K, N, seed),
          stop(sprintf("Unknown mode '%s'. Use BEAST_XML_MODES=smbp,dense or one mode.", mode), call. = FALSE)
        )
        file <- write_beast_xml(case, K, N, mode, out_dir, iteration_count, seed)
        written <- c(written, file)
        cat("wrote ", file, "\n", sep = "")
      }
    }
  }
  invisible(written)
}

if (bt_is_main_script(script_file)) {
  main()
}
