rewards_data <- function(tree_rewards) {
  tree_data <- tree_rewards@data
  
  ac_columns <- grep("^AC[0-9]+_R$", names(tree_data), value = TRUE)
  ac_columns_95 <- grep("^AC[0-9]+_R_0.95_HPD$", names(tree_data), value = TRUE)
  ac_columns_05 <- grep("^AC[0-9]+_R_0.05_HPD$", names(tree_data), value = TRUE)
  drop_cols <- grep("AC", names(tree_data), value = TRUE)
  drop_cols <- setdiff(drop_cols, c(ac_columns, ac_columns_95, ac_columns_05))
  tree_data <- tree_data |> dplyr::select(-all_of(drop_cols))
  
  names(tree_data) <- sub("_R$", "", names(tree_data))
  ac_columns <- sub("_R", "", ac_columns)
  empty_values <- 0
  tree_data <- tree_data %>%
    rowwise() %>%
    mutate(
      # Create a temporary indicator for whether this row is "empty"
      is_empty_row = all(is.na(c_across(all_of(ac_columns)))),
      aircommunity = {
        values <- c_across(all_of(ac_columns))
        if (is_empty_row) { # Use the indicator
          "AC1"
        } else {
          ac_columns[which.max(values)]
        }
      }
    ) %>%
    ungroup()
  
  if (sum(tree_data$is_empty_row) != 1) {
    cat("There should be one empty value (the root) but there are: ", sum(tree_data$is_empty_row), "\n")
  }
  
  trunk_nodes <- nodepath(phy = as.phylo(tree), from = 820, to = 975)
  
  a <- tree_data |>  filter(node %in% trunk_nodes) |> dplyr::select(all_of(ac_columns))
  a[] <- lapply(a, as.numeric)
  summary_stats <- a %>%
    pivot_longer(
      cols = everything(), # Select all columns
      names_to = "variable",
      values_to = "value"
    ) %>%
    group_by(variable) %>%
    summarise(
      median_val = sum(value, na.rm = TRUE),
      lower_ci = quantile(value, 0.05, na.rm = TRUE),
      upper_ci = quantile(value, 0.95, na.rm = TRUE),
      .groups = 'drop' # Drop the grouping after summarising
    )
  summary_stats <- summary_stats |> rename(aircommunity = variable) 
  
  if (any(summary_stats$aircommunity %in% names(ac_to_letter))) { #ac_to_letter defined by the dictionary above
    summary_stats$aircommunity <- ac_to_letter[summary_stats$aircommunity]
  }
  summary_stats$aircommunity <- factor(summary_stats$aircommunity, levels = sort(ac_to_letter, decreasing = FALSE))
  return(summary_stats)
}



rewards_data_trunk <- function(tree_rewards) {
  tree_data <- tree_rewards@data
  trunk_nodes <- nodepath(phy = as.phylo(tree), from = 821, to = 975)
  
  tree_data <- tree_data |>  filter(node %in% trunk_nodes)
  
  ac_columns <- grep("^AC[0-9]+_R$", names(tree_data), value = TRUE)
  ac_columns_hpd <- grep("^AC[0-9]+_R_0.95_HPD$", names(tree_data), value = TRUE)
  
  all_ac_cols <- c(ac_columns, ac_columns_hpd)
  
  tree_data <- tree_data |> dplyr::select(all_of(all_ac_cols))
  names(tree_data) <- sub("_R$", "", names(tree_data))
  names(tree_data) <- sub("_R_", "_", names(tree_data))
  ac_columns <- sub("_R", "", ac_columns)
  ac_columns_hpd <- sub("_R_", "_", ac_columns_hpd)
  
  for (col_name in ac_columns) {
    
    col_name <- paste0(col_name, "_0.95_HPD")
    if (col_name %in% names(tree_data)) {
      
      # Step 1: Replace NA elements with c(0, 0) in the list-column
      tree_data[[col_name]] <- lapply(tree_data[[col_name]], function(x) {
        if (is.na(x[1])) {
          c(0, 0)
        } else {
          x
        }
      })
      
      # Extract the base "ACi" part of the column name (e.g., "AC1" from "AC1_R_0.95_HPD")
      base_ac_name <- sub("_0.95_HPD", "", col_name)
      
      # Step 2: Create new 'lower_HPD_ACi' and 'upper_HPD_ACi' columns
      # We use base R's `$` or `[[]]` for dynamic column assignment in a loop
      tree_data[[paste0("lower_HPD_", base_ac_name)]] <- sapply(tree_data[[col_name]], function(x) x[1])
      tree_data[[paste0("upper_HPD_", base_ac_name)]] <- sapply(tree_data[[col_name]], function(x) x[2])
      
      # Step 3: Remove the original '_R_0.95_HPD' column
      # Use !!sym() to dynamically refer to the column name string within dplyr::select
      tree_data <- tree_data |> dplyr::select(-!!dplyr::sym(col_name))
      
    } else {
      message(paste0("Warning: Column '", col_name, "' not found in tree_data. Skipping."))
    }
  }
  
  n_nodes <- nrow(tree_data)
  tree_data[] <- lapply(tree_data, as.numeric)
  intermediate_data <- colSums(tree_data)
  
  summ_stats <- tibble(
    aircommunity = character(),
    median = numeric(),
    lower_HPD = numeric(),
    upper_HPD = numeric()
  )
  
  # Loop through AC1 to AC14
  for (i in 1:14) {
    ac_name <- paste0("AC", i)
    lower_hpd_col_name <- paste0("lower_HPD_", ac_name)
    upper_hpd_col_name <- paste0("upper_HPD_", ac_name)
    median_value <- intermediate_data[[ac_name]]
    summ_stats <- rbind(summ_stats, 
                        tibble(
                          aircommunity = ac_name,
                          median = median_value,
                          lower_HPD = intermediate_data[[lower_hpd_col_name]],
                          upper_HPD = intermediate_data[[upper_hpd_col_name]]
                        ))
  }
  
  if (any(summ_stats$aircommunity %in% names(ac_to_letter))) { #ac_to_letter defined by the dictionary above
    summ_stats$aircommunity <- ac_to_letter[summ_stats$aircommunity]
  }
  summ_stats$aircommunity <- factor(summ_stats$aircommunity, levels = sort(ac_to_letter, decreasing = FALSE))
  return(summ_stats)
}


loadTreeRewards <- function(basicPath, useBurnIn = FALSE, useThinning = FALSE,
                            useCleaner = FALSE,
                            thinningValue = 10, 
                            headerLinesN = 1650, tailLinesN = 1000) {
  xmlFilePath <- paste0(basicPath, "airCommunitiesMM_rewards.trees")
  temp_path <-  paste0(basicPath, "temp.trees")
  lines <- readLines(xmlFilePath)
  cat("Total number of lines in the file:", length(lines), "\n")
  if (useBurnIn) {
    treesLines <- tail(lines, tailLinesN)
    if (useThinning) {
      treesLines <- treesLines[seq(1, length(treesLines), by = thinningValue)]
    } 
    cat("Analysing", length(treesLines), "trees after burn-in and thinning.\n")
    lines <- c(head(lines, headerLinesN), treesLines)
  }
  if (useCleaner) {
    patterns <- str_c("AC", rep(1:14, each = 2), "_", rep(c("to", "from"), times = 14), "=\\d+\\.?\\d*,")
    full_pattern <- str_c("(", str_c(patterns, collapse = "|"), ")")
    lines <- str_replace_all(lines, full_pattern, "")
    for (i in 1:2) {
      for (j in c("to", "from")) {
        pattern <- str_c(paste0("AC", i, "_", j), "=\\d+\\.?\\d*,")
        lines <- str_replace_all(lines, pattern, "")
      }
    }
  }
  length(lines)
  writeLines(lines, temp_path)
  tailer <- c("END;", "END;")
  cat(tailer, file = temp_path, sep = "\n", append = TRUE)
  full_beast_trees <- read.beast(temp_path)
  return(full_beast_trees)
}