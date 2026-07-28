# ============================================================================
# MAIC analysis for HIV host factor integration
# Input: tab-delimited gene lists (category | label | RANKED/UNRANKED | NA | genes...)
# Output: ranked genes, list weights, iteration history
# ============================================================================
# 1. User settings
# Input file.
INPUT_FILE <- "resistance_lists_corrected.txt"

# Folder used to save results.
OUTPUT_DIR <- "MAIC_resistance_final_all_ranked"

# Maximum number of genes retained from each ranked list.
MAX_INPUT_LEN <- 2000L

# MAIC iteration settings.
STABILITY <- 0.01
MAX_ITERATIONS <- 100L
VERBOSE <- TRUE

# If curve fitting fails, a nearly flat descending fallback is used.
FAILED_RANKED_LIST_DESCENT_STEP <- 1e-8


# ============================================================================
# 2. Publication-level category map

CATEGORY_MAP <- data.frame(
  list_label = c(
    "Hiatt_PrimaryT_Base_R_validated",

    "OhAinle_THP1_IFN_R_primary",
    "OhAinle_THP1_ZAPKO_IFN_R_primary_pval",
    "OhAinle_THP1_ZAPKO_VSVG_IFN_R_primary_rank",

    "OhAinle2020_THP1_IFN_N74D_R_primary",
    "OhAinle2020_THP1_IFN_P90A_R_primary",
    "OhAinle2020_THP1_IFN_WT_R_primary",
    "OhAinle2020_THP1_noIFN_N74D_R_primary",

    "Rathore_PrimaryT_Base_R_primary",
    "Zhang_MT4_CRISPRa_R_primary_stars"
  ),

  category = c(
    "Hiatt",

    "OhAinle2018",
    "OhAinle2018",
    "OhAinle2018",

    "OhAinle2020",
    "OhAinle2020",
    "OhAinle2020",
    "OhAinle2020",

    "Rathore",
    "Zhang"
  ),

  stringsAsFactors = FALSE
)


# ============================================================================
# 3. Basic helper functions

normalise_category <- function(x, unknown_index = NULL) {
  x <- trimws(as.character(x))

  if (!nzchar(x)) {
    if (is.null(unknown_index)) {
      stop("A blank category was found.")
    }
    return(sprintf("UNKNOWN-%03d", unknown_index))
  }

  gsub(" +", "-", toupper(x))
}


split_maic_line <- function(line) {
  # Preferred format: tab-delimited.
  columns <- strsplit(line, "\t", fixed = TRUE)[[1]]

  # Fallback: any whitespace delimiter.
  if (length(columns) < 4L) {
    columns <- strsplit(trimws(line), "[[:space:]]+")[[1]]
  }

  columns
}


fallback_rank_weights <- function(
    list_weight,
    n_entities,
    descent_step = FAILED_RANKED_LIST_DESCENT_STEP
) {
  if (n_entities == 0L) {
    return(numeric(0))
  }

  base <- list_weight + (descent_step * n_entities) / 2

  base - (seq_len(n_entities) - 1L) * descent_step
}


# ============================================================================
# 4. Read and validate the input file

read_maic_input <- function(
    file_path,
    max_input_len = 2000L,
    category_map = NULL
) {
  if (!file.exists(file_path)) {
    stop("Input file does not exist: ", file_path)
  }

  raw_lines <- readLines(
    file_path,
    warn = FALSE,
    encoding = "UTF-8"
  )

  valid_lines <- character(0)

  for (line in raw_lines) {
    # Stop at a separator line beginning with at least five hyphens.
    if (grepl("^-----", line)) {
      break
    }

    # Skip blank lines and comments.
    if (grepl("^\\s*$", line) || startsWith(line, "#")) {
      next
    }

    valid_lines <- c(valid_lines, line)
  }

  if (length(valid_lines) == 0L) {
    stop("No valid gene lists were found in the input file.")
  }

  metadata_rows <- vector("list", length(valid_lines))
  membership_rows <- vector("list", length(valid_lines))

  seen_labels <- character(0)
  unknown_counter <- 0L

  for (i in seq_along(valid_lines)) {
    columns <- split_maic_line(valid_lines[[i]])

    if (length(columns) < 4L) {
      stop(
        "Line ", i,
        " contains fewer than four required columns:\n",
        valid_lines[[i]]
      )
    }

    category_raw <- trimws(columns[[1]])
    list_label <- trimws(columns[[2]])
    list_type <- toupper(trimws(columns[[3]]))

    if (!nzchar(list_label)) {
      stop("Line ", i, " has a blank list label.")
    }

    if (list_label %in% seen_labels) {
      stop("Duplicated list label: ", list_label)
    }

    seen_labels <- c(seen_labels, list_label)

    is_ranked <- if (list_type %in% c("RANKED", "R")) {
      TRUE
    } else if (list_type %in% c("UNRANKED", "U")) {
      FALSE
    } else {
      stop(
        "Invalid list type on line ", i, ": ", list_type,
        ". Use RANKED or UNRANKED."
      )
    }

    genes <- if (length(columns) >= 5L) {
      trimws(columns[5:length(columns)])
    } else {
      character(0)
    }

    genes <- genes[nzchar(genes)]

    # Retain only the first occurrence of duplicated genes in a list.
    genes <- genes[!duplicated(genes)]

    if (is_ranked && length(genes) > max_input_len) {
      genes <- genes[seq_len(max_input_len)]
    }

    if (!nzchar(category_raw)) {
      unknown_counter <- unknown_counter + 1L
      category <- normalise_category("", unknown_counter)
    } else {
      category <- normalise_category(category_raw)
    }

    metadata_rows[[i]] <- data.frame(
      input_order = i,
      category = category,
      list_label = list_label,
      is_ranked = is_ranked,
      n_genes = length(genes),
      stringsAsFactors = FALSE
    )

    if (length(genes) > 0L) {
      membership_rows[[i]] <- data.frame(
        list_label = list_label,
        gene = genes,
        rank = seq_along(genes),
        stringsAsFactors = FALSE
      )
    }
  }

  metadata <- do.call(rbind, metadata_rows)

  memberships <- if (all(vapply(
    membership_rows,
    is.null,
    logical(1)
  ))) {
    data.frame(
      list_label = character(0),
      gene = character(0),
      rank = integer(0),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(
      rbind,
      membership_rows[!vapply(
        membership_rows,
        is.null,
        logical(1)
      )]
    )
  }

  # Override input categories using CATEGORY_MAP, if supplied.
  if (!is.null(category_map)) {
    required_columns <- c("list_label", "category")

    if (!all(required_columns %in% colnames(category_map))) {
      stop("CATEGORY_MAP must contain list_label and category columns.")
    }

    if (anyDuplicated(category_map$list_label)) {
      stop("CATEGORY_MAP contains duplicated list labels.")
    }

    matched <- match(
      metadata$list_label,
      category_map$list_label
    )

    has_override <- !is.na(matched)

    if (any(has_override)) {
      replacement <- category_map$category[matched[has_override]]

      if (any(!nzchar(trimws(replacement)))) {
        stop("CATEGORY_MAP contains blank categories.")
      }

      metadata$category[has_override] <- vapply(
        replacement,
        normalise_category,
        character(1)
      )
    }

    missing_in_input <- setdiff(
      category_map$list_label,
      metadata$list_label
    )

    if (length(missing_in_input) > 0L) {
      warning(
        "The following CATEGORY_MAP labels were not found in the input file: ",
        paste(missing_in_input, collapse = ", ")
      )
    }

    unmapped_input <- metadata$list_label[
      !(metadata$list_label %in% category_map$list_label)
    ]

    if (length(unmapped_input) > 0L) {
      warning(
        "The following input lists were not found in CATEGORY_MAP and retained ",
        "their original input categories: ",
        paste(unmapped_input, collapse = ", ")
      )
    }
  }

  if (any(metadata$n_genes == 0L)) {
    warning(
      "Empty lists were found: ",
      paste(
        metadata$list_label[metadata$n_genes == 0L],
        collapse = ", "
      )
    )
  }

  list(
    metadata = metadata,
    memberships = memberships
  )
}


# ============================================================================
# 5. Build internal list objects

build_list_objects <- function(maic_input) {
  metadata <- maic_input$metadata
  memberships <- maic_input$memberships

  list_objects <- vector("list", nrow(metadata))

  for (i in seq_len(nrow(metadata))) {
    current_label <- metadata$list_label[[i]]

    current_members <- memberships[
      memberships$list_label == current_label,
      ,
      drop = FALSE
    ]

    if (nrow(current_members) > 0L) {
      current_members <- current_members[
        order(current_members$rank),
        ,
        drop = FALSE
      ]

      genes <- current_members$gene
    } else {
      genes <- character(0)
    }

    list_objects[[i]] <- list(
      input_order = metadata$input_order[[i]],
      category = metadata$category[[i]],
      list_label = current_label,
      is_ranked = metadata$is_ranked[[i]],
      genes = genes,

      # Initial values.
      weight = 1,
      rank_weights = rep(1, length(genes)),

      fit_ok = NA,
      fit_message = NA_character_,
      fit_parameters = setNames(
        rep(NA_real_, 4),
        c("a", "b", "c", "d")
      ),
      fit_objective = NA_real_
    )
  }

  list_objects
}


# ============================================================================
# 6. Create the list-to-gene contribution matrix

make_contribution_matrix <- function(list_objects, all_genes) {
  contribution_matrix <- matrix(
    NA_real_,
    nrow = length(all_genes),
    ncol = length(list_objects),
    dimnames = list(
      all_genes,
      vapply(
        list_objects,
        `[[`,
        character(1),
        "list_label"
      )
    )
  )

  for (j in seq_along(list_objects)) {
    current_list <- list_objects[[j]]

    if (length(current_list$genes) == 0L) {
      next
    }

    gene_index <- match(
      current_list$genes,
      all_genes
    )

    contribution_matrix[gene_index, j] <-
      current_list$rank_weights
  }

  contribution_matrix
}


# ============================================================================
# 7. Calculate MAIC gene scores

calculate_gene_scores <- function(
    list_objects,
    all_genes,
    return_details = FALSE
) {
  contribution_matrix <- make_contribution_matrix(
    list_objects,
    all_genes
  )

  list_categories <- vapply(
    list_objects,
    `[[`,
    character(1),
    "category"
  )

  categories <- unique(list_categories)

  gene_scores <- setNames(
    numeric(length(all_genes)),
    all_genes
  )

  winner_matrix <- if (return_details) {
    matrix(
      NA_character_,
      nrow = length(all_genes),
      ncol = length(categories),
      dimnames = list(all_genes, categories)
    )
  } else {
    NULL
  }

  category_score_matrix <- if (return_details) {
    matrix(
      0,
      nrow = length(all_genes),
      ncol = length(categories),
      dimnames = list(all_genes, categories)
    )
  } else {
    NULL
  }

  for (category_index in seq_along(categories)) {
    category_name <- categories[[category_index]]

    list_columns <- which(
      list_categories == category_name
    )

    category_matrix <- contribution_matrix[
      ,
      list_columns,
      drop = FALSE
    ]

    has_contribution <- rowSums(
      !is.na(category_matrix)
    ) > 0L

    category_best <- numeric(length(all_genes))

    if (any(has_contribution)) {
      category_best[has_contribution] <- apply(
        category_matrix[
          has_contribution,
          ,
          drop = FALSE
        ],
        1,
        max,
        na.rm = TRUE
      )
    }

    gene_scores <- gene_scores + category_best

    if (return_details) {
      category_score_matrix[, category_index] <-
        category_best

      rows_with_data <- which(has_contribution)

      for (gene_index in rows_with_data) {
        values <- category_matrix[gene_index, ]

        maximum_value <- max(
          values,
          na.rm = TRUE
        )

        # In case of a tie, retain the first list encountered.
        winner_local <- which(
          !is.na(values) &
            values == maximum_value
        )[1]

        winner_matrix[gene_index, category_index] <-
          colnames(category_matrix)[winner_local]
      }
    }
  }

  output <- list(scores = gene_scores)

  if (return_details) {
    output$contribution_matrix <- contribution_matrix
    output$winner_matrix <- winner_matrix
    output$category_score_matrix <- category_score_matrix
  }

  output
}


# ============================================================================
# 8. Fit rank-dependent weights for ranked lists

fit_exponential_rank_weights <- function(
    y_vector,
    list_weight,
    descent_step = FAILED_RANKED_LIST_DESCENT_STEP
) {
  n_entities <- length(y_vector)

  fallback <- function(message) {
    list(
      weights = fallback_rank_weights(
        list_weight = list_weight,
        n_entities = n_entities,
        descent_step = descent_step
      ),
      ok = FALSE,
      message = message,
      parameters = setNames(
        rep(NA_real_, 4),
        c("a", "b", "c", "d")
      ),
      objective = NA_real_
    )
  }

  if (n_entities == 0L) {
    return(list(
      weights = numeric(0),
      ok = TRUE,
      message = "Empty list",
      parameters = setNames(
        rep(NA_real_, 4),
        c("a", "b", "c", "d")
      ),
      objective = NA_real_
    ))
  }

  if (n_entities < 4L) {
    return(fallback("Fewer than four ranked entities"))
  }

  if (!all(is.finite(y_vector))) {
    return(fallback("Non-finite target values"))
  }

  if (!is.finite(list_weight) || list_weight <= 0) {
    return(fallback("Non-positive or non-finite list weight"))
  }

  x_original <- seq_len(n_entities)
  x_scaled <- x_original / n_entities

  sigma <- list_weight *
    (n_entities - x_original + 1) /
    n_entities

  sigma[
    !is.finite(sigma) |
      sigma <= .Machine$double.eps
  ] <- .Machine$double.eps

  y_scaled <- y_vector / sigma

  y_first <- y_vector[[1]]
  y_last <- y_vector[[n_entities]]

  amplitude <- max(
    y_first - y_last,
    0.01
  )

  tail_level <- max(
    y_last,
    1e-6
  )

  decline <- max(
    y_first - y_last,
    1e-6
  )

  # Multiple starting points improve convergence.
  start_list <- list(
    c(
      a = amplitude,
      b = 1,
      c = max(
        decline * 0.05,
        1e-8 * n_entities
      ),
      d = tail_level
    ),

    c(
      a = amplitude,
      b = 5,
      c = max(
        decline * 0.10,
        1e-8 * n_entities
      ),
      d = tail_level
    ),

    c(
      a = amplitude,
      b = 20,
      c = max(
        decline * 0.01,
        1e-8 * n_entities
      ),
      d = tail_level
    ),

    c(
      a = 3,
      b = 0.01 * n_entities,
      c = 0.0001 * n_entities,
      d = list_weight
    )
  )

  lower_bounds <- c(
    a = 0,
    b = 0,
    c = 1e-8 * n_entities,
    d = 0
  )

  upper_bounds <- c(
    a = Inf,
    b = Inf,
    c = Inf,
    d = Inf
  )

  best_fit <- NULL
  best_objective <- Inf
  error_messages <- character(0)

  # --------------------------------------------------------------------------
  # 8.1 First attempt: nls with the port algorithm
  # --------------------------------------------------------------------------

  for (start_values in start_list) {
    fit <- tryCatch(
      suppressWarnings(
        nls(
          y_scaled ~
            (
              a * exp(-b * x_scaled) -
                c * x_scaled +
                d
            ) / sigma,

          start = as.list(start_values),

          algorithm = "port",

          lower = lower_bounds,
          upper = upper_bounds,

          control = nls.control(
            maxiter = 2000,
            tol = 1e-8,
            minFactor = 1e-12,
            warnOnly = TRUE
          )
        )
      ),

      error = function(e) e
    )

    if (inherits(fit, "error")) {
      error_messages <- c(
        error_messages,
        conditionMessage(fit)
      )

      next
    }

    coefficients <- coef(fit)

    fitted_weights <- (
      coefficients[["a"]] *
        exp(-coefficients[["b"]] * x_scaled) -
        coefficients[["c"]] * x_scaled +
        coefficients[["d"]]
    )

    if (
      length(fitted_weights) != n_entities ||
        any(!is.finite(fitted_weights))
    ) {
      next
    }

    weighted_residuals <- (
      fitted_weights - y_vector
    ) / sigma

    objective <- sum(weighted_residuals^2)

    converged <- isTRUE(
      fit$convInfo$isConv
    )

    if (
      converged &&
        is.finite(objective) &&
        objective < best_objective
    ) {
      best_fit <- list(
        fit = fit,
        coefficients = coefficients,
        fitted_weights = fitted_weights,
        objective = objective
      )

      best_objective <- objective
    }
  }

  if (!is.null(best_fit)) {
    coefficients <- best_fit$coefficients

    # Convert b and c back to the original rank scale for reporting.
    original_parameters <- c(
      a = unname(coefficients[["a"]]),
      b = unname(coefficients[["b"]] / n_entities),
      c = unname(coefficients[["c"]] / n_entities),
      d = unname(coefficients[["d"]])
    )

    return(list(
      weights = best_fit$fitted_weights,
      ok = TRUE,
      message = "Converged with nls-port",
      parameters = original_parameters,
      objective = best_fit$objective
    ))
  }

  # --------------------------------------------------------------------------
  # 8.2 Second attempt: nlminb
  # --------------------------------------------------------------------------

  objective_function <- function(parameters) {
    a <- parameters[["a"]]
    b <- parameters[["b"]]
    c_value <- parameters[["c"]]
    d <- parameters[["d"]]

    predicted <- (
      a * exp(-b * x_scaled) -
        c_value * x_scaled +
        d
    )

    if (any(!is.finite(predicted))) {
      return(.Machine$double.xmax / 100)
    }

    weighted_residuals <- (
      predicted - y_vector
    ) / sigma

    value <- sum(weighted_residuals^2)

    if (!is.finite(value)) {
      return(.Machine$double.xmax / 100)
    }

    value
  }

  best_nlminb <- NULL
  best_nlminb_objective <- Inf

  for (start_values in start_list) {
    fit <- tryCatch(
      nlminb(
        start = start_values,
        objective = objective_function,
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(
          iter.max = 5000,
          eval.max = 10000,
          rel.tol = 1e-10,
          x.tol = 1e-10
        )
      ),

      error = function(e) e
    )

    if (inherits(fit, "error")) {
      error_messages <- c(
        error_messages,
        conditionMessage(fit)
      )

      next
    }

    if (
      fit$convergence == 0 &&
        is.finite(fit$objective) &&
        fit$objective < best_nlminb_objective
    ) {
      best_nlminb <- fit
      best_nlminb_objective <- fit$objective
    }
  }

  if (!is.null(best_nlminb)) {
    coefficients <- best_nlminb$par

    fitted_weights <- (
      coefficients[["a"]] *
        exp(-coefficients[["b"]] * x_scaled) -
        coefficients[["c"]] * x_scaled +
        coefficients[["d"]]
    )

    if (
      length(fitted_weights) == n_entities &&
        all(is.finite(fitted_weights))
    ) {
      original_parameters <- c(
        a = unname(coefficients[["a"]]),
        b = unname(coefficients[["b"]] / n_entities),
        c = unname(coefficients[["c"]] / n_entities),
        d = unname(coefficients[["d"]])
      )

      return(list(
        weights = fitted_weights,
        ok = TRUE,
        message = "Converged with nlminb",
        parameters = original_parameters,
        objective = best_nlminb$objective
      ))
    }
  }

  # --------------------------------------------------------------------------
  # 8.3 Fallback only after both fitting methods fail
  # --------------------------------------------------------------------------

  unique_errors <- unique(error_messages)

  fallback(
    paste0(
      "Both nls-port and nlminb failed",
      if (length(unique_errors) > 0L) {
        paste0(
          ": ",
          paste(
            head(unique_errors, 3),
            collapse = " | "
          )
        )
      } else {
        ""
      }
    )
  )
}


# ============================================================================
# 9. Run the MAIC iterations
# ============================================================================

run_maic <- function(
    maic_input,
    stability = 0.01,
    max_iterations = 100L,
    verbose = TRUE
) {
  list_objects <- build_list_objects(maic_input)

  all_genes <- unique(
    unlist(
      lapply(
        list_objects,
        `[[`,
        "genes"
      ),
      use.names = FALSE
    )
  )

  if (length(all_genes) == 0L) {
    stop("No genes were found in the input lists.")
  }

  iteration_history <- vector(
    "list",
    max_iterations
  )

  converged <- FALSE
  final_iteration <- 0L

  for (iteration in seq_len(max_iterations)) {
    # Step 1: update gene scores using current rank weights.
    current_gene_result <- calculate_gene_scores(
      list_objects = list_objects,
      all_genes = all_genes,
      return_details = FALSE
    )

    gene_scores <- current_gene_result$scores

    maximum_delta <- 0

    current_history_rows <- vector(
      "list",
      length(list_objects)
    )

    # Step 2: update each list.
    for (list_index in seq_along(list_objects)) {
      current_list <- list_objects[[list_index]]

      old_weight <- current_list$weight

      current_scores <- gene_scores[
        current_list$genes
      ]

      if (length(current_scores) > 0L) {
        new_weight <- sqrt(
          mean(current_scores)
        )
      } else {
        new_weight <- 0
      }

      current_list$weight <- new_weight

      if (current_list$is_ranked) {
        if (length(current_scores) > 0L) {
          y_vector <- sqrt(
            cumsum(current_scores) /
              seq_along(current_scores)
          )
        } else {
          y_vector <- numeric(0)
        }

        fit_result <- fit_exponential_rank_weights(
          y_vector = y_vector,
          list_weight = new_weight
        )

        current_list$rank_weights <-
          fit_result$weights

        current_list$fit_ok <-
          fit_result$ok

        current_list$fit_message <-
          fit_result$message

        current_list$fit_parameters <-
          fit_result$parameters

        current_list$fit_objective <-
          fit_result$objective
      } else {
        # Unranked list: all genes receive the same contribution.
        current_list$rank_weights <- rep(
          new_weight,
          length(current_list$genes)
        )

        current_list$fit_ok <- NA

        current_list$fit_message <-
          "Unranked list"

        current_list$fit_parameters <- setNames(
          rep(NA_real_, 4),
          c("a", "b", "c", "d")
        )

        current_list$fit_objective <- NA_real_
      }

      delta <- abs(
        new_weight - old_weight
      )

      maximum_delta <- max(
        maximum_delta,
        delta
      )

      top_contribution <- if (
        length(current_list$rank_weights) > 0L
      ) {
        current_list$rank_weights[[1]]
      } else {
        NA_real_
      }

      mean_contribution <- if (
        length(current_list$rank_weights) > 0L
      ) {
        mean(current_list$rank_weights)
      } else {
        NA_real_
      }

      last_contribution <- if (
        length(current_list$rank_weights) > 0L
      ) {
        current_list$rank_weights[[
          length(current_list$rank_weights)
        ]]
      } else {
        NA_real_
      }

      current_history_rows[[list_index]] <- data.frame(
        iteration = iteration,
        category = current_list$category,
        list_label = current_list$list_label,
        is_ranked = current_list$is_ranked,
        base_weight = current_list$weight,
        top_rank_contribution = top_contribution,
        mean_rank_contribution = mean_contribution,
        last_rank_contribution = last_contribution,
        delta = delta,
        max_delta_this_iteration = maximum_delta,
        fit_ok = current_list$fit_ok,
        fit_message = current_list$fit_message,
        stringsAsFactors = FALSE
      )

      list_objects[[list_index]] <- current_list
    }

    iteration_history[[iteration]] <- do.call(
      rbind,
      current_history_rows
    )

    final_iteration <- iteration

    if (verbose) {
      cat(
        sprintf(
          "Iteration %d: maximum list-weight change = %.8f\n",
          iteration,
          maximum_delta
        )
      )
    }

    if (maximum_delta <= stability) {
      converged <- TRUE
      break
    }
  }

  iteration_history <- do.call(
    rbind,
    iteration_history[
      seq_len(final_iteration)
    ]
  )

  # Recalculate final gene scores after the final list update.
  final_gene_result <- calculate_gene_scores(
    list_objects = list_objects,
    all_genes = all_genes,
    return_details = TRUE
  )

  final_scores <- final_gene_result$scores

  contribution_matrix <-
    final_gene_result$contribution_matrix

  winner_matrix <-
    final_gene_result$winner_matrix

  contributors <- apply(
    winner_matrix,
    1,
    function(winners) {
      available <- !is.na(winners)

      if (!any(available)) {
        return("")
      }

      paste(
        paste0(
          colnames(winner_matrix)[available],
          ": ",
          winners[available]
        ),
        collapse = ", "
      )
    }
  )

  n_lists <- rowSums(
    !is.na(contribution_matrix)
  )

  list_categories <- vapply(
    list_objects,
    `[[`,
    character(1),
    "category"
  )

  n_categories <- vapply(
    seq_along(all_genes),
    function(gene_index) {
      present_lists <- which(
        !is.na(
          contribution_matrix[gene_index, ]
        )
      )

      length(
        unique(
          list_categories[present_lists]
        )
      )
    },
    integer(1)
  )

  raw_contributions <- contribution_matrix
  raw_contributions[is.na(raw_contributions)] <- 0

  gene_order <- order(
    -as.numeric(final_scores),
    names(final_scores)
  )

  ordered_genes <- names(final_scores)[
    gene_order
  ]

  gene_summary <- data.frame(
    gene = ordered_genes,
    MAIC_score = as.numeric(
      final_scores[ordered_genes]
    ),
    rank = seq_along(ordered_genes),
    n_lists = n_lists[
      match(ordered_genes, all_genes)
    ],
    n_categories = n_categories[
      match(ordered_genes, all_genes)
    ],
    contributors = contributors[
      match(ordered_genes, all_genes)
    ],
    stringsAsFactors = FALSE
  )

  ordered_contributions <- raw_contributions[
    match(ordered_genes, all_genes),
    ,
    drop = FALSE
  ]

  contribution_df <- as.data.frame(
    ordered_contributions,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  gene_results_wide <- cbind(
    data.frame(
      gene = ordered_genes,
      stringsAsFactors = FALSE
    ),

    contribution_df,

    data.frame(
      MAIC_score = gene_summary$MAIC_score,
      rank = gene_summary$rank,
      n_lists = gene_summary$n_lists,
      n_categories = gene_summary$n_categories,
      contributors = gene_summary$contributors,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )

  list_weight_rows <- vector(
    "list",
    length(list_objects)
  )

  for (list_index in seq_along(list_objects)) {
    current_list <- list_objects[[list_index]]

    weights <- current_list$rank_weights

    list_weight_rows[[list_index]] <- data.frame(
      input_order = current_list$input_order,
      category = current_list$category,
      list_label = current_list$list_label,
      is_ranked = current_list$is_ranked,
      n_genes = length(current_list$genes),

      base_weight = current_list$weight,

      top_rank_contribution =
        if (length(weights) > 0L) {
          weights[[1]]
        } else {
          NA_real_
        },

      mean_rank_contribution =
        if (length(weights) > 0L) {
          mean(weights)
        } else {
          NA_real_
        },

      last_rank_contribution =
        if (length(weights) > 0L) {
          weights[[length(weights)]]
        } else {
          NA_real_
        },

      fit_ok = current_list$fit_ok,
      fit_message = current_list$fit_message,

      fit_a =
        current_list$fit_parameters[["a"]],

      fit_b =
        current_list$fit_parameters[["b"]],

      fit_c =
        current_list$fit_parameters[["c"]],

      fit_d =
        current_list$fit_parameters[["d"]],

      fit_objective =
        current_list$fit_objective,

      stringsAsFactors = FALSE
    )
  }

  list_weights <- do.call(
    rbind,
    list_weight_rows
  )

  list_weights <- list_weights[
    order(list_weights$input_order),
    ,
    drop = FALSE
  ]

  maximum_top_weight <- suppressWarnings(
    max(
      list_weights$top_rank_contribution,
      na.rm = TRUE
    )
  )

  if (
    is.finite(maximum_top_weight) &&
      maximum_top_weight != 0
  ) {
    list_weights$normalised_top_rank_contribution <-
      list_weights$top_rank_contribution /
      maximum_top_weight
  } else {
    list_weights$normalised_top_rank_contribution <-
      NA_real_
  }

  list(
    gene_results = gene_summary,
    gene_results_wide = gene_results_wide,
    list_weights = list_weights,
    iteration_history = iteration_history,
    list_objects = list_objects,
    final_contribution_matrix = contribution_matrix,
    final_winner_matrix = winner_matrix,
    converged = converged,
    iterations = final_iteration,
    stability = stability
  )
}


# ============================================================================
# 10. Export all results
# ============================================================================

write_maic_results <- function(
    maic_results,
    maic_input,
    output_dir
) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.table(
    maic_results$gene_results_wide,
    file = file.path(
      output_dir,
      "maic_raw.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  write.csv(
    maic_results$gene_results,
    file = file.path(
      output_dir,
      "maic_ranked_genes.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    maic_results$list_weights,
    file = file.path(
      output_dir,
      "maic_list_weights.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    maic_results$iteration_history,
    file = file.path(
      output_dir,
      "maic_iteration_history.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    maic_input$metadata,
    file = file.path(
      output_dir,
      "maic_input_summary.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  saveRDS(
    maic_results,
    file = file.path(
      output_dir,
      "maic_complete_results.rds"
    )
  )

  run_summary <- data.frame(
    converged = maic_results$converged,
    iterations = maic_results$iterations,
    stability_threshold = maic_results$stability,
    n_lists = nrow(maic_input$metadata),
    n_categories = length(
      unique(maic_input$metadata$category)
    ),
    n_unique_genes = nrow(
      maic_results$gene_results
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    run_summary,
    file = file.path(
      output_dir,
      "maic_run_summary.csv"
    ),
    row.names = FALSE
  )
}


# ============================================================================
# 11. Run the full analysis
# ============================================================================

maic_input <- read_maic_input(
  file_path = INPUT_FILE,
  max_input_len = MAX_INPUT_LEN,
  category_map = CATEGORY_MAP
)

cat("\n===== Input summary =====\n")
print(maic_input$metadata)

cat(
  "Number of lists:",
  nrow(maic_input$metadata),
  "\n"
)

cat(
  "Number of publication-level categories:",
  length(unique(maic_input$metadata$category)),
  "\n"
)

maic_results <- run_maic(
  maic_input = maic_input,
  stability = STABILITY,
  max_iterations = MAX_ITERATIONS,
  verbose = VERBOSE
)

write_maic_results(
  maic_results = maic_results,
  maic_input = maic_input,
  output_dir = OUTPUT_DIR
)

cat("\n===== Run summary =====\n")
cat(
  "Converged:",
  maic_results$converged,
  "\n"
)

cat(
  "Iterations:",
  maic_results$iterations,
  "\n"
)

cat("\n===== Fit status =====\n")
print(
  maic_results$list_weights[
    ,
    c(
      "category",
      "list_label",
      "n_genes",
      "fit_ok",
      "fit_message",
      "fit_objective"
    )
  ]
)

cat("\n===== Rank-weight ranges =====\n")
print(
  maic_results$list_weights[
    ,
    c(
      "list_label",
      "top_rank_contribution",
      "mean_rank_contribution",
      "last_rank_contribution"
    )
  ]
)

cat("\n===== Top 20 genes =====\n")
print(
  head(
    maic_results$gene_results,
    20
  )
)

cat("\nResults saved to:\n")
cat(
  normalizePath(
    OUTPUT_DIR,
    mustWork = FALSE
  ),
  "\n"
)

