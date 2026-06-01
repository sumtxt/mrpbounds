# Helper functions for mrpbounds package
# These functions are internal and not exported to users

# Validate and normalize input tables
# @param table_list list of data frames
# @param prob_col name of probability column
# @param tol tolerance for normalization
# @return normalized table list
validate_and_normalize_tables <- function(table_list, prob_col, tol) {
  if (!is.list(table_list)) {
    stop("table_list must be a list of data frames")
  }

  normalized_tables <- list()

  for (i in seq_along(table_list)) {
    table <- table_list[[i]]

    if (!is.data.frame(table)) {
      stop(paste("Element", i, "is not a data frame"))
    }

    if (!prob_col %in% names(table)) {
      stop(paste("Probability column '", prob_col, "' not found in table", i))
    }

    if (!is.numeric(table[[prob_col]])) {
      stop(paste("Probability column in table", i, "must be numeric"))
    }

    if (any(table[[prob_col]] < 0)) {
      stop(paste("Negative probabilities found in table", i))
    }

    prob_sum <- sum(table[[prob_col]])
    if (abs(prob_sum - 1) > tol) {
      warning(paste("Probabilities in table", i, "sum to", prob_sum, "- normalizing to 1"))
      table[[prob_col]] <- table[[prob_col]] / prob_sum
    }

    normalized_tables[[i]] <- table
  }

  return(normalized_tables)
}

# Create full table from cartesian product of all variables in table_list
# @param table_list list of normalized data frames
# @param prob_col name of probability column
# @return data frame with all possible combinations of variable levels
create_full_table_from_tables <- function(table_list, prob_col) {
  all_var_data <- list()

  for (table in table_list) {
    var_cols <- setdiff(names(table), prob_col)

    for (var_name in var_cols) {
      var_values <- table[[var_name]]

      if (var_name %in% names(all_var_data)) {
        all_var_data[[var_name]] <- c(all_var_data[[var_name]], var_values)
      } else {
        all_var_data[[var_name]] <- var_values
      }
    }
  }

  unique_var_data <- lapply(all_var_data, function(x) {
    if (is.factor(x)) {
      unique_vals <- unique(x)
      factor(unique_vals, levels = levels(x))
    } else {
      sort(unique(x))
    }
  })

  full_table <- expand.grid(unique_var_data, stringsAsFactors = FALSE)

  return(full_table)
}

# Extract variable information from all tables
# @param table_list list of normalized data frames
# @param prob_col name of probability column
# @return list with variable names and levels
extract_variable_info <- function(table_list, prob_col) {
  all_vars <- list()

  for (table in table_list) {
    var_cols <- setdiff(names(table), prob_col)

    for (var_name in var_cols) {
      levels <- sort(unique(table[[var_name]]))

      if (var_name %in% names(all_vars)) {
        existing_levels <- all_vars[[var_name]]
        if (!identical(levels, existing_levels)) {
          stop(paste(
            "Inconsistent levels for variable", var_name,
            "across tables:", toString(existing_levels), "vs", toString(levels)
          ))
        }
      } else {
        all_vars[[var_name]] <- levels
      }
    }
  }

  return(all_vars)
}

# Build constraint system from input tables
# @param table_list list of normalized data frames
# @param var_info variable information
# @param full_table full cartesian product structure
# @param prob_col name of probability column
# @param tol tolerance for validation
# @return list with Aeq (constraint matrix) and beq (constraint bounds)
build_constraint_system <- function(table_list, var_info, full_table, prob_col, tol) {
  var_names <- names(var_info)
  K <- nrow(full_table)

  constraint_rows <- list()
  constraint_values <- c()

  constraint_rows[["sum_to_one"]] <- rep(1, K)
  constraint_values <- c(constraint_values, 1)

  for (table_idx in seq_along(table_list)) {
    table <- table_list[[table_idx]]
    table_vars <- setdiff(names(table), prob_col)

    if (length(table_vars) == 1) {
      var_name <- table_vars[1]

      for (row_idx in 1:(nrow(table) - 1)) {
        level_value <- table[row_idx, var_name]
        prob_value <- table[row_idx, prob_col]

        matching_cells <- which(full_table[, var_name] == level_value)

        constraint_row <- rep(0, K)
        constraint_row[matching_cells] <- 1

        constraint_name <- paste0(var_name, "_", level_value)
        constraint_rows[[constraint_name]] <- constraint_row
        constraint_values <- c(constraint_values, prob_value)
      }
    } else {
      for (row_idx in 1:(nrow(table) - 1)) {
        prob_value <- table[row_idx, prob_col]

        matching_condition <- rep(TRUE, K)
        for (var_name in table_vars) {
          level_value <- table[row_idx, var_name]
          matching_condition <- matching_condition & (full_table[, var_name] == level_value)
        }

        matching_cells <- which(matching_condition)

        constraint_row <- rep(0, K)
        constraint_row[matching_cells] <- 1

        cell_name <- paste(paste0(table_vars, "=", table[row_idx, table_vars]), collapse = "_")
        constraint_name <- paste0("joint_", cell_name)
        constraint_rows[[constraint_name]] <- constraint_row
        constraint_values <- c(constraint_values, prob_value)
      }
    }
  }

  validate_constraint_consistency(constraint_rows, constraint_values, var_info, tol)

  Aeq <- do.call(rbind, constraint_rows)
  beq <- constraint_values

  return(list(Aeq = Aeq, beq = beq, constraint_names = names(constraint_rows)))
}

# Validate consistency of constraints
# @param constraint_rows list of constraint row vectors
# @param constraint_values vector of constraint bounds
# @param var_info variable information
# @param tol tolerance for consistency check
validate_constraint_consistency <- function(constraint_rows, constraint_values, var_info, tol) {
  var_constraints <- list()

  for (constraint_name in names(constraint_rows)) {
    if (grepl("^joint_", constraint_name)) {
      joint_part <- sub("^joint_", "", constraint_name)
      vars_in_constraint <- stringr::str_extract_all(joint_part, "[A-Za-z]+")[[1]]
    } else if (constraint_name == "sum_to_one") {
      next
    } else {
      vars_in_constraint <- c(stringr::str_extract(constraint_name, "^[^_]+"))
    }

    for (var_name in vars_in_constraint) {
      if (is.null(var_constraints[[var_name]])) {
        var_constraints[[var_name]] <- list()
      }
      var_constraints[[var_name]][[constraint_name]] <- list(
        row = constraint_rows[[constraint_name]],
        value = constraint_values[which(names(constraint_rows) == constraint_name)]
      )
    }
  }

  for (var_name in names(var_constraints)) {
    constraints <- var_constraints[[var_name]]

    marginal_constraints <- constraints[grepl(paste0("^", var_name, "_"), names(constraints))]
    joint_constraints <- constraints[grepl("^joint_", names(constraints))]

    if (length(marginal_constraints) > 0 && length(joint_constraints) > 0) {
      for (marginal_name in names(marginal_constraints)) {
        level_value <- as.numeric(stringr::str_extract(marginal_name, "\\d+$"))
        expected_marginal <- marginal_constraints[[marginal_name]]$value

        actual_marginal <- 0
        for (joint_name in names(joint_constraints)) {
          joint_constraint <- joint_constraints[[joint_name]]
          overlap <- sum(joint_constraint$row * marginal_constraints[[marginal_name]]$row)
          if (overlap > 0) {
            actual_marginal <- actual_marginal + joint_constraint$value
          }
        }

        if (abs(actual_marginal - expected_marginal) > tol) {
          stop(paste(
            "Inconsistent constraints for variable", var_name, "level", level_value,
            ": marginal =", expected_marginal, "but joints sum to", actual_marginal
          ))
        }
      }
    }
  }
}

# Solve linear programming problem using lpSolve
# @param objective objective vector (coefficients to minimize)
# @param Aeq equality constraint matrix
# @param beq equality constraint bounds
# @param bounds_lower lower bounds for variables
# @param bounds_upper upper bounds for variables
# @param sense "minimize" or "maximize"
# @return list with status, solution, and objective_value
solve_linear_program <- function(objective, Aeq, beq, bounds_lower, bounds_upper, sense = "minimize") {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required for optimization. Install with: install.packages('lpSolve')")
  }

  n_vars <- length(objective)
  n_constraints <- nrow(Aeq)

  constraint_matrix <- rbind(Aeq, -Aeq)
  constraint_bounds <- c(beq, -beq)
  constraint_dirs <- rep("<=", 2 * n_constraints)

  if (sense == "maximize") {
    objective <- -objective
  }

  result <- tryCatch(
    {
      if (all(is.infinite(bounds_upper))) {
        lpSolve::lp(
          direction = "min",
          objective.in = objective,
          const.mat = constraint_matrix,
          const.dir = constraint_dirs,
          const.rhs = constraint_bounds
        )
      } else {
        lpSolve::lp(
          direction = "min",
          objective.in = objective,
          const.mat = constraint_matrix,
          const.dir = constraint_dirs,
          const.rhs = constraint_bounds,
          bounds = list(
            lower = list(ind = 1:n_vars, val = bounds_lower),
            upper = list(ind = 1:n_vars, val = bounds_upper)
          )
        )
      }
    },
    error = function(e) {
      return(list(status = 2, message = e$message))
    }
  )

  if (result$status == 0) {
    solution <- result$solution
    obj_value <- result$objval

    if (sense == "maximize") {
      obj_value <- -obj_value
    }

    return(list(
      status = "optimal",
      solution = solution,
      objective_value = obj_value
    ))
  } else {
    status_msg <- switch(as.character(result$status),
      "2" = "no feasible solution",
      "3" = "unbounded solution",
      "infeasible"
    )

    return(list(
      status = status_msg,
      solution = NULL,
      objective_value = NULL
    ))
  }
}
