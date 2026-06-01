#' Optimize Cell Bounds for Multidimensional Contingency Tables
#'
#' This function computes the exact minimum and maximum possible values for each
#' cell in a multidimensional contingency table given marginal and joint table
#' constraints.
#'
#' @param table_list A list of data frames, where each data frame represents either:
#'   \itemize{
#'     \item A marginal distribution with columns: variable_name, prob_col
#'     \item A joint table with columns: var1_name, var2_name, ..., prob_col
#'   }
#'   Variable levels should be encoded as integers starting from 0 or 1.
#' @param prob_col Name of the probability column in input tables (default: "prob")
#' @param tol Tolerance for numerical precision (default: 1e-8)
#' @param verbose Logical, whether to print progress information (default: FALSE)
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item One column for each variable with integer values representing categories
#'     \item A "min_bound" column containing the minimum possible probability for each cell
#'     \item A "max_bound" column containing the maximum possible probability for each cell
#'     \item A "range" column containing the difference (max_bound - min_bound)
#'   }
#'
#' @details
#' The function formulates the constraint satisfaction problem as a linear program
#' and solves for the minimum and maximum value of each cell subject to:
#' \itemize{
#'   \item All probabilities are non-negative
#'   \item All probabilities sum to 1
#'   \item All marginal constraints are satisfied
#'   \item All joint table constraints are satisfied
#' }
#' 
#'
#' @examples
#' # Example 1: Using joint table A×B and margin C
#' joint_ab <- data.frame(
#'   A = c(0, 0, 1, 1),
#'   B = c(0, 1, 0, 1), 
#'   prob = c(0.1, 0.2, 0.3, 0.4)
#' )
#' margin_c <- data.frame(
#'   C = c(0, 1),
#'   prob = c(0.6, 0.4)
#' )
#' tables <- list(joint_ab, margin_c)
#' bounds <- tab_bounds(tables)
#' print(bounds)
#'
#' # Example 2: Overlapping tables A×B and B×C
#' ab_table <- data.frame(A = c(0,0,1,1), B = c(0,1,0,1), prob = c(0.2,0.3,0.3,0.2))
#' bc_table <- data.frame(B = c(0,0,1,1), C = c(0,1,0,1), prob = c(0.15,0.35,0.25,0.25))
#' tables <- list(ab_table, bc_table)
#' bounds <- tab_bounds(tables)
#' 
#' # Find cells with tightest constraints
#' tight_cells <- bounds[bounds$range < 0.01, ]
#' print(tight_cells)
#'
#' @export
tab_bounds <- function(table_list, prob_col = "prob", tol = 1e-8, verbose = FALSE) {
  
  # Basic input validation
  if (!is.list(table_list) || length(table_list) == 0) {
    stop("table_list must be a non-empty list of data frames")
  }
  
  # Validate and normalize input tables (reuse existing function)
  table_list <- validate_and_normalize_tables(table_list, prob_col, tol)
  
  # Extract all variables and their levels from all tables
  var_info <- extract_variable_info(table_list, prob_col)
  var_names <- names(var_info)
  var_levels <- sapply(var_info, length)
  
  # Total number of cells in the full multidimensional table
  K <- prod(var_levels)
  
  # Create the full cartesian product structure using actual levels
  level_lists <- lapply(var_info, function(levels) levels)
  full_table <- expand.grid(level_lists)
  names(full_table) <- var_names
  
  # Build constraint system from input tables
  constraint_system <- build_constraint_system(table_list, var_info, full_table, prob_col, tol)
  
  # Extract constraint matrix and bounds
  Aeq <- constraint_system$Aeq
  beq <- constraint_system$beq
  
  if (verbose) {
    cat("Optimizing bounds for", K, "cells with", nrow(Aeq), "constraints\n")
  }
  
  # Check if system has a unique solution (no degrees of freedom)
  nullspace_result <- pracma::nullspace(Aeq)
  if (is.null(nullspace_result) || ncol(nullspace_result) == 0) {
    # Unique solution - min and max are the same
    unique_solution <- as.numeric(pracma::pinv(Aeq) %*% beq)
    full_table$min_bound <- unique_solution
    full_table$max_bound <- unique_solution
    full_table$range <- 0
    
    if (verbose) {
      cat("Unique solution found - all cells have fixed values\n")
    }
    
    return(full_table)
  }
  
  # Solve optimization problems for each cell
  min_bounds <- numeric(K)
  max_bounds <- numeric(K)
  
  for (i in 1:K) {
    if (verbose && i %% max(1, K %/% 10) == 0) {
      cat("Processing cell", i, "of", K, "\n")
    }
    
    # Create objective vector (minimize/maximize cell i)
    obj_min <- rep(0, K)
    obj_min[i] <- 1  # minimize cell i
    
    obj_max <- rep(0, K)
    obj_max[i] <- -1  # maximize cell i (minimize negative)
    
    # Solve for minimum
    min_result <- solve_linear_program(
      objective = obj_min,
      Aeq = Aeq,
      beq = beq,
      bounds_lower = rep(0, K),  # non-negativity
      bounds_upper = rep(Inf, K),
      sense = "minimize"
    )
    
    # Solve for maximum  
    max_result <- solve_linear_program(
      objective = obj_max,
      Aeq = Aeq,
      beq = beq,
      bounds_lower = rep(0, K),  # non-negativity
      bounds_upper = rep(Inf, K),
      sense = "minimize"  # minimizing negative = maximizing
    )
    
    if (min_result$status == "optimal") {
      min_bounds[i] <- min_result$solution[i]
    } else {
      warning(paste("Failed to find minimum for cell", i))
      min_bounds[i] <- NA
    }
    
    if (max_result$status == "optimal") {
      max_bounds[i] <- -max_result$objective_value  # convert back from negative
    } else {
      warning(paste("Failed to find maximum for cell", i))
      max_bounds[i] <- NA
    }
  }
  
  # Add results to output table
  full_table$min_bound <- min_bounds
  full_table$max_bound <- max_bounds
  full_table$range <- max_bounds - min_bounds
  
  if (verbose) {
    cat("Optimization complete\n")
    cat("Average range:", mean(full_table$range, na.rm = TRUE), "\n")
    cat("Cells with range < tol:", sum(full_table$range < tol, na.rm = TRUE), "\n")
  }
  
  return(full_table)
}