#' Fit Multidimensional Contingency Tables using Iterative Proportional Fitting
#'
#' This function estimates multidimensional contingency table cell probabilities
#' using the Iterative Proportional Fitting Procedure (IPFP) via the mipfp package.
#' It takes the same input format as tab_bounds() but returns a single table instead of 
#' bounds. 
#'
#' @param table_list A list of data frames, where each data frame represents either:
#'   \itemize{
#'     \item A marginal distribution with columns: variable_name, prob_col
#'     \item A joint table with columns: var1_name, var2_name, ..., prob_col
#'   }
#'   Variable levels should be encoded as integers starting from 0 or 1.
#' @param prob_col Name of the probability column in input tables (default: "prob")
#' @param tol Tolerance for numerical precision (default: 1e-8)
#' @param method Estimation method: "ipfp" (default), "ml", "chi2", or "lsq"
#' @param max_iter Maximum number of iterations for IPFP (default: 1000)
#' @param verbose Logical, whether to print progress information (default: FALSE)
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item One column for each variable with integer values representing categories
#'     \item A "mipfp_prob" column containing the fitted probabilities for each cell
#'     \item A "mipfp_count" column containing the fitted counts (if different from probabilities)
#'   }
#'
#' @details
#' The function converts the input table_list format to the format expected by mipfp,
#' runs the iterative proportional fitting procedure, and converts the results back
#' to the standard format.
#' 
#' The IPF algorithm iteratively adjusts cell values to match all marginal constraints
#' simultaneously, converging to the maximum likelihood estimate under the assumption
#' of independence.
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
#' result <- ipfp_table(tables)
#' print(result)
#'
#' @export
#' @importFrom utils tail
ipfp_table <- function(table_list, prob_col = "prob", tol = 1e-8, 
                          method = "ipfp", max_iter = 1000, verbose = FALSE) {

  # Check if mipfp package is available
  if (!requireNamespace("mipfp", quietly = TRUE)) {
    stop("Package 'mipfp' is required but not installed. Please install it with: install.packages('mipfp')")
  }
  
  # Basic input validation
  if (!is.list(table_list) || length(table_list) == 0) {
    stop("table_list must be a non-empty list of data frames")
  }
  
  if (!method %in% c("ipfp", "ml", "chi2", "lsq")) {
    stop("method must be one of: 'ipfp', 'ml', 'chi2', 'lsq'")
  }
  
  # Validate and normalize input tables (reuse existing function)
  table_list <- validate_and_normalize_tables(table_list, prob_col, tol)
  
  # Extract all variables and their levels from all tables
  var_info <- extract_variable_info(table_list, prob_col)
  var_names <- names(var_info)
  var_levels <- sapply(var_info, length)
  
  if (verbose) {
    cat("Variables found:", toString(var_names), "\n")
    cat("Variable levels:", toString(var_levels), "\n")
  }
  
  # Total number of cells in the full multidimensional table
  K <- prod(var_levels)
  
  # Create the full cartesian product structure using actual levels
  level_lists <- lapply(var_info, function(levels) levels)
  full_table <- expand.grid(level_lists)
  names(full_table) <- var_names
  
  if (verbose) {
    cat("Full table dimensions:", toString(var_levels), "=", K, "cells\n")
  }
  
  # Convert to mipfp format
  mipfp_data <- convert_to_mipfp_format(table_list, var_info, var_names, prob_col, verbose)
  
  # Create seed array (uniform start)
  seed <- array(1, dim = var_levels)
  
  if (verbose) {
    cat("Created seed array with dimensions:", toString(dim(seed)), "\n")
    cat("Number of constraints:", length(mipfp_data$target_list), "\n")
  }
  
  # Run mipfp estimation
  if (method == "ipfp") {
    # Use Iterative Proportional Fitting
    result <- mipfp::Ipfp(
      seed = seed,
      target.list = mipfp_data$target_list,
      target.data = mipfp_data$target_data,
      iter = max_iter,
      tol = tol,
      tol.margins = tol,
      print = verbose
    )
    
    if (!result$conv && verbose) {
      warning("IPFP did not converge within ", max_iter, " iterations")
    }
    
    fitted_probs <- as.vector(result$p.hat)
    fitted_counts <- as.vector(result$x.hat)
    
  } else {
    # Use alternative estimation methods
    result <- mipfp::ObtainModelEstimates(
      seed = seed,
      target.list = mipfp_data$target_list,
      target.data = mipfp_data$target_data,
      method = method,
      tol.margins = tol
    )
    
    fitted_probs <- as.vector(result$p.hat)
    fitted_counts <- as.vector(result$x.hat)
  }
  
  if (verbose) {
    cat("Estimation completed successfully\n")
    cat("Sum of fitted probabilities:", sum(fitted_probs), "\n")
    if (method == "ipfp") {
      cat("Converged:", result$conv, "\n")
      if (!is.null(result$error.margins)) {
        cat("Max margin error:", max(abs(result$error.margins)), "\n")
      }
    }
  }
  
  # Add results to output table
  full_table$mipfp_prob <- fitted_probs
  full_table$mipfp_count <- fitted_counts
  
  # Add convergence information as attributes
  if (method == "ipfp") {
    attr(full_table, "converged") <- result$conv
    attr(full_table, "iterations") <- length(result$evol.stp.crit)
    attr(full_table, "final_criterion") <- tail(result$evol.stp.crit, 1)
    if (!is.null(result$error.margins)) {
      attr(full_table, "margin_errors") <- result$error.margins
    }
  }
  attr(full_table, "method") <- method
  
  return(full_table)
}

# Convert table_list format to mipfp format
# 
# Internal function to convert the mrpbounds table_list format to the
# target_list and target_data format expected by mipfp functions.
#
# @param table_list list of input tables
# @param var_info variable information from extract_variable_info
# @param var_names vector of variable names
# @param prob_col name of probability column
# @param verbose logical for progress output
# @return list with target_list and target_data components
convert_to_mipfp_format <- function(table_list, var_info, var_names, prob_col, verbose = FALSE) {
  
  target_list <- list()
  target_data <- list()
  
  for (i in seq_along(table_list)) {
    table <- table_list[[i]]
    table_vars <- setdiff(names(table), prob_col)
    
    if (verbose) {
      cat("Processing table", i, "with variables:", toString(table_vars), "\n")
    }
    
    # Map variable names to dimension indices
    var_dims <- match(table_vars, var_names)
    
    if (length(table_vars) == 1) {
      # Marginal constraint
      var_name <- table_vars[1]
      var_levels <- var_info[[var_name]]
      
      # Create target array for this marginal
      target_array <- numeric(length(var_levels))
      
      for (j in 1:nrow(table)) {
        level_value <- table[j, var_name]
        level_index <- match(level_value, var_levels)
        if (is.na(level_index)) {
          stop(paste("Level", level_value, "not found in variable", var_name))
        }
        target_array[level_index] <- table[j, prob_col]
      }
      
      target_list[[length(target_list) + 1]] <- var_dims[1]
      target_data[[length(target_data) + 1]] <- target_array
      
      if (verbose) {
        cat("  Marginal for dimension", var_dims[1], ":", toString(round(target_array, 3)), "\n")
      }
      
    } else {
      # Joint constraint
      # Determine dimensions for this joint table
      joint_dims <- sapply(table_vars, function(v) length(var_info[[v]]))
      
      # Create multi-dimensional target array
      target_array <- array(0, dim = joint_dims)
      
      # Fill the target array
      for (j in 1:nrow(table)) {
        # Get indices for this row
        indices <- numeric(length(table_vars))
        for (k in seq_along(table_vars)) {
          var_name <- table_vars[k]
          level_value <- table[j, var_name]
          level_index <- match(level_value, var_info[[var_name]])
          if (is.na(level_index)) {
            stop(paste("Level", level_value, "not found in variable", var_name))
          }
          indices[k] <- level_index
        }
        
        # Set the probability value
        target_array[matrix(indices, nrow = 1)] <- table[j, prob_col]
      }
      
      target_list[[length(target_list) + 1]] <- var_dims
      target_data[[length(target_data) + 1]] <- target_array
      
      if (verbose) {
        cat("  Joint for dimensions", toString(var_dims), "with", prod(joint_dims), "cells\n")
      }
    }
  }
  
  return(list(target_list = target_list, target_data = target_data))
}