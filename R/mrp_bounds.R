#' Compute MRP bounds for predicted quantities
#'
#' This function computes bounds for multilevel regression with poststratification
#' (MRP) estimates by finding the minimum and maximum possible values of a
#' target quantity subject to constraints from marginal and joint probability
#' tables.
#'
#' @param model A fitted model object (e.g., from \code{glm}) that can be used
#'   with \code{predict()} to generate predictions on new data.
#' @param table_list A list of data frames, each containing probability
#'   constraints. Each data frame should have columns for categorical variables
#'   and a probability column.
#' @param addon_data Optional data frame merged with the full poststratification
#'   table before prediction. Use when the model requires variables not in
#'   \code{table_list}. Default is \code{NULL}.
#' @param by Character vector of variable names to use as the join key when
#'   merging \code{addon_data} with the full table. Required when
#'   \code{addon_data} is provided. Default is \code{NULL}.
#' @param prob_col Character string specifying the name of the probability
#'   column in all tables. Default is "prob".
#' @param tol Numeric tolerance for optimization and constraint validation.
#'   Default is 1e-8.
#' @param verbose Logical; if \code{TRUE} prints the full poststratification
#'   table and constraint matrix. Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{predict()}.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{min}: A list with \code{status} and \code{objective_value} for
#'     the minimum bound
#'   \item \code{max}: A list with \code{status} and \code{objective_value} for
#'     the maximum bound
#' }
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Validating and normalizing the input probability tables
#'   \item Creating a full multidimensional table from all variable combinations
#'   \item Building a constraint system from the input tables
#'   \item Using linear programming to find bounds on the expected value of
#'     model predictions
#' }
#'
#' The optimization problem is a linear program:
#' \deqn{\min/\max \sum_{i} \hat{y}_i \cdot p_i}
#' subject to the linear constraints defined by \code{table_list}, where
#' \eqn{\hat{y}_i} are the model predictions and \eqn{p_i} are the cell
#' probabilities. Since this is a linear programming problem, the solution
#' is guaranteed to be globally optimal.
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' svy <- data.frame(
#'   A = rbinom(1000, 1, 0.3),
#'   B = rbinom(1000, 1, 0.1),
#'   C = rbinom(1000, 1, 0.5)
#' )
#' model <- glm(A ~ B + C, data = svy, family = binomial())
#'
#' # Define constraints
#' joint_ab <- data.frame(A = c(0, 0, 1, 1), B = c(0, 1, 0, 1),
#'                        prob = c(0.1, 0.2, 0.3, 0.4))
#' margin_c <- data.frame(C = c(0, 1), prob = c(0.6, 0.4))
#' tables <- list(joint_ab, margin_c)
#'
#' # Compute bounds
#' bounds <- mrp_bounds(model, tables)
#' }
#'
#' @export
#' @importFrom stats predict
mrp_bounds <- function(model, table_list, addon_data=NULL, by=NULL, 
  prob_col = "prob", tol = 1e-8, verbose=FALSE, ...) {

  # Input validation
  if (!is.list(table_list)) {
    stop("table_list must be a list of data frames")
  }

  # Validate and normalize input tables
  table_list <- validate_and_normalize_tables(table_list, prob_col, tol)

  # Extract all variables and their levels from all tables
  var_info <- extract_variable_info(table_list, prob_col)
  var_names <- names(var_info)

  full_table <- create_full_table_from_tables(table_list, prob_col)

  if(verbose) {
    print(full_table)
  }

  # Build constraint system from input tables
  constraint_system <- build_constraint_system(table_list, var_info,
                                               full_table, prob_col, tol)

  # Extract constraint matrix and bounds
  aeq_all <- constraint_system$Aeq
  beq_all <- constraint_system$beq

  if(verbose) {
    print(aeq_all)
  }

  # Precompute predictions
  if(!(is.null(addon_data))){
    df <- merge(full_table, addon_data, by=by)
  } else {
    df <- full_table
  }
  yhat <- predict(model, newdata = df, ... , type = "response")
  

  # Solve using linear programming (lpSolve)
  # The problem is: min/max sum(yhat * x) subject to Aeq * x = beq, 0 <= x <= 1
  # for lp, the documtation says: Note that every variable is assumed to be >= 0!
  
  n_vars <- nrow(full_table)
  
  # For minimization (r = 1)
  result_min <- lpSolve::lp(
    direction = "min",
    objective.in = yhat,
    const.mat = aeq_all,
    const.dir = rep("=", nrow(aeq_all)),
    const.rhs = beq_all,
    all.bin = FALSE,
    all.int = FALSE
  )
  
  if(verbose) {
    cat("Minimization result:\n")
    cat("Status:", result_min$status, "\n")
    cat("Objective:", result_min$objval, "\n")
  }
  
  # For maximization (r = -1)
  result_max <- lpSolve::lp(
    direction = "max",
    objective.in = yhat,
    const.mat = aeq_all,
    const.dir = rep("=", nrow(aeq_all)),
    const.rhs = beq_all,
    all.bin = FALSE,
    all.int = FALSE
  )
  
  if(verbose) {
    cat("Maximization result:\n")
    cat("Status:", result_max$status, "\n")
    cat("Objective:", result_max$objval, "\n")
  }

  # Initialize results storage
  mrp_bounds <- data.frame(
    quantity = c("mrp"),
    min_value = numeric(1),
    max_value = numeric(1)
  )

  # Verify constraints are satisfied
  verify_constraints <- function(solution, result_name) {
    if(length(solution) != n_vars) {
      warning(result_name, ": Solution length mismatch")
      return(FALSE)
    }
    
    # Check equality constraints: Aeq * x = beq
    constraint_values <- as.vector(aeq_all %*% solution)
    constraint_violations <- abs(constraint_values - beq_all)
    max_violation <- max(constraint_violations)
    
    if(verbose) {
      cat(result_name, "constraint verification:\n")
      cat("  Max constraint violation:", max_violation, "\n")
      cat("  Tolerance:", tol, "\n")
    }
    
    if(max_violation > tol) {
      warning(result_name, ": Constraint violation detected (max: ", 
              round(max_violation, 8), ", tolerance: ", tol, ")")
      return(FALSE)
    }
    
    # Check bounds: 0 <= x <= 1
    if(any(solution < -tol) || any(solution > 1 + tol)) {
      bound_violations <- pmax(0 - solution, solution - 1, 0)
      max_bound_violation <- max(bound_violations)
      warning(result_name, ": Bound violation detected (max: ", 
              round(max_bound_violation, 8), ")")
      return(FALSE)
    }
    
    if(verbose) {
      cat("  All constraints satisfied\n")
    }
    return(TRUE)
  }
  
  # Process and store results
  mrp_bounds$min_value[1] <- result_min$objval
  mrp_bounds$max_value[1] <- result_max$objval
  
  # Check for optimization errors
  if(result_min$status != 0) {
    warning("Minimization failed with status: ", result_min$status)
  } else {
    verify_constraints(result_min$solution, "Minimization")
  }
  
  if(result_max$status != 0) {
    warning("Maximization failed with status: ", result_max$status)
  } else {
    verify_constraints(result_max$solution, "Maximization")
  }

  return(mrp_bounds)
}
