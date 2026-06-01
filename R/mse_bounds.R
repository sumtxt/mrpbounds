#' Compute Tight Bounds on Mean Squared Error
#'
#' @description
#' Computes the minimum and maximum possible values of the mean squared error
#' between bounded variables x_i and known values z_i.
#'
#' @param x_min Numeric vector of lower bounds for x_i
#' @param x_max Numeric vector of upper bounds for x_i  
#' @param z Numeric vector of known values z_i
#'
#' @return Named numeric vector with elements \code{mse_min} and \code{mse_max}
#'
#' @details
#' For each term (x_i - z_i)^2, the function finds its minimum and maximum over
#' x_i ∈ [x_min_i, x_max_i]. Since (x_i - z_i)^2 is a convex parabola with 
#' minimum at x_i = z_i, there are three cases:
#' 
#' \strong{Case 1: z_i ∈ [x_min_i, x_max_i]} (z_i is inside the interval)
#' \itemize{
#'   \item Minimum: 0 (achieved at x_i = z_i)
#'   \item Maximum: max\{(x_min_i - z_i)^2, (x_max_i - z_i)^2\}
#'     (at farther endpoint)
#' }
#' 
#' \strong{Case 2: z_i < x_min_i} (z_i is below the interval)
#' \itemize{
#'   \item Minimum: (x_min_i - z_i)^2 (achieved at x_i = x_min_i)
#'   \item Maximum: (x_max_i - z_i)^2 (achieved at x_i = x_max_i)
#' }
#' 
#' \strong{Case 3: z_i > x_max_i} (z_i is above the interval)
#' \itemize{
#'   \item Minimum: (x_max_i - z_i)^2 (achieved at x_i = x_max_i)
#'   \item Maximum: (x_min_i - z_i)^2 (achieved at x_i = x_min_i)
#' }
#' 
#' The bounds are tight because they can be achieved by choosing each x_i 
#' independently to minimize or maximize its corresponding term.
#'
#' @examples
#' # Simple example
#' x_min <- c(0, -1, 2)
#' x_max <- c(2, 1, 5)
#' z <- c(1, 0, 6)
#' mse_bounds(x_min, x_max, z)
#' 
#' # Case where all z values are within bounds
#' mse_bounds(c(0, 0, 0), c(10, 10, 10), c(5, 3, 7))
#'
#' @export
mse_bounds <- function(x_min, x_max, z) {
  # Check input dimensions
  if (length(x_min) != length(x_max) || length(x_min) != length(z)) {
    stop("x_min, x_max, and z must have the same length")
  }
  
  # Check that x_min <= x_max
  if (any(x_min > x_max)) {
    stop("x_min must be <= x_max for all elements")
  }
  
  n <- length(z)
  min_terms <- numeric(n)
  max_terms <- numeric(n)
  
  for (i in 1:n) {
    if (z[i] >= x_min[i] && z[i] <= x_max[i]) {
      # Case 1: z_i is inside the interval
      min_terms[i] <- 0
      max_terms[i] <- max((x_min[i] - z[i])^2, (x_max[i] - z[i])^2)
    } else if (z[i] < x_min[i]) {
      # Case 2: z_i is below the interval
      min_terms[i] <- (x_min[i] - z[i])^2
      max_terms[i] <- (x_max[i] - z[i])^2
    } else {
      # Case 3: z_i is above the interval
      min_terms[i] <- (x_max[i] - z[i])^2
      max_terms[i] <- (x_min[i] - z[i])^2
    }
  }
  
  # Calculate mean squared error bounds
  mse_min <- mean(min_terms)
  mse_max <- mean(max_terms)
  
  # Return as named vector
  return(c(mse_min = mse_min, mse_max = mse_max))
}