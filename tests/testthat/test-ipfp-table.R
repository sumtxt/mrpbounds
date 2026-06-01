test_that("ipfp_table basic functionality works", {
  skip_if_not_installed("mipfp")
  
  # Simple joint table with margin
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.1, 0.2, 0.3, 0.4)
  )
  margin_c <- data.frame(
    C = c(0, 1),
    prob = c(0.6, 0.4)
  )
  tables <- list(joint_ab, margin_c)
  
  # Test basic functionality
  result <- ipfp_table(tables)
  
  # Basic checks
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 8)  # 2x2x2 table
  expect_true(all(c("A", "B", "C", "mipfp_prob", "mipfp_count") %in% names(result)))
  
  # Check that probabilities are valid
  expect_true(all(result$mipfp_prob >= 0))
  expect_equal(sum(result$mipfp_prob), 1, tolerance = 1e-6)
  
  # Check that counts are non-negative
  expect_true(all(result$mipfp_count >= 0))
  
  # Check variable levels are correct
  expect_true(all(result$A %in% c(0, 1)))
  expect_true(all(result$B %in% c(0, 1)))
  expect_true(all(result$C %in% c(0, 1)))
})

test_that("ipfp_table handles overlapping constraints", {
  skip_if_not_installed("mipfp")
  
  # Overlapping joint tables A×B and B×C
  ab_table <- data.frame(
    A = c(0, 0, 1, 1), 
    B = c(0, 1, 0, 1), 
    prob = c(0.2, 0.3, 0.3, 0.2)
  )
  bc_table <- data.frame(
    B = c(0, 0, 1, 1), 
    C = c(0, 1, 0, 1), 
    prob = c(0.15, 0.35, 0.25, 0.25)
  )
  tables <- list(ab_table, bc_table)
  
  # Should not throw an error
  expect_no_error(result <- ipfp_table(tables, verbose = FALSE))
  
  # Check structure
  expect_equal(nrow(result), 8)  # 2x2x2 table
  expect_true(all(c("A", "B", "C", "mipfp_prob") %in% names(result)))
  
  # All probabilities should be valid
  expect_true(all(is.finite(result$mipfp_prob)))
  expect_true(all(result$mipfp_prob >= 0))
  expect_equal(sum(result$mipfp_prob), 1, tolerance = 1e-6)
})

test_that("ipfp_table gives consistent solution for unique case", {
  skip_if_not_installed("mipfp")
  
  # Simple case where solution should be unique and known
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.1, 0.2, 0.3, 0.4)  # Fully specified joint distribution
  )
  tables <- list(joint_ab)
  
  # Get IPF solution
  result <- ipfp_table(tables, verbose = FALSE)
  
  # For a fully specified joint table, IPF should recover the original probabilities
  expect_equal(sum(result$mipfp_prob), 1, tolerance = 1e-6)
  
  # Sort both vectors for comparison (order might differ)
  expected_probs <- sort(c(0.1, 0.2, 0.3, 0.4))
  actual_probs <- sort(result$mipfp_prob)
  expect_equal(actual_probs, expected_probs, tolerance = 1e-6)
  
  # Check convergence attributes
  expect_true(attr(result, "converged"))
  expect_equal(attr(result, "method"), "ipfp")
})

test_that("ipfp_table marginal constraints work correctly", {
  skip_if_not_installed("mipfp")
  
  # Test with only marginal constraints
  margin_a <- data.frame(A = c(0, 1), prob = c(0.3, 0.7))
  margin_b <- data.frame(B = c(0, 1), prob = c(0.4, 0.6))
  tables <- list(margin_a, margin_b)
  
  result <- ipfp_table(tables, verbose = FALSE)
  
  # Check that marginals are satisfied
  margin_a_actual <- aggregate(result$mipfp_prob, 
                              by = list(result$A), 
                              FUN = sum)
  expect_equal(margin_a_actual$x, c(0.3, 0.7), tolerance = 1e-6)
  
  margin_b_actual <- aggregate(result$mipfp_prob, 
                              by = list(result$B), 
                              FUN = sum)
  expect_equal(margin_b_actual$x, c(0.4, 0.6), tolerance = 1e-6)
  
  # For independent margins, expect product probabilities
  expected <- c(0.3*0.4, 0.3*0.6, 0.7*0.4, 0.7*0.6)
  expect_equal(sort(result$mipfp_prob), sort(expected), tolerance = 1e-6)
})

test_that("ipfp_table works with different estimation methods", {
  skip_if_not_installed("mipfp")
  
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.25, 0.25, 0.25, 0.25)
  )
  tables <- list(joint_ab)
  
  # Test different methods
  methods <- c("ipfp", "ml", "chi2", "lsq")
  
  for (method in methods) {
    result <- ipfp_table(tables, method = method, verbose = FALSE)
    
    # All methods should give valid probability distributions
    expect_true(all(result$mipfp_prob >= 0), 
                info = paste("Method:", method))
    expect_equal(sum(result$mipfp_prob), 1, tolerance = 1e-6,
                info = paste("Method:", method))
    expect_equal(attr(result, "method"), method,
                info = paste("Method:", method))
    
    # For uniform input, all methods should give uniform output
    expect_true(all(abs(result$mipfp_prob - 0.25) < 1e-6),
                info = paste("Method:", method))
  }
})

test_that("ipfp_table input validation works", {
  skip_if_not_installed("mipfp")
  
  # Test with empty list
  expect_error(ipfp_table(list()), "non-empty list")
  
  # Test with invalid method
  valid_table <- data.frame(A = c(0, 1), prob = c(0.5, 0.5))
  expect_error(
    ipfp_table(list(valid_table), method = "invalid"),
    "method must be one of"
  )
  
  # Test with invalid probability column
  invalid_table <- data.frame(A = c(0, 1), wrong_col = c(0.5, 0.5))
  expect_error(
    ipfp_table(list(invalid_table), prob_col = "prob"),
    "not found"
  )
  
  # Test with negative probabilities
  negative_table <- data.frame(A = c(0, 1), prob = c(-0.1, 1.1))
  expect_error(
    ipfp_table(list(negative_table)),
    "Negative probabilities"
  )
})

test_that("ipfp_table verbose output works", {
  skip_if_not_installed("mipfp")
  
  joint_ab <- data.frame(A = c(0, 1), B = c(0, 1), prob = c(0.3, 0.7))
  tables <- list(joint_ab)
  
  # Test that verbose mode produces output
  expect_output(
    result <- ipfp_table(tables, verbose = TRUE),
    "Variables found"
  )
  
  expect_s3_class(result, "data.frame")
  
  # Test silent mode
  expect_silent(ipfp_table(tables, verbose = FALSE))
})

test_that("ipfp_table handles 3-dimensional tables", {
  skip_if_not_installed("mipfp")
  
  # Create 3D example with joint table and margin
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.2, 0.3, 0.3, 0.2)
  )
  margin_c <- data.frame(
    C = c(0, 1, 2),
    prob = c(0.5, 0.3, 0.2)
  )
  tables <- list(joint_ab, margin_c)
  
  result <- ipfp_table(tables, verbose = FALSE)
  
  # Check dimensions
  expect_equal(nrow(result), 12)  # 2x2x3 table
  expect_true(all(c("A", "B", "C", "mipfp_prob") %in% names(result)))
  
  # Check variable levels
  expect_true(all(result$A %in% c(0, 1)))
  expect_true(all(result$B %in% c(0, 1)))
  expect_true(all(result$C %in% c(0, 1, 2)))
  
  # Check probability validity
  expect_true(all(result$mipfp_prob >= 0))
  expect_equal(sum(result$mipfp_prob), 1, tolerance = 1e-6)
  
  # Check that C marginal is satisfied
  margin_c_actual <- aggregate(result$mipfp_prob, 
                              by = list(result$C), 
                              FUN = sum)
  expect_equal(margin_c_actual$x, c(0.5, 0.3, 0.2), tolerance = 1e-6)
})

test_that("ipfp_table convergence attributes are set correctly", {
  skip_if_not_installed("mipfp")
  
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.1, 0.2, 0.3, 0.4)
  )
  tables <- list(joint_ab)
  
  # Test IPFP method (has convergence info)
  result_ipfp <- ipfp_table(tables, method = "ipfp", verbose = FALSE)
  
  expect_true("converged" %in% names(attributes(result_ipfp)))
  expect_true("iterations" %in% names(attributes(result_ipfp)))
  expect_true("method" %in% names(attributes(result_ipfp)))
  expect_equal(attr(result_ipfp, "method"), "ipfp")
  expect_true(is.logical(attr(result_ipfp, "converged")))
  expect_true(is.numeric(attr(result_ipfp, "iterations")))
  
  # Test other methods (no convergence info but still have method)
  result_ml <- ipfp_table(tables, method = "ml", verbose = FALSE)
  expect_equal(attr(result_ml, "method"), "ml")
})

test_that("ipfp_table handles tolerance parameter", {
  skip_if_not_installed("mipfp")
  
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.25, 0.25, 0.25, 0.25)
  )
  tables <- list(joint_ab)
  
  # Test with different tolerance values
  result_tight <- ipfp_table(tables, tol = 1e-12, verbose = FALSE)
  result_loose <- ipfp_table(tables, tol = 1e-4, verbose = FALSE)
  
  # Both should be valid
  expect_equal(sum(result_tight$mipfp_prob), 1, tolerance = 1e-10)
  expect_equal(sum(result_loose$mipfp_prob), 1, tolerance = 1e-3)
  
  # For this simple case, results should be very similar
  expect_equal(result_tight$mipfp_prob, result_loose$mipfp_prob, tolerance = 1e-3)
})

test_that("ipfp_table comparison with tab_bounds", {
  skip_if_not_installed("mipfp")
  
  # Use same example as tab_bounds tests
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.1, 0.2, 0.3, 0.4)
  )
  margin_c <- data.frame(
    C = c(0, 1),
    prob = c(0.6, 0.4)
  )
  tables <- list(joint_ab, margin_c)
  
  # Get both results
  mipfp_result <- ipfp_table(tables, verbose = FALSE)
  bounds_result <- tab_bounds(tables, verbose = FALSE)
  
  # IPF solution should be within the optimization bounds
  for (i in 1:nrow(mipfp_result)) {
    expect_gte(mipfp_result$mipfp_prob[i], bounds_result$min_bound[i] - 1e-6)
    expect_lte(mipfp_result$mipfp_prob[i], bounds_result$max_bound[i] + 1e-6)
  }
  
  # Both should have same table structure (excluding result columns)
  mipfp_vars <- mipfp_result[, c("A", "B", "C")]
  bounds_vars <- bounds_result[, c("A", "B", "C")]
  expect_equal(mipfp_vars, bounds_vars)
})