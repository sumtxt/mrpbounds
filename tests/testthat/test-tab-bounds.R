test_that("tab_bounds basic functionality works", {
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
  
  # Test optimization
  bounds <- tab_bounds(tables)
  
  # Basic checks
  expect_s3_class(bounds, "data.frame")
  expect_equal(nrow(bounds), 8)  # 2x2x2 table
  expect_true(all(c("A", "B", "C", "min_bound", "max_bound", "range") %in% names(bounds)))
  
  # Check bounds are valid
  expect_true(all(bounds$min_bound >= 0, na.rm = TRUE))
  expect_true(all(bounds$max_bound >= bounds$min_bound, na.rm = TRUE))
  expect_true(all(bounds$range >= 0, na.rm = TRUE))
  
  # Check that sum of max bounds is >= 1 (since probabilities must sum to 1)
  expect_gte(sum(bounds$max_bound, na.rm = TRUE), 1)
})

test_that("tab_bounds handles overlapping constraints", {
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
  expect_no_error(bounds <- tab_bounds(tables, verbose = FALSE))
  
  # Check structure
  expect_equal(nrow(bounds), 8)  # 2x2x2 table
  expect_true(all(c("min_bound", "max_bound", "range") %in% names(bounds)))
  
  # All bounds should be valid
  expect_true(all(is.finite(bounds$min_bound)))
  expect_true(all(is.finite(bounds$max_bound)))
  expect_true(all(bounds$min_bound >= 0))
  expect_true(all(bounds$max_bound >= bounds$min_bound))
})

test_that("tab_bounds gives reasonable bounds", {
  # Simple case for testing bounds
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.25, 0.25, 0.25, 0.25)  # Uniform joint distribution
  )
  tables <- list(joint_ab)
  
  # Get optimization bounds
  bounds <- tab_bounds(tables, verbose = FALSE)
  
  # For a uniform distribution, expect tight bounds (unique solution)
  expect_true(all(bounds$range < 1e-6))
  expect_true(all(bounds$min_bound == bounds$max_bound))
  
  # Check that probabilities sum to 1
  expect_equal(sum(bounds$min_bound), 1, tolerance = 1e-6)
})

test_that("tab_bounds handles unique solutions", {
  # Create a system with a unique solution (fully determined)
  # Use a single joint table with all combinations specified
  joint_ab <- data.frame(
    A = c(0, 0, 1, 1),
    B = c(0, 1, 0, 1), 
    prob = c(0.1, 0.2, 0.3, 0.4)
  )
  tables <- list(joint_ab)
  
  bounds <- tab_bounds(tables, verbose = FALSE)
  
  # When there's a unique solution, min_bound should equal max_bound
  expect_true(all(abs(bounds$max_bound - bounds$min_bound) < 1e-6))
  expect_true(all(bounds$range < 1e-6))
  
  # The unique solution should sum to 1 and match the input probabilities (in some order)
  expect_equal(sum(bounds$min_bound), 1, tolerance = 1e-6)
  expect_setequal(round(bounds$min_bound, 6), c(0.1, 0.2, 0.3, 0.4))
})

test_that("tab_bounds input validation", {
  # Test with invalid inputs
  expect_error(tab_bounds(list()), "non-empty list")
  
  # Test with invalid probability column
  invalid_table <- data.frame(A = c(0, 1), wrong_col = c(0.5, 0.5))
  expect_error(
    tab_bounds(list(invalid_table), prob_col = "prob"),
    "not found"
  )
  
  # Test with negative probabilities
  negative_table <- data.frame(A = c(0, 1), prob = c(-0.1, 1.1))
  expect_error(
    tab_bounds(list(negative_table)),
    "Negative probabilities"
  )
})

test_that("tab_bounds verbose output works", {
  joint_ab <- data.frame(A = c(0, 1), B = c(0, 1), prob = c(0.3, 0.7))
  tables <- list(joint_ab)
  
  # Test that verbose mode doesn't throw errors
  expect_output(
    bounds <- tab_bounds(tables, verbose = TRUE),
    "Optimization complete"
  )
  
  expect_s3_class(bounds, "data.frame")
})