# Tests for tab_bounds against analytical 2x2 bounds

library(mrpbounds)

# Analytical bounds for 2x2 table cells given marginals
bounds_2x2_cells <- function(row1_margin, col1_margin) {
  list(
    p00 = c(
      lower = max(0, row1_margin + col1_margin - 1),
      upper = min(row1_margin, col1_margin)
    ),
    p01 = c(
      lower = max(0, row1_margin - col1_margin),
      upper = min(row1_margin, 1 - col1_margin)
    ),
    p10 = c(
      lower = max(0, col1_margin - row1_margin),
      upper = min(1 - row1_margin, col1_margin)
    ),
    p11 = c(
      lower = max(0, 1 - row1_margin - col1_margin),
      upper = min(1 - row1_margin, 1 - col1_margin)
    )
  )
}

create_2x2_margin_tables <- function(row1_margin, col1_margin) {
  list(
    data.frame(row = c(0, 1), prob = c(row1_margin, 1 - row1_margin)),
    data.frame(col = c(0, 1), prob = c(col1_margin, 1 - col1_margin))
  )
}

test_that("tab_bounds matches exact bounds_2x2_cells for 2x2 tables", {
  tol <- 1e-6

  # Symmetric margins
  exact <- bounds_2x2_cells(0.5, 0.5)
  num <- tab_bounds(create_2x2_margin_tables(0.5, 0.5), verbose = FALSE)
  num <- num[order(num$row, num$col), ]
  expect_equal(num$min_bound[1], as.numeric(exact$p00["lower"]), tolerance = tol)
  expect_equal(num$max_bound[1], as.numeric(exact$p00["upper"]), tolerance = tol)
  expect_equal(num$min_bound[2], as.numeric(exact$p01["lower"]), tolerance = tol)
  expect_equal(num$max_bound[2], as.numeric(exact$p01["upper"]), tolerance = tol)
  expect_equal(num$min_bound[3], as.numeric(exact$p10["lower"]), tolerance = tol)
  expect_equal(num$max_bound[3], as.numeric(exact$p10["upper"]), tolerance = tol)
  expect_equal(num$min_bound[4], as.numeric(exact$p11["lower"]), tolerance = tol)
  expect_equal(num$max_bound[4], as.numeric(exact$p11["upper"]), tolerance = tol)

  # Asymmetric margins
  exact <- bounds_2x2_cells(0.3, 0.4)
  num <- tab_bounds(create_2x2_margin_tables(0.3, 0.4), verbose = FALSE)
  num <- num[order(num$row, num$col), ]
  expect_equal(num$min_bound[1], as.numeric(exact$p00["lower"]), tolerance = tol)
  expect_equal(num$max_bound[1], as.numeric(exact$p00["upper"]), tolerance = tol)
  expect_equal(num$min_bound[2], as.numeric(exact$p01["lower"]), tolerance = tol)
  expect_equal(num$max_bound[2], as.numeric(exact$p01["upper"]), tolerance = tol)
  expect_equal(num$min_bound[3], as.numeric(exact$p10["lower"]), tolerance = tol)
  expect_equal(num$max_bound[3], as.numeric(exact$p10["upper"]), tolerance = tol)
  expect_equal(num$min_bound[4], as.numeric(exact$p11["lower"]), tolerance = tol)
  expect_equal(num$max_bound[4], as.numeric(exact$p11["upper"]), tolerance = tol)

  # Extreme margins
  exact <- bounds_2x2_cells(0.1, 0.2)
  num <- tab_bounds(create_2x2_margin_tables(0.1, 0.2), verbose = FALSE)
  num <- num[order(num$row, num$col), ]
  expect_equal(num$min_bound[1], as.numeric(exact$p00["lower"]), tolerance = tol)
  expect_equal(num$max_bound[1], as.numeric(exact$p00["upper"]), tolerance = tol)
  expect_equal(num$min_bound[2], as.numeric(exact$p01["lower"]), tolerance = tol)
  expect_equal(num$max_bound[2], as.numeric(exact$p01["upper"]), tolerance = tol)
  expect_equal(num$min_bound[3], as.numeric(exact$p10["lower"]), tolerance = tol)
  expect_equal(num$max_bound[3], as.numeric(exact$p10["upper"]), tolerance = tol)
  expect_equal(num$min_bound[4], as.numeric(exact$p11["lower"]), tolerance = tol)
  expect_equal(num$max_bound[4], as.numeric(exact$p11["upper"]), tolerance = tol)

  # Edge case with some bounds exactly 0
  exact <- bounds_2x2_cells(0.2, 0.1)
  num <- tab_bounds(create_2x2_margin_tables(0.2, 0.1), verbose = FALSE)
  num <- num[order(num$row, num$col), ]
  expect_equal(num$min_bound[1], as.numeric(exact$p00["lower"]), tolerance = tol)
  expect_equal(num$max_bound[1], as.numeric(exact$p00["upper"]), tolerance = tol)
  expect_equal(num$min_bound[2], as.numeric(exact$p01["lower"]), tolerance = tol)
  expect_equal(num$max_bound[2], as.numeric(exact$p01["upper"]), tolerance = tol)
  expect_equal(num$min_bound[3], as.numeric(exact$p10["lower"]), tolerance = tol)
  expect_equal(num$max_bound[3], as.numeric(exact$p10["upper"]), tolerance = tol)
  expect_equal(num$min_bound[4], as.numeric(exact$p11["lower"]), tolerance = tol)
  expect_equal(num$max_bound[4], as.numeric(exact$p11["upper"]), tolerance = tol)
})
