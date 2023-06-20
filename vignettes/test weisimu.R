## install "weisimu" package
library(devtools)
install_github("Dorayaya/weisimu")
library(weisimu)

## use "testthat" to test package
library(testthat)

## Define test cases for simtrim function

test_simtrim <- function() {
  
  # Test case 1: Check if the function returns the correct list of results
  result <- simtrim(n = 100, shape = 2, scale = 1, p = 0.1, S = 100)
  expected_keys <- c("mu_0", "bias_mu_0", "var_mu_0", "MSE_mu_0",
                     "mu_p", "bias_mu_p", "var_mu_p", "MSE_mu_p")
  expect_equal(names(result), expected_keys)
  
  # Test case 2: Check if the bias of mu_0 is zero
  result_mu_0 <- simtrim(n = 200, shape = 2, scale = 1, p = 0.1, S = 1000)
  expect_true(abs(result_mu_0$bias_mu_0) == 0)
  
  # Test case 3: Check if the bias of mu_p is close to zero for large sample size
  result_large_n <- simtrim(n = 5000, shape = 2, scale = 1, p = 0.1, S = 1000)
  expect_true(abs(result_large_n$bias_mu_0) < 0.1)
  
  # Test case 4: check if mu_p and mu_0 are close when p is close to 0
  result <- simtrim(n = 5000, shape = 3, scale = 2, p = 0.0001, S = 1000)
  expect_true(all(abs(result$mu_p - result$mu_0) < 1e-3), label = "mu_p is very close to mu_0")
  
  # Test case 5: check error messages
  expect_error(simtrim(n = "invalid", shape = 2, scale = 1, p = 0.1, S = 100), "All parameters must be numeric")
  expect_error(simtrim(n = 100, shape = -2, scale = 1, p = 0.1, S = 100), "require shape > 0")
  expect_error(simtrim(n = 100, shape = 2, scale = -1, p = 0.1, S = 100), "require scale > 0")
  expect_error(simtrim(n = 100, shape = 2, scale = 1, p = 0.1, S = -100), "require S > 0")
  
  ## check some statistical properties
  # Test case 6: check if S->\infty, does simulation error reduce?
  sim_result <- simtrim(n = 1000, shape = 2, scale = 1, p = 0.1, S = 500)
  sim_result_large_S <- simtrim(n = 1000, shape = 2, scale = 1, p = 0.1, S = 10000)
  expect_true(sim_result_large_S$MSE_mu_p < sim_result$MSE_mu_p)
  
  # Test case 7: for given values of shape and scale, check if the mean/var/MSE of \hat mu_0 match the theoretical values
  sim_result_large_S <- simtrim(n = 2000, shape = 2, scale = 1, p = 0.1, S = 5000)
  shape <- 2
  scale <- 1
  n <- 2000
  true_mean <- scale * gamma(1 + 1/shape)
  true_variance <- (scale^2) * (gamma(1 + 2/shape) - (gamma(1 + 1/shape))^2) / n
  true_MSE <- 0^2 + true_variance
  expect_true(abs(sim_result_large_S$mu_0 - true_mean) < 1e-5)
  expect_true(abs(sim_result_large_S$var_mu_0 - true_variance) < 1e-5)
  expect_true(abs(sim_result_large_S$MSE_mu_0 - true_MSE) < 1e-5)

}

## Define test cases for simtrim_by function
test_simtrim_by <- function() {
  # Test case 1: Check if the function returns the expected output format
  result <- simtrim_by(n = c(100, 200), shape = 2, scale = 1, p = 0.1, S = 100, plot = FALSE)
  expect_is(result, "array")
  expect_true(all(names(result) %in% c("MSE_mu_0", "MSE_mu_p", "n")))
  
  # Test case 2: Check if the function produces a plot without any errors
  result_plot <- simtrim_by(n = c(100, 200), shape = 2, scale = 1, p = 0.1, S = 100, plot = TRUE)
  plot <- recordPlot()
  expect_false(is.null(plot))
  
  # Test case 3: Check if the function works when shape is a vector
  result <- simtrim_by(n = 100, shape = c(1, 2, 3), scale = 1, p = 0.1, S = 1000)
  expect_is(result, "array")
  expect_true(all(names(result) %in% c("MSE_mu_0", "MSE_mu_p", "n")))
  plot <- recordPlot()
  expect_false(is.null(plot))
  
  # Test case 4: Check if the function works when scale is a vector
  result <- simtrim_by(n = 100, shape = 2, scale = c(1, 2, 3), p = 0.05, S = 1000)
  expect_is(result, "array")
  expect_true(all(names(result) %in% c("MSE_mu_0", "MSE_mu_p", "n")))
  plot <- recordPlot()
  expect_false(is.null(plot))
  
  # Test case 6: Check if the function works when trimming proportion is a vector
  result <- simtrim_by(n = 100, shape = 2, scale = 1, p = c(0.01, 0.05, 0.1), S = 1000)
  expect_is(result, "array")
  expect_true(all(names(result) %in% c("MSE_mu_0", "MSE_mu_p", "n")))
  plot <- recordPlot()
  expect_false(is.null(plot))
  
  # Test case 7: Test error messages
  expect_error(simtrim_by(n = "invalid", shape = 2, scale = 1, p = 0.1, S = 100), "n, shape, scale, p, S must be numeric")
  expect_error(simtrim_by(n = c(50, 100), shape = -2, scale = 1, p = 0.1, S = 100), "require shape > 0")
  expect_error(simtrim_by(n = c(10,100), shape = c(1,2), scale = -1, p = 0.1, S = 100), "only one of \\(n, shape, scale, p\\) may be vector of length>1")
  expect_error(simtrim_by(n = 100, shape = 2, scale = 1, p = 0.1, S = 100), "one of \\(n, shape, scale, p\\) must be vector of length>1")
  expect_error(simtrim_by(n = c(0,1), shape = 2, scale = -1, p = 0.1, S = 100), "require n >\\= 2")
  expect_error(simtrim_by(n = c(100,200), shape = 1, scale = -1, p = 0.1, S = 100), "require scale > 0")
  expect_error(simtrim_by(n = c(100,200), shape = 1, scale = 1, p = 0.1, S = -100), "require S > 0")
  expect_error(simtrim_by(n = c(100,200), shape = 1, scale = 1, p = -0.1, S = 100), "'probs' outside \\[0,1\\]")
}


# Run the tests
test_that("simtrim function tests", {test_simtrim()})

test_that("simtrim_by function tests", {test_simtrim_by()})
