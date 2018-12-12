library(testthat)
library(ars)
source('ars.R')

test_that("test chisq(3)", {
  out <- ars(function(x){dchisq(x, df = 3)}, n=5000, lb=0)
  reject = TRUE
  if(ks.test(out, rchisq(length(out), 3))$p.value < 0.05){
    reject = FALSE
  }
  expect(reject==TRUE, failure_message = "ERROR: output does not follow the expected distribution. Check distribution parameters, lb, and ub.")
})

test_that("test chisq(6)", {
  out <- ars(function(x){dchisq(x, df = 6)}, n=5000, lb=0)
  reject = TRUE
  if(ks.test(out, rchisq(length(out), 6))$p.value < 0.05){
    reject = FALSE
  }
  expect(reject==TRUE, failure_message = "ERROR: output does not follow the expected distribution. Check distribution parameters, lb, and ub.")
})

test_that("test beta(2, 2)", {
  out <- ars(function(x){dbeta(x, shape1 = 2, shape2 = 2)}, n=5000, lb=0, ub=1)
  reject = TRUE
  if(ks.test(out, rbeta(length(out), 2, 2))$p.value < 0.05){
    reject = FALSE
  }
  expect(reject==TRUE, failure_message = "ERROR: output does not follow the expected distribution. Check distribution parameters, lb, and ub.")
})

test_that("test chisq(6) data is not from chisq(5)", {
  out <- ars(function(x){dchisq(x, df = 6)}, n=5000, lb=0)
  reject = TRUE
  if(ks.test(out, rchisq(length(out), 5))$p.value < 0.05){
    reject = FALSE
  }
  expect(reject==FALSE, failure_message = "ERROR: output should not follow the expected distribution. Check distribution parameters, lb, and ub.")
})


test_that("test non log-concave chisq(1)", {
  expect_error(ars(function(x){dchisq(x, df = 1)}, n=5000, lb=0, ub=1))
})

test_that("test non log-concave gamma(0.5, 1)", {
  expect_error(ars(function(x){dgamma(x, 0.5, 1)}, n=5000, lb=0, ub=10))
})

test_that("test non log-concave x^2", {
  expect_error(ars(function(x){x^2}, n=5000))
})


