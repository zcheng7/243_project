library(testthat)
library(ars)



test_that("Test empirical distributions", {
  
  #1
  print("Now testing ars chisq(3) and rchisq(3): ")
  out <- ars(function(x){dchisq(x, df = 3)}, n=5000, lb=0)
  pv <- ks.test(out, rchisq(length(out), 3))$p.value
  print("KS test estimate p-value is: ")
  print(pv)
  cat("\n\n")
  if(pv > 0.05){
    print("TEST PASSED for Chisq(3).")
  } else {
    print("!! TEST FAILED: output does not follow the expected distribution. Check distribution parameters, lb, and ub.")
  }
  
  #2
  print("Now testing ars chisq(6) and rchisq(6): ")
  out <- ars(function(x){dchisq(x, df = 6)}, n=5000, lb=0)
  pv <- ks.test(out, rchisq(length(out), 6))$p.value
  print("KS test estimate p-value is: ")
  print(pv)
  cat("\n\n")
  if(pv > 0.05){
    print("TEST PASSED for Chisq(6).")
  } else {
    print("!! TEST FAILED: output does not follow the expected distribution. Check distribution parameters, lb, and ub.")
  }
  
  #3
  print("Now testing ars exp(1) and rexp(1): ")
  out <- ars(function(x){dexp(x, 1)}, n=5000)
  pv <- ks.test(out, rexp(length(out), 1))$p.value
  print("KS test estimate p-value is: ")
  print(pv)
  cat("\n\n")
  if(pv > 0.05){
    print("** TEST PASSED for Chisq(3).")
  } else {
    print("!! TEST FAILED: output does not follow the expected distribution. Check distribution parameters, lb, and ub.")
  }
  
  #4
  print("Now testing ars beta(2, 2) and rbeta(2, 2): ")
  out <- ars(function(x){dbeta(x, shape1 = 2, shape2 = 2)}, n=5000, lb=0, ub=1)
  pv <- ks.test(out, rbeta(length(out), 2, 2))$p.value
  print("KS test estimate p-value is: ")
  print(pv)
  cat("\n\n")
  if(pv > 0.05){
    print("TEST PASSED for Beta(2, 2).")
  } else {
    print("!! TEST FAILED: output does not follow the expected distribution. Check distribution parameters, lb, and ub.")
  }

  #5
  print("Now testing ars chisq(5) and rchisq(6): ")
  out <- ars(function(x){dchisq(x, df = 5)}, n=5000, lb=0)
  pv <- ks.test(out, rchisq(length(out), 6))$p.value
  print("KS test estimate p-value is: ")
  print(pv)
  cat("\n\n")
  if(pv < 0.05){
    print("** TEST PASSED for chisq(5) and rchisq(6).")
  } else {
    print("!! TEST FAILED: output should not follow a different distribution. Check distribution parameters, lb, and ub.")
  }
  
  #6
  print("Now testing non-log-concave chisq(1): ")
  expect_error(ars(function(x){dchisq(x, df = 1)}, n=5000))
  print("** TEST PASSED for non-log-concave Chisq(1).")
  
  #7
  print("Now testing non-log-concave gamma(0.5, 1): ")
  expect_error(ars(function(x){dgamma(x, 0.5, 1)}, n=5000))
  print("** TEST PASSED for non-log-concave Gamma(0.5, 1).")
  
  #8
  print("Now testing non-log-concave x^2: ")
  expect_error(ars(function(x){x^2}, n=5000))
  print("** TEST PASSED for non-log-concave x^2.")

})

