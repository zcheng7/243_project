library(assertthat)
library(testthat)

## construct zk using formula from the paper
## input: point x, logfunction h, lowerbound, upper bound.
## output: z_k
generate_zk <- function(x,h,lb,ub){  
  
  ## tests for inputs
  assert_that(is.numeric(x), see_if(length(x)>0), msg = "ERROR: x is not a numeric in function generate_zk")
  assert_that(is.function(h), msg = "ERROR: h is not a function in function generate_zk")
  assert_that(is.numeric(lb), msg = "ERROR: lb is not numeric in function generate_zk")
  assert_that(is.numeric(ub), msg = "ERROR: ub is not numeric in function generate_zk")
  assert_that(see_if(lb<ub), msg = "ERROR: lb should be smaller than ub in function generate_zk")
  
  k = length(x)
  x_j1 <- c(1e-8,x)
  x_j <- c(x,1e-8)
  zk <- (h(x_j1)-h(x_j)-x_j1*derive(x_j1,h)+x_j*derive(x_j,h))/(derive(x_j,h)-derive(x_j1,h))
  tmp <- zk
  tmp[1] <- lb
  tmp[length(x)+1] <- ub
  return(tmp)
}

## function find gradient at given point
## input: point x, function
## output: gradient of function at x.
derive <- function(x, f){
  
  # tests for inputs
  assert_that(is.numeric(x), see_if(length(x)>0), msg = "ERROR: x is not a numeric in function derive")
  assert_that(is.function(f), msg = "ERROR: f is not a function in function derive")
  
  return((f(x + 1e-8) - f(x - 1e-8))/2e-8)
}

# When no bounds are given, find starting values on each side of the mode.
# input: logfunction h, lowerbound, upperbound, mode.
# output: lowerbound, upperbound.
find_two_points <- function(h, lb, ub){
  mode <- optim(par=0, fn = h, method = "L-BFGS-B", lower = lb, upper = ub, control=list(fnscale=-1))$par
  ## tests for inputs
  assert_that(is.function(h), msg = "ERROR: h is not a function in function find_two_points")
  assert_that(is.numeric(lb), msg = "ERROR: lb is not numeric in function find_two_points")
  assert_that(is.numeric(ub), msg = "ERROR: ub is not numeric in function find_two_points")
  assert_that(see_if(lb<ub), msg = "ERROR: lb should be smaller than ub in function find_two_points")
  
  
  ## If unbounded below, step from mode until derivative = 0, or 100
  lb_iterate <- ub_iterate <- 1
  if(lb == -Inf){
    lb_test <- mode - 10 
    lb_slope <- derive(lb_test, h)
    #print(lb_slope)
    while(lb_test > -Inf && abs(lb_slope) >= 1e-5 && lb_iterate < 100){
      lb_test <- lb_test - 10
      lb_slope <- derive(lb_test, h)
      lb_iterate <- lb_iterate + 1
      #print(lb_iterate)
    }
    lb <- lb_test
  }
  ## If unbounded above, step from mode until derivative = 0 or 100 steps
  if(ub == Inf){
    ub_test <- mode + 10 
    ub_slope <- derive(ub_test, h)
    #print(ub_slope)
    while(ub_test < Inf && abs(ub_slope) >= 1e-5 && ub_iterate < 100){
      #print(ub_iterate)
      ub_test <- ub_test + 10
      ub_slope <- derive(ub_test, h)
      ub_iterate <- ub_iterate + 1
    }
    ub <- ub_test
    #print(ub)
  }
  ## Provide warning if too many iterations
  if(lb_iterate == 100 | ub_iterate == 100){
    print ("Starting values could not be found after 100 iterations, 100th value being used.")
  }
  ## Stop function if function maximum is not between bounds
  lb_slope <- derive(lb, h)
  ub_slope <- derive(ub, h)
  if(lb_slope < 0 | ub_slope > 0){
    stop("Starting values could not be identified. Please make check bounds and estimate of function mode", call. = FALSE)
  }
  return(c(lb, ub))
}


## Using function bounds, construct starting basis x-values to create z-values and envelope
## input: logfunctionh, lowerbound, upperbound
## output: three x values for initiation of the algorithm
init_basis_x <- function(h, lb, ub){
  
  ## tests for inputs
  assert_that(is.function(h), msg = "ERROR: h is not a function in function init_basis_x")
  assert_that(is.numeric(lb), msg = "ERROR: lb is not numeric in function init_basis_x")
  assert_that(is.numeric(ub), msg = "ERROR: ub is not numeric in function init_basis_x")
  assert_that(see_if(lb<ub), msg = "ERROR: lb should be smaller than ub in function init_basis_x")
  
  delta <- (ub - lb)/500
  max <- optimize(f = h, interval = c(lb, ub), lower = lb, upper = ub, maximum = TRUE)$maximum
  #taking care of exp case
  if (abs(max - ub) < .0001) {
    right_point <- max
    mid <- max - .5*delta
    left_point <- (max - delta)
    X_init <- c(left_point,mid,right_point)
  } else if (abs(max - lb) < .0001) {
    right_point <- (max + delta)
    mid <- max + .5*delta
    left_point <- max
    X_init <- c(left_point,right_point)
  } else {
    right_point <- (max + delta)
    left_point <- (max - delta)
    X_init <- c(left_point,max,right_point)
  }
  return(X_init)
}


## find lower hull for certain points
## input: point x, h_k, gradient at h_k, x_k
## output: lower hull for point x
lower_function <- function(x,h_k,dh_k,x_k){
  
  ## tests for inputs
  assert_that(is.numeric(x), msg = "ERROR: x is not numeric in function lower_function")
  assert_that(is.numeric(h_k), msg = "ERROR: h_k is not numeric in function lower_function")
  assert_that(is.numeric(dh_k), msg = "ERROR: dh_k is not numeric in function lower_function")
  assert_that(is.numeric(x_k), msg = "ERROR: x_k is not numeric in function lower_function")  
  
  if(x < min(x_k) | x > max(x_k)){
    return(-Inf)
  }else{
    i <- max(which(x >= x_k))
    lower_func <- ((x_k[i + 1] - x) * x_k[i] + (x - x_k[i]) * h_k[i + 1])/(x_k[i + 1] - x_k[i])
    return(lower_func)
    
  }
}

## find upper hull for certain point. And exponential of it.
## input: point x, h_k, gradient at h_k, x_k, z_k
## output: upper hull at point x for first function, and exponential upper hull at point x for second function.
upper_function <- function(x,h_k,dh_k,x_k,z_k){
  
  ## tests for inputs
  assert_that(is.numeric(x), msg = "ERROR: x is not numeric in function upper_function")
  assert_that(is.numeric(h_k), msg = "ERROR: h_k is not numeric in function upper_function")
  assert_that(is.numeric(dh_k), msg = "ERROR: dh_k is not numeric in function upper_function")
  assert_that(is.numeric(x_k), msg = "ERROR: x_k is not numeric in function upper_function") 
  assert_that(is.numeric(z_k), msg = "ERROR: z_k is not numeric in function upper_function") 
  
  i <- min(which(x < z_k)-1)
  upper_func <- h_k[i] + (x - x_k[i]) * dh_k[i]
  return(upper_func)
}

exp_upper <- function(x,h_k,dh_k,x_k,z_k){
  
  ## tests for inputs
  assert_that(is.numeric(x), msg = "ERROR: x is not numeric in function exp_upper")
  assert_that(is.numeric(h_k), msg = "ERROR: h_k is not numeric in function exp_upper")
  assert_that(is.numeric(dh_k), msg = "ERROR: dh_k is not numeric in function exp_upper")
  assert_that(is.numeric(x_k), msg = "ERROR: x_k is not numeric in function exp_upper") 
  assert_that(is.numeric(z_k), msg = "ERROR: z_k is not numeric in function exp_upper") 
  
  exp_up <- exp(upper_function(x,h_k,dh_k,x_k,z_k))
  return(exp_up)
}

## drawing sample from the envelope using inverse cdf
## input: point u, cumulative density area, h_k, x_k, gradient of h_k, z_k, probability density area
## output: sample points drawn from envelop area
generate_sample <- function(u,cum_env,h_k,x_k,dh_k,z_k,area){
  
  ## tests for inputs
  assert_that(is.numeric(u), msg = "ERROR: u is not numeric in function generate_sample")
  assert_that(is.numeric(cum_env), msg = "ERROR: cum_env is not numeric in function generate_sample")
  assert_that(is.numeric(h_k), msg = "ERROR: h_k is not numeric in function generate_sample") 
  assert_that(is.numeric(x_k), msg = "ERROR: x_k is not numeric in function generate_sample") 
  assert_that(is.numeric(dh_k), msg = "ERROR: dh_k is not numeric in function generate_sample")
  assert_that(is.numeric(z_k), msg = "ERROR: z_k is not numeric in function generate_sample") 
  assert_that(is.numeric(x_k), msg = "ERROR: x_k is not numeric in function generate_sample")
  
  j <- max(which(u > cum_env))
  if(dh_k[j] == 0){
    x <- runif(1, z_k[j],z_k[j+1])
    return(x)
  }else{
    
    # Sample from uniform random
    w = runif(1)
    
    # Rescale sample value w to area of the selected segment, since area under segment is not equal to 1
    w_sc = w*(1/area[j])*exp(h_k[j] - x_k[j]*dh_k[j])*(exp(dh_k[j]*z_k[j+1]) -exp(dh_k[j]*z_k[j]))
    
    # Use inverse CDF of selected segment to generate a sample
    x = (1/dh_k[j])*log(w_sc*area[j]/(exp(h_k[j] - x_k[j]*dh_k[j])) + exp(z_k[j]*dh_k[j]))
  }
  return(x)
}


## adaptive rejection test
## input: sample point x, h_k, gradient of h_k, x_k, z_k, logfunction h
## output: indicator of which of the three scenario does sample point x lies in
rej_test <- function(x, h_k, dh_k, x_k, z_k, h){
  
  ## tests for inputs
  assert_that(is.numeric(x), msg = "ERROR: x is not numeric in function rej_test")
  assert_that(is.numeric(h_k), msg = "ERROR: h_k is not numeric in function rej_test")
  assert_that(is.numeric(dh_k), msg = "ERROR: dh_k is not numeric in function rej_test")
  assert_that(is.numeric(x_k), msg = "ERROR: x_k is not numeric in function rej_test") 
  assert_that(is.numeric(z_k), msg = "ERROR: z_k is not numeric in function rej_test")
  assert_that(is.function(h), msg = "ERROR: h is not a function in function rej_test")
  
  # Generate random seed
  w = runif(1)
  
  # squeeze and reject tests indicator for adding point in boolean form
  squeeze = FALSE
  accept = FALSE
  add = FALSE
  
  # get rejection point for squeeze and accept test
  lower_test = exp(lower_function(x,h_k,dh_k,x_k) - upper_function(x,h_k,dh_k,x_k,z_k))
  
  if(w <= lower_test){
    
    squeeze = TRUE
    accept = TRUE
    
  }else{
    function_test = exp(h(x) - upper_function(x,h_k,dh_k,x_k,z_k))
    if (w <= function_test){
      
      squeeze = FALSE
      accept = TRUE
      
    }else{
      accept = FALSE
    }
  } 
  # Determine whether to add point to abscissae
  if(squeeze * accept == FALSE) add = TRUE
  
  # Return boolean indicating whether to accept candidate sample point
  return(list(rej = squeeze + accept,add = add))
}

#########################################################################
#' Adaptive Rejection Sampling (ars) Function
#'
#' Adaptive rejection sampling function that creates a sample
#' from a density function based on the Gilks(1992) paper.
#'
#' @param g probability density function, as a function of x (e.g. g <- function(x) dnorm(x,0,1))
#' @param n number of values the final sample should contain
#' @param lb lower bound to evaluate the density on
#' @param ub upper bound to evaluate the density on
#' @return sample of with specified (n) elements
#' @export

##########################################################################

ars <- function(g, n, lb = -Inf, ub = Inf){
  options(warn=-1)
  ## tests for inputs
  assert_that(is.function(g), msg = "ERROR: g is not a function in function ars. Try different lb and/or ub.")
  assert_that(is.numeric(n), msg = "ERROR: n is not numeric in function ars")
  
  #log of the original function
    h <- function(x){
    return(log(g(x)))
  }

  
  #find starting x_k
  if(lb == -Inf | ub == Inf){
    init_bound <- find_two_points(g, lb, ub)
  } else{
    init_bound <- c(lb, ub)
  }
  
  
  
  
  
  
    # test for log concavity 
    xt <- seq(init_bound[1], init_bound[2],length.out = 200)
    dhk.test <- derive(xt, h)
    
    #print(dhk.test)
    
    concavity = TRUE
    iter = 1
    while(concavity & iter < length(dhk.test)-1){
      assert_that(see_if(is.na(dhk.test[iter])==FALSE), msg = "ERROR: Change boundaries and try again.")
      if(dhk.test[iter] >= dhk.test[iter+1]) {iter <- iter+1}
      else {concavity = FALSE}
    }
    assert_that(see_if(concavity==TRUE), msg = "ERROR: g is not log-concave in function ars. Check function g and boundaries.")
  
  
  
  
  
  
  x_k <- init_basis_x(h, init_bound[1], init_bound[2])
  
  #initialize output variable
  new_sample <- NULL
  
  #iterate until we have enough points
  while(length(new_sample) < n){
    
    #intersection points
    z_k <- generate_zk(x_k,h,init_bound[1],init_bound[2])
    indexna <- which(is.na(z_k))
    x_k[indexna] <- NaN
    x_k <- na.omit(x_k)
    z_k <- na.omit(z_k)
    #h_k & h_k'
    h_k <- h(x_k)
    dh_k <- derive(x_k, h)
    #print(dh_k)
    
    #cumulative envelop 
    #Calculate areas under exponential upper bound function for normalization purposes
    
    area <- unlist(sapply(2:length(z_k),function(i){integrate(exp_upper,z_k[i-1],z_k[i],h_k,dh_k,x_k,z_k)})[1,])
    cum <- sum(area)
    #Normalize
    envelop <- area/cum
    cum_env <- cumsum(envelop)
    cum_env <- c(0,cum_env)
    
    #Sampling: Generate seeds for Inverse CDF method
    seeds <- runif(100)
    x_sample <- sapply(seeds, generate_sample,
                       cum_env = cum_env,
                       h_k = h_k,
                       x_k = x_k,
                       dh_k = dh_k,
                       z_k = z_k,
                       area = area)
    
    #Rejection testing
    test_result <- sapply(x_sample, rej_test,
                          h_k = h_k,
                          x_k = x_k,
                          dh_k = dh_k,
                          z_k = z_k,
                          h = h)
    keep_sample <- test_result[1,]
    add_to_x_k <- test_result[2,]
    
    # update accpeted points to sample
    x_keep <- x_sample[keep_sample > 0]
    new_sample <- c(x_keep, new_sample)
    
    # update x_k
    x_k_updated = x_sample[add_to_x_k > 0]
    x_k = sort(c(x_k, x_k_updated))
    
    
  }
  
  # test for output
  assert_that(see_if(min(new_sample)>=lb), msg = "ERROR: output sample is smaller than lower bound. Check boundaries.")
  assert_that(see_if(max(new_sample)<=ub), msg = "ERROR: output sample is larger than upper bound. Check boundaries.")
  
  return(new_sample)
  options(warn=0)
}




G <- function(x) dchisq(x, 3)
x <- ars(G,n = 1000)
hist(x)
# h0: our output and a random sample from the known distribution are drawn from the same distribution
# if p value from ks test > 0.05, then we fail to reject the h0
reject = TRUE
if(ks.test(x, rchisq(1000, 3))$p.value > 0.05){
  reject = FALSE
}
