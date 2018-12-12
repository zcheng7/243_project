\name{ars}
\alias{ars}
\title{Adapative Rejection Sampling}
\usage{
ars(G, n, lb, ub)
}
\arguments{
\item{G}{function g(x) that is proportional to the target density for sampling. This function must be log-concave}

\item{n}{The number of points to sample}

\item{lb}{The lower bound for the function g(x). By default, the function assumes no lower bound (lb = -Inf).}

\item{ub}{The upper bound for the function g(x). By default, the function assumes no upper bound (ub = Inf)}

}
\value{
\code{ars} returns a list containing the following components:
\item{new_sample}{A vector of length n conatining the sampled points}
}
\description{
Adaptive rejection sampling algorithm for functions that are log-concave. Function will return a vector of points of desired length, using rejection sampling to approximate the target distribution. Algorithm is based on Gilks & Wild, 1992.
}
\examples{
## Example 1: Standard Normal
g <- function(x) dnorm(x)
samp <- ars(g, n = 500)
density(samp)


## Example 2: Chi-Squared Distribution with Lower Bound
g <- function(x) dchisq(x, df = 5)
samp <- ars(g, n = 500, lb = 0)
density(samp)
}
\author{
Package developed by Zichen Chen, Isabel Zhou, and Aidan Mullan based on Gilks & Wild, 1992.
}
\references{
Gilks, W.R., P. Wild. (1992) Adaptive Rejection Sampling for Gibbs Sampling, Applied Statistics 41:337–348.
}
