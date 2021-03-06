% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_Cross_Dep.R
\name{polyexp}
\alias{polyexp}
\title{Bivariate Polynomial Approximation}
\usage{
polyexp(xi, xj, k)
}
\arguments{
\item{xi}{first variable term}

\item{xj}{second variable term}

\item{k}{order of polynomial approximation}
}
\value{
vector of length \eqn{k(k+1)/2}
}
\description{
Take a k-order polynomial approximation of a bivariate \eqn{\phi(xi,xj)}. The expansion excludes 
terms that are solely in the first variable term because coefficients on them are not identified.
}
\examples{
polyexp(1.2,3.1,3)
polyexp(1,0,1) # an order one expansion useful for dummy variables
}
