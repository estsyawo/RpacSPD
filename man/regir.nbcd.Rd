% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions_Model_Criteria.R
\name{regir.nbcd}
\alias{regir.nbcd}
\title{Regression - nbcd class}
\usage{
\method{regir}{nbcd}(object, ...)
}
\arguments{
\item{object}{a list of Y - outcome variable and Xmat - design matrix of class "nbcd"}

\item{...}{additional parameters to be passed to \code{\link[MASS]{glm.nb}}}
}
\value{
fitted model object
}
\description{
A negative binomial regression implementation for the "nbcd" class. It uses \code{\link[MASS]{glm.nb}}
}
