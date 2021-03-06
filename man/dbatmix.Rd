% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BatMixtures.R
\name{dbatmix}
\alias{dbatmix}
\alias{dbatmix_pmat}
\alias{rinvbatmix}
\title{The Batschelet-type Mixture distribution}
\usage{
dbatmix(x, dbat_fun = dinvbat, mus = c(-pi/2, 0, pi/2), kps = c(8, 8,
  8), lams = c(-0.5, 0, 0.5), alphs = c(0.3, 0.4, 0.3), log = FALSE)

dbatmix_pmat(x, dbat_fun = dinvbat, pmat = cbind(mu = c(-pi/2, 0,
  pi/2), kp = c(8, 8, 8), lam = c(-0.5, 0, 0.5), alph = c(0.3, 0.4, 0.3)),
  log = FALSE)

rinvbatmix(n, mus = c(-pi/2, 0, pi/2), kps = c(8, 8, 8),
  lams = c(-0.5, 0, 0.5), alphs = c(0.3, 0.4, 0.3))
}
\arguments{
\item{x}{An angle in radians.}

\item{dbat_fun}{A function that provides the pdf function of the desired}

\item{mus}{A vector of component mean directions.}

\item{kps}{A vector of nonnegative concentration parameters.}

\item{lams}{A vector of shape parameters.}

\item{alphs}{A vector of component weights. Batschelet density. By default,
inverse Batschelet.}

\item{log}{Logical; indicates if the logarithm of the density must be return.}

\item{pmat}{A \code{k*4} matrix of parameters for \code{dbatmix_pmat}.}

\item{n}{The number of data points to sample.}
}
\value{
For \code{dbatmix} and \code{dbatmix_pmat}, vector of probabilities.
  For \code{rinvbatmix}, a vector of random variates from the mixture of the
  inverse Batschelet distribution.
}
\description{
The Batschelet-type Mixture distribution
}
\section{Functions}{
\itemize{
\item \code{dbatmix_pmat}: A version that takes a parameter matrix as input.

\item \code{rinvbatmix}: Random variate generation for Inverse Batschelet mixtures.
}}

\examples{
x <- rinvbat(10)
dbatmix(x)

curve(dbatmix(x), -pi, pi)
curve(dbatmix(x, dbat_fun = dpowbat), -pi, pi)

}
