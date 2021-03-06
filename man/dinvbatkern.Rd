% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InvBatschelet.R
\name{dinvbatkern}
\alias{dinvbatkern}
\title{Kernel of the von-Mises based symmetric inverse Batschelet distribution}
\usage{
dinvbatkern(x, mu = 0, kp = 1, lam = 0, log = FALSE)
}
\arguments{
\item{x}{An angle in radians.}

\item{mu}{A mean direction, in radians.}

\item{kp}{Numeric, \eqn{> 0,}the concentration parameter.}

\item{lam}{The shape parameter (peakedness), -1 < \code{lam} < 1.}

\item{log}{Logical; whether to return the log of the probability or not.}
}
\value{
The unnormalized density value of the inverse Batschelet distribution.
}
\description{
Kernel of the von-Mises based symmetric inverse Batschelet distribution
}
