\name{tb}
\alias{tb}
\title{Teitz-Bart algorithm applied to Spatial* and Spatial*DataFrame objects - report p-median set}
\usage{
  tb(swdf1, swdf2, p, metric, verbose = FALSE)
}
\arguments{
  \item{swdf1}{- first Spatial* or Spatial*DataFrame
  objects}

  \item{swdf2}{- second Spatial* or Spatial*DataFrame
  objects (if omitted, defaults to the same value as
  \code{swdf1})}

  \item{p}{- either a guess at the initial \eqn{p}-median
  set of a single integer indicating the size of the set
  (which is then chosen randomly)}

  \item{metric}{- the distance matrix (defaults to
  Euclidean computed via \code{euc.dists(swdf1,swdf2)} if
  not supplied)}

  \item{verbose}{- if TRUE print out each swap in the
  algorithm (default is FALSE)}
}
\value{
  Set of point indices for \eqn{p}-median (may be local
  optimum)
}
\description{
  Teitz-Bart algorithm applied to Spatial* and
  Spatial*DataFrame objects - report p-median set
}
\examples{
data(meuse)
coordinates(meuse) <- ~x+y
tb(meuse,p=5)
}

