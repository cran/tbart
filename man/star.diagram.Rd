\name{star.diagram}
\alias{star.diagram}
\title{Creates the lines for a 'star diagram'}
\usage{
  star.diagram(swdf1, swdf2, alloc)
}
\arguments{
  \item{swdf1}{- first Spatial* or Spatial*DataFrame
  objects}

  \item{swdf2}{- second Spatial* or Spatial*DataFrame
  objects (if omitted, defaults to the same value as
  \code{swdf1})}

  \item{alloc}{- a list saying which coordinate in swdf2 is
  allocated to each point in swdf1}
}
\description{
  Creates the lines for a 'star diagram'
}
\examples{
data(meuse)
coordinates(meuse) <- ~x+y
allocations.list <- allocate(meuse,p=5)
star.lines <- star.diagram(meuse,alloc=allocations.list)
plot(star.lines)
}

