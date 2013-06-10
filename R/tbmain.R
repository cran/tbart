.complement <- function(ivec, imax) {
  result <- rep(TRUE,imax)
  result[ivec] <- FALSE
  return(which(result))
}


.tb <- function(d,guess,verbose=FALSE) {
  config <- guess
  n <- ncol(d)
  repeat {
    old.config <- config
    config <- .bestswap(d,config,.complement(config,n))
    if (verbose) {
      cat("Configuration: ",config)
      score <- .dtotal(d,config)
      cat("  Score:",score,"\n")
    }
    if (all(old.config==config)) break 
  }
  return(config)
}

##' Teitz-Bart algorithm applied to a 'raw' distance matrix 
##'  
##' @param  d - A distance matrix (not necessarily Euclidean)
##' @param  guess - a guess at the set of p points constituting the \eqn{p}-median
##' @param  verbose - if TRUE print out each swap in the algorithm (default is FALSE)
##' @return Set of point indices for \eqn{p}-median (may be local optimum)
##' 
##' @examples 
##' x1 <- rnorm(100)
##' y1 <- rnorm(100)
##' d <- as.matrix(dist(cbind(x1,y1)))
##' tb.raw(d,c(1,2))
##' @export
##' 


tb.raw <- .tb

##' Euclidean distances from a Spatial* or Spatial*DataFrame object
##'  
##' @param  swdf1 - First Spatial*DataFrame object
##' @param  swdf2 - Second Spatial*DataFrame object (if omitted,  defaults to the same value as \code{swdf1})
##' @param  scale - allows re-scaling eg: value of 1000 means distances in km if coordinates of \code{swdf1}/\code{swdf2} in meters.
##' @return Distance matrix (if \code{swdf1} or \code{swdf2} not SpatialPoints*, distances are based on points obtained from \code{coordinates} function)
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' euc.dists(meuse,scale=1000)
##'
##' @export
##' 


euc.dists <- function(swdf1,swdf2,scale) {
  if (missing(swdf2)) swdf2 <- swdf1
  xy1 <- coordinates(swdf1)
  xy2 <- coordinates(swdf2)
  if (!missing(scale)) {
    xy1 <- xy1/scale
    xy2 <- xy2/scale
  }
  return(.dmat(xy1[,1],xy2[,1],xy1[,2],xy2[,2])) 
}


##' Teitz-Bart algorithm applied to Spatial* and Spatial*DataFrame objects - report p-median set
##'  
##' @param  swdf1 - first Spatial* or Spatial*DataFrame objects
##' @param  swdf2 - second Spatial* or Spatial*DataFrame objects (if omitted,  defaults to the same value as \code{swdf1})
##' @param  p - either a guess at the initial \eqn{p}-median set of a single integer indicating the size of the set (which is then chosen randomly)
##' @param  metric - the distance matrix (defaults to Euclidean computed via \code{euc.dists(swdf1,swdf2)} if not supplied)
##' @param  verbose - if TRUE print out each swap in the algorithm (default is FALSE)
##' @return Set of point indices for \eqn{p}-median (may be local optimum)
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' tb(meuse,p=5)
##' @export
##' 

tb <- function(swdf1,swdf2,p,metric,verbose=FALSE) {
  if (missing(swdf2))  swdf2 <- swdf1
  n.choices <- nrow(coordinates(swdf2)) 
  if (length(p) == 1) p <- sample(n.choices,p)
  if (missing(metric)) metric <- euc.dists(swdf1,swdf2)
  result <- .tb(metric,p,verbose)
  return(result)
}

##' Teitz-Bart algorithm applied to Spatial* and Spatial*DataFrame objects - report allocations
##'  
##' @param  swdf1 - first Spatial* or Spatial*DataFrame objects
##' @param  swdf2 - second Spatial* or Spatial*DataFrame objects (if omitted,  defaults to the same value as \code{swdf1})
##' @param  p - either a guess at the initial \eqn{p}-median set of a single integer indicating the size of the set (which is then chosen randomly)
##' @param  metric - the distance matrix (defaults to Euclidean computed via \code{euc.dists(swdf1,swdf2)} if not supplied)
##' @param  verbose - if TRUE print out each swap in the algorithm (default is FALSE)
##' @return List of nearest neigbour indices for each element from the \eqn{p}-median set
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' allocate(meuse,p=5)
##' 
##'
##' 
##' require(RColorBrewer)
##' require(GISTools)
##' data(georgia)
##' allocations.list <- allocate(georgia2,p=5)
##' zones <- gUnaryUnion(georgia2,allocations.list)
##' plot(zones,col=brewer.pal(5,"Accent"))
##' plot(georgia2,border=rgb(0,0,0,0.1),add=TRUE)
##' points(coordinates(georgia2)[allocations.list,],pch=16,cex=2,col=rgb(1,0.5,0.5,0.1))
##' 
##' @export
##' 

allocate <-function(swdf1,swdf2,p,metric,verbose=FALSE) {
  if (missing(swdf2))  swdf2 <- swdf1
  n.choices <- nrow(coordinates(swdf2)) 
  if (length(p) == 1) p <- sample(n.choices,p)
  if (missing(metric)) metric <- euc.dists(swdf1,swdf2)
  indices <- .tb(metric,p,verbose)
  nni <- .rviss(metric,indices)
  return(nni)
}

##' Creates the lines for a 'star diagram' 
##'  
##' @param  swdf1 - first Spatial* or Spatial*DataFrame objects
##' @param  swdf2 - second Spatial* or Spatial*DataFrame objects (if omitted,  defaults to the same value as \code{swdf1})
##' @param  alloc - a list saying which coordinate in swdf2 is allocated to each point in swdf1
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' allocations.list <- allocate(meuse,p=5)
##' star.lines <- star.diagram(meuse,alloc=allocations.list)
##' plot(star.lines)
##' @export
##' 


star.diagram <- function(swdf1,swdf2,alloc) {
  if (missing(swdf2))  swdf2 <- swdf1
  co1 <- coordinates(swdf1)
  co2 <- coordinates(swdf2)
  result <- vector(nrow(co1),mode='list')
  for (i in 1:nrow(co1))  result[[i]] <- Lines(list(Line(cbind(c(co1[i,1],co2[alloc[i],1]),c(co1[i,2],co2[alloc[i],2])))),ID=sprintf("Star%d",i))
  sl <- SpatialLines(result)
  sldf <- SpatialLinesDataFrame(sl,data.frame(allocate=alloc),match.ID=FALSE)
  return(sldf)}



