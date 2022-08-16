\name{smooth_per}
\alias{smooth_per}
\title{Function to smooth the raw wavelet periodogram using the default \code{mvLSW} routine.}

\usage{
smooth_per(RawPer, type = "all", kernel.name="daniell", optimize = FALSE, kernel.param =
                 NULL, smooth.Jset = NULL)
}

\arguments{
\item{RawPer}{Raw wavelet periodogram that is to be smoothed, can be either a 4D array or a \code{mvLSW} object.}

\item{type}{Determines the type of smoothing to be performed, if \code{"all"} then the same smoothing kernel is applied to all levels, if \code{"by.level"} then a different smoothing kernel is applied to each level.}

\item{kernel.name}{Name of the smoothing kernel to be applied.}

\item{optimize}{Should the smoothing be optimized. If \code{FALSE} then smoothing is carried out with \code{kernel.name} and \code{kernel.param}.}

\item{kernel.param}{Value of the smoothing kernel parameter to be applied.}

\item{smooth.Jset}{Vector indicating which levels should be used in the calculation of the optimal kernel parameter. By default all levels are used.}
}

\description{
This function smooths the raw wavelet periodogram, similar to the \code{mvEWS} function in the \pkg{mvLSW} package, except acting on the raw periodogram directly.  See
\code{\link[mvLSW]{mvEWS}} for more details.  Note: this function is not really intended to be used separately, but internally within the \code{spec_estimation} function.
}

\value{
Returns a \code{mvLSW} object containing the smoothed EWS of a multivariate locally stationary time series.
}
\references{
Taylor, S.A.C., Park, T.A. and Eckley, I. (2019) Multivariate locally stationary wavelet analysis with the mvLSW R package.
_Journal of Statistical Software_ *90*(11) pp. 1-16, doi:10.18637/jss.v090.i11.\cr\cr

Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating time-evolving partial coherence between signals via multivariate
locally stationary wavelet processes. _IEEE Transactions on Signal Processing_ *62*(20) pp. 5240-5250.

}

\seealso{
	\code{\link[mvLSW]{mvEWS}}, \code{\link{spec_estimation}}
}
\author{
Rebecca Wilson
}
\examples{
## Sample bivariate locally stationary time series

set.seed(1)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

# form periodogram
tmp = apply(X, 2, function(x){haarWT(x)$D})
D = array(t(tmp), dim = c(2, 2^8, 8))

#sqrv <- function(d) return( d %*% t(d) )

#RawPer = array(apply(D, c(2, 3), sqrv), dim = c(2, 2, 2^8, 8))
RawPer = array(apply(D, c(2, 3), tcrossprod), dim = c(2, 2, 2^8, 8))
RawPer = aperm(RawPer, c(1, 2, 4, 3))

# now smooth

smoothper = smooth_per(RawPer)

}
