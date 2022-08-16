\name{spec_estimation}
\alias{spec_estimation}
\title{Function to estimate the Local Wavelet Spectral matrix for a multivariate locally stationary time series containing missing values}

\usage{
spec_estimation(data, interp = "linear")
}

\arguments{
\item{data}{Input multivariate time series, matrix of dimension TxP where P is the number of channels and T is the length of the series.}
\item{interp}{Method of interpolation of NAs in spectrum.  Can be \code{"linear"} or \code{"spline"}; see \code{\link[imputeTS]{na_interpolation}} for more detals.  See also note below.}
}

\description{
This function estimates the LWS matrix for a multivariate locally stationary time series containing missing values. If the input time series does not contain missing values then spectral estimation is carried out using routines from the \pkg{mvLSW} package.
}
\value{
Returns a \code{mvLSW} object containing the estimated LWS matrix.
}

\seealso{
\code{\link{correct_per}}, \code{\link{smooth_per}}, \code{\link[mvLSW]{mvEWS}}, \code{\link[imputeTS]{na_interpolation}}
}
\note{
	For some series with a lot of missing values, the linear interpolation will result in zero periodogram values (due to the form of the Haar filters).  This may not be desirable, so a higher order (spline) interpolation function may be better. 
}
\examples{

## Sample bivariate locally stationary time series

set.seed(1)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

# create some missing values, taking care to provide some data at the start of the series 

missing.index = sort(sample(10:2^8, 30))

X[missing.index, ] <-NA

# estimate the spectrum

spec = spec_estimation(X)

}
