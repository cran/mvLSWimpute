\name{pdef}
\alias{pdef}
\title{Function to regularise the LWS matrix.}

\usage{
pdef(spec, W = 1e-10)
}

\arguments{
\item{spec}{EWS matrix that is to be regularised, can be either a 4D array or a \code{mvLSW} object.}

\item{W}{Tolerance in applying matrix regularisation to ensure each EWS matrix to be strictly positive definite. This is \code{1e-10} by default.}

}

\description{
This function regularises each EWS matrix to ensure that they are strictly positive definite, similar to the \code{mvEWS} function in the \pkg{mvLSW} package, except acting on a (bias-corrected) periodogram directly.  See
\code{\link[mvLSW]{mvEWS}} for more details.  Note: this function is not really intended to be used separately, but internally within the \code{spec_estimation} function.
}
\value{
Returns a \code{mvLSW} object containing the regularised EWS of a multivariate locally stationary time series.
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

set.seed(1)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

# form periodogram
tmp = apply(X, 2, function(x){haarWT(x)$D})

D = array(t(tmp), dim = c(2, 2^8, 8))

RawPer = array(apply(D, c(2, 3), tcrossprod), dim = c(2, 2, 2^8, 8))
RawPer = aperm(RawPer, c(1, 2, 4, 3))

# now correct

correctedper = correct_per(RawPer)

# now regularize

newper = pdef(correctedper)
}
