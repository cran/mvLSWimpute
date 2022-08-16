\name{form_lacv_forward}
\alias{form_lacv_forward}
\alias{form_lacv_backward}
\title{Function to form the local autocovariance array for the forecasting / backcasting step.}

\usage{
form_lacv_forward(spectrum, index, p.len = 2)
form_lacv_backward(spectrum, index, p.len = 2)
}

\arguments{
\item{spectrum}{Local wavelet spectral matrix for which we wish to form the local autocovariance array.}

\item{index}{Time index of the missing data which we wish to impute.}

\item{p.len}{Number of terms to include in the clipped predictor when forecasting / backcasting.}

}

\description{
This function generates the local autocovariance (LACV) array that is used in the forecasting / backcasting step to form the prediction equations.
}

\details{
In order to form the one-step ahead predictor for use in the imputation algorithm of Wilson et al. (2021), one needs the local autocovariance (LACV).  This is computed using the relationship between the LACV and the local wavelet spectrum (LWS).  See equations (4) and (5) in Wilson et al. (2021) for more details.
}

\value{
Returns the local autocovariance array that can be used as an input to the \code{pred_eq_forward} or \code{pred_eq_backward} function.
}
\references{
Wilson, R. E., Eckley, I. A., Nunes, M. A. and Park, T. (2021) A wavelet-based approach for imputation in nonstationary multivariate time series.
_Statistics and Computing_ *31* Article 18, doi:10.1007/s11222-021-09998-2.\cr\cr

Taylor, S.A.C., Park, T.A. and Eckley, I. (2019) Multivariate locally stationary wavelet analysis with the mvLSW R package.
_Journal of Statistical Software_ *90*(11) pp. 1-16, doi:10.18637/jss.v090.i11.\cr\cr

Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating time-evolving partial coherence between signals via multivariate
locally stationary wavelet processes. _IEEE Transactions on Signal Processing_ *62*(20) pp. 5240-5250.

}
\seealso{
\code{\link{pred_eq_forward}}, \code{\link{pred_eq_backward}}
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

# create some missing values, taking care to provide some data at the start of the series 

missing.index = sort(sample(10:2^8, 30))

X[missing.index, ] <-NA

# estimate the spectrum

spec = spec_estimation(X)

out <- form_lacv_forward(spec$spectrum, missing.index[1], p.len=2)

}
