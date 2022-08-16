\name{haarWT}
\alias{haarWT}
\title{Function to apply the (univariate) Haar wavelet transform}

\usage{
haarWT(data)
}

\arguments{
\item{data}{Input univariate time series.}

}

\description{
This function applies the (univariate) Haar wavelet transform. For a time series containing missing values, the wavelet coefficients are generating and any NAs remain intact.
}
\value{
Returns a list containing the following elements:

\item{C}{Matrix containing the smooth coefficients for the transform.}
\item{D}{Matrix containing the detail coefficients for the transform.}
}

\examples{

set.seed(1)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

# compute the haar wavelet coefficients of the first time series component:

Xwt1 = haarWT(X[, 1])
}
