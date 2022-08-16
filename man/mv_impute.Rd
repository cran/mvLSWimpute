\name{mv_impute}
\alias{mv_impute}
\title{Function to apply the mvLSWimpute method and impute missing values in a multivariate locally stationary time series}

\usage{
mv_impute(data, p = 2, type = "forward", index = NULL)
}

\arguments{
\item{data}{Input multivariate time series, matrix of dimension TxP where P is the number of channels and T is the length of the series.}

\item{p}{The number of terms to include in the clipped predictor when carrying out one step ahead forecasting/backcasting.}

\item{type}{The type of imputation to carry out, either \code{"forward"} or \code{"forward-backward"} }

\item{index}{The set of time indices containing missing values, this is \code{NULL} by default and will be determined from the input series.}

}

\description{
This function applies the mvLSWimpute method to impute missing values in a multivariate locally stationary time series. The imputation can be based on forecasts only or use information from both a forecasting and backcasting step.
}
\value{
Returns a list containing the following elements:

\item{ImputedData}{Matrix containing the imputed time series.}
\item{missing.index}{Vector containing the set of time indices that have missing values.}
}
\references{
Wilson, R. E., Eckley, I. A., Nunes, M. A. and Park, T. (2021) A wavelet-based approach	for imputation in nonstationary multivariate time series.
_Statistics and Computing_ *31* Article 18, doi:10.1007/s11222-021-09998-2.
}
\author{
Rebecca Wilson
}
\note{
As with other time series imputation methods, \code{mv_impute} requires some data values at the start of the series.  In this case, we need 5 time points.
}
\examples{
set.seed(1)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

# create some fake missing data, taking care not to have missingness hear the start of the series
missing.index = sort(sample(10:2^8, 30))

X[missing.index, ] <- NA

newdata = mv_impute(X)

}
