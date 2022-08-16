\name{pred_eq_forward}
\alias{pred_eq_forward}
\alias{pred_eq_backward}
\title{Function to form the prediction equations for the forecasting / backcasting step.}

\usage{
pred_eq_forward(lacv.array, p = 2, index)
pred_eq_forward(lacv.array, p = 2, index)
}

\arguments{
\item{lacv.array}{The local autocovariance array from which we want to form the prediction equations, can be obtained as the output of the \code{form_lacv_forward} or \code{form_lacv_backward} functions.}

\item{p}{Number of terms to include in the clipped predictor when forecasting / backcasting.}

\item{index}{Time index of the missing data which we wish to impute.}

}

\description{
This function generates the prediction equations (B matrix and RHS matrix) for one step ahead prediction.
}

\details{
The one-step ahead predictor is formed as a linear combination of the time series.  The coefficients involved in optimal predictor (in the sense of minimising the mean square prediction error) are obtained by
solving a matrix equation formed using parts of the (estimated) local autocovariance array.  This function forms the matrices involved in the equation used to find the optimal linear predictor.  See equation (6) in Wilson et al. (2021) or Section 3.3 in Fryzlewicz et al. (2003) for more details. 
}


\value{
Returns a list containing the following elements:
\item{B}{The left-hand side of the matrix equation to compute the optimal one-step ahead predictor, which is essentially used to approximate the MSPE for a particular set of coefficients used in a predictor.}
\item{RHS}{The right hand side of the matrix equation used to compute the optimal one-step ahead predictor.}
}
\references{
Wilson, R. E., Eckley, I. A., Nunes, M. A. and Park, T. (2021) A wavelet-based approach for imputation in nonstationary multivariate time series.
_Statistics and Computing_ *31* Article 18, doi:10.1007/s11222-021-09998-2.\cr\cr

Fryzlewicz, P. van Bellegem, S. and von Sachs, R. (2003) Forecasting non-stationary time series by wavelet process modelling. _Annals of the Institute of Statistical Mathematics_ *55* (4), pp. 737-764.
}
\seealso{
\code{\link{form_lacv_forward}}, \code{\link{pred_eq_backward}}
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

# obtain the LACV

lacvfor <- form_lacv_forward(spec$spectrum, missing.index[1], p.len = 2)

# form matrix equation terms
mspeterms = pred_eq_forward(lacvfor, p = 2, missing.index[1])

# compute the optimal coefficients in the linear predictor:
predcoeffs = solve(mspeterms$B, mspeterms$RHS)

}
