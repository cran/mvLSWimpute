spec_estimation <-
function(data, interp = "linear"){
  if (!is.matrix(data))
    stop("Invalid input argument")
  if (is.ts(data)) {
    TIME <- time(data)
  }
  else if (is.zoo(data) || is.xts(data)) {
    TIME <- time(data)
    data <- as.ts(data)
  }
  else if (is.matrix(data)) {
    TIME <- 1:nrow(data)
    data <- as.ts(data)
  }
  else {
    stop("Invalid input argument")
  }
  if (!is.numeric(data))
    stop("Input data must be numeric")
  if (ncol(data) == 1)
    stop("Input data is univariate.")
  if (log2(nrow(data))%%1 != 0)
    stop("Invalid length of input data.")
  P <- dim(data)[2]
  T <- n <- dim(data)[1]
  J <- log2(T)

  ## define convenience functions for later

  # Function to generate the wavelet coefficient arrays for the multivariate signal
   raw_coefficients <- function(data){
     P <- ncol(data)
     T <- nrow(data)
     J <- log2(T)

     # Apply haar wavelet transform function to each channel of the signal
     WT_coefficients = apply(data, 2, function(x) list(haarWT(x)))

     # Form matrices containing the smooth coefficients for each channel of the signal C and the detail coefficients D
     C = t(sapply(WT_coefficients, function(x) as.numeric(t(x[[1]]$C))))
     D = t(sapply(WT_coefficients, function(x) as.numeric(t(x[[1]]$D))))

     # Reshape dimensions so that C and D are now arrays of dimension P x J+1 x T and P x J x T respectively
     C = aperm(array(t(C), dim = c(T, J+1, P)), c(3, 2, 1))
     D = aperm(array(t(D), dim = c(T, J, P)), c(3, 2, 1))

     return(list(C = C, D = D))

   }

  # interpolation function

  interpfn = function(rawcoeffs, opt = "linear"){
	
    # For the scales up to and including Jstar, linearly interpolate each level of the periodogram
    for (p in 1:P){
       for (j in 1:Jstar){
          rawcoeffs$C[p, j+1, ] <- na_interpolation(rawcoeffs$C[p, j+1, ], option = opt)
          rawcoeffs$D[p, j, ] <- na_interpolation(rawcoeffs$D[p, j , ], option = opt)
       }
     }

     rawcoeffs
   }

  # haar filter functions
    C_int = function(C){
      for (j in (Jstar+1):J){
        h <- rep(0, ((2^(j-1))+1))
        h[1] <- h[((2^(j-1))+1)] <- 1
        C[((j*n)+1):((j+1)*n)] <- 2^{-1/2}*filter(x = C[(((j-1)*n)+1):((j)*n)], filter = h, method = "convolution", sides = 2, circular = TRUE)
      }

      for (j in (Jstar+1):J){
        C[((j*n)+1):((j+1)*n)] <- binhf::shift(C[((j*n)+1):((j+1)*n)], dir = "left", places = (2^(j-2)))
      }
      return(C)}

    D_int = function(C,D){
      for (p in 1:P){
        for (j in (Jstar+1):J){
          g <- rep(0, ((2^(j-1))+1))
          g[1] <- -1
          g[((2^(j-1))+1)] <- 1
          D[(((j-1)*n)+1):((j)*n), p] <- 2^{-1/2}*filter(x = C[(((j-1)*n)+1):((j)*n), p], filter = g, method = "convolution", sides = 2, circular = TRUE)
        }
      }
      for(p in 1:P){
        for (j in (Jstar+1):J){
          D[(((j-1)*n)+1):((j)*n), p] <- binhf::shift(D[(((j-1)*n)+1):((j)*n),p], dir = "left", places = (2^(j-2)))
        }
      }
      return(D)
    }

   #################### 
   ## now do computation...
   ####################

  rawcoeffs <- raw_coefficients(data)

  # Determine the coarsest scale Jstar that contains some wavelet coefficient information

  Jrange = apply(apply(rawcoeffs$D, c(1, 2), function(x) any(!is.na(x))), 2, prod)
  Jstar = max(which(Jrange==1))

  rawcoeffs = interpfn(rawcoeffs, opt = interp)

  # Reshape the interpolated arrays into matrices to allow us to fill the empty coarser levels

   C = apply(rawcoeffs$C, 1, function(x) as.numeric(t(x)))
   D = apply(rawcoeffs$D, 1, function(x) as.numeric(t(x)))

  # Fill the coarser levels by applying the filter equations to the coarsest scale containing information
  if(Jstar<J){
    C = apply(C, FUN = C_int, MARGIN = 2)
    D = D_int(C, D)
  }

  # Reshape the detail coefficient matrix into an array
  D = aperm(array((D), dim=c(T, J, P)), c(3, 2, 1))

  # Form the raw wavelet periodogram from the empirical wavelet coefficient vector

  D <- apply(D, 2:3, tcrossprod)
  D <- array(D, dim = c(P, P, J, T))

    # Smooth and correct the raw wavelet periodogram before checking that each matrix is positive semi definite
    periodogram <- smooth_per(as.mvLSW(D, filter.number = 1, family = "DaubExPhase"), smooth.Jset = 1:J)
    periodogram <- correct_per(periodogram)
    periodogram <- pdef(periodogram, W = 1e-10)

  return(periodogram)
}
