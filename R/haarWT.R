# Function to carry out the (univariate) Haar wavelet transform, keeping any NAs intact
haarWT <- function(data){

  nlevels <- nlevelsWT(data)
  n <- length(data)
  total <- n*(nlevels+1)

  C <- rep(0,length.out = total)
  D <- rep(0,length.out = (nlevels*n))

  C[1:n]<-data

  for (j in 1:nlevels){
    h <- rep(0, ((2^(j-1))+1))
    h[1] <- h[((2^(j-1))+1)] <- 1
    C[((j*n)+1):((j+1)*n)] <- 2^{-1/2}*filter(x = C[(((j-1)*n)+1):((j)*n)], filter = h, method = "convolution", sides = 2, circular = TRUE)
  }

  for (j in 2:nlevels){
    C[((j*n)+1):((j+1)*n)] <- binhf::shift(C[((j*n)+1):((j+1)*n)], dir = "left", places = (2^(j-1))-1)
  }

  for (j in 1:nlevels){
    g <- rep(0,((2^(j-1))+1))
    g[1] <- -1
    g[((2^(j-1))+1)] <- 1
    D[(((j-1)*n)+1):((j)*n)] <- 2^{-1/2}*filter(x = C[(((j-1)*n)+1):((j)*n)], filter = g, method = "convolution", sides = 2, circular = TRUE)
  }

  for (j in 2:nlevels){
    D[(((j-1)*n)+1):((j)*n)] <- binhf::shift(D[(((j-1)*n)+1):((j)*n)], dir = "left", places = (2^(j-2)))
  }

  C <- matrix(C, nrow = (nlevels+1), byrow = TRUE)
  D <- matrix(D, nrow = nlevels, byrow = TRUE)

  return(list(C = C, D = D))

}
