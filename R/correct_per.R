# Function to correct the raw wavelet periodogram using the inverse of the autocorrelation wavelet matrix for the Haar wavelet

correct_per <- function(RawPer){

  if (!is.mvLSW(RawPer)){
    RawPer = as.mvLSW(RawPer, filter.number = 1, family = "DaubExPhase")
  }
  P <- dim(RawPer$spectrum)[1]
  J <- dim(RawPer$spectrum)[3]
  T <- dim(RawPer$spectrum)[4]

  # Calculate inverse of the autocorrelation wavelet matrix for haar wavelets
  A <- ipndacw(-J, filter.number = 1, family = "DaubExPhase")
  Ainv <- solve(A)
  acorrect <- function(I){
    return(Ainv %*% I)
  }

  cI <- apply(RawPer$spectrum, c(1, 2), acorrect)
  cI <- aperm(cI, c(2, 3, 1))
  CorrectPer <- array(cI, dim = c(P, P, J, T))
  return(as.mvLSW(CorrectPer, filter.number = 1, family = "DaubExPhase"))
}
