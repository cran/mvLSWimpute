# Function to form the lacv array for the forecasting step, where index is the index of the missing data
# and p.len is the number of terms to consider in the clipped predictor

form_lacv_forward <- function(spectrum, index, p.len = 2){
  P <- dim(spectrum)[1]

  # TT is the last observed timepoint, we wish to predict the value at time TT+1
  TT <- index-1
  filter <- wavethresh::filter.select(filter.number = 1, family = "DaubExPhase")
  Jp <- ceiling(logb(TT, 2))
  Nh <- length(filter$H != 0)
  K <- (2^(1:Jp) - 1) * (Nh - 1) + 1
  Jt <- max(which(K <= TT))

  # Fill the clipped spectrum with the spectral values from time TT-p.len+1 to TT (i.e. last p.len observed values)
  Smat <- array(NA, dim = c(P*P, length((TT-p.len+1):(TT+1)), Jt))
  for(p in 1:P){
    for (q in 1:P){
      for (j in 1:Jt) {
        Smat[(p-1) * P+1+(q-1), 1:(p.len+1), j] <- spectrum[p, q, j, ((TT-p.len+1):(TT+1))]
      }
    }
  }

  Psi <- PsiJmat(-Jt, filter.number = 1, family = "DaubExPhase")
  nc <- ncol(Psi)
  L <- (nc - 1)/2
  dimnames(Psi) <- list(NULL, c(-L:0, 1:L))

  # Form the auto and cross-covariance functions and save as an array
  lacv <- array(NA, dim = c((P*P), length((TT-p.len+1):(TT+1)), nc))
  for (i in 1:(P*P)){
    lacv[i, , ]  <- Smat[i, , ] %*% Psi
  }
  return(lacv)
}
