# Function to smooth the raw wavelet periodogram using the smoothing function contained in the MvLSW package

smooth_per <- function(RawPer, type = "all", kernel.name = "daniell", optimize = FALSE, kernel.param = NULL, smooth.Jset = NULL){
  if (!is.mvLSW(RawPer)){
    RawPer = as.mvLSW(RawPer, filter.number = 1, family = "DaubExPhase")
  }
  P <- dim(RawPer$spectrum)[1]
  J <- dim(RawPer$spectrum)[3]
  T <- dim(RawPer$spectrum)[4]

  if(is.null(kernel.param)){
    kernel.param = floor(sqrt(T))
  }

  if(is.null(smooth.Jset)){
    smooth.Jset=1:J
  }

  Smooth_EWS=utils::getFromNamespace("Smooth_EWS", "mvLSW")

  SmoothPer <- Smooth_EWS(RawPer, kernel.name, kernel.param, optimize, type, level = smooth.Jset)

  return(SmoothPer)
}
