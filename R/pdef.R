# Function to ensure the spectrum is positive definite.

pdef <- function(spec, W = 1e-10){
  if(!is.mvLSW(spec)){
    spec = as.mvLSW(spec, filter.number = 1, family = "DaubExPhase")
  }

  AdjPositiveDef = getFromNamespace("AdjPositiveDef","mvLSW")

  spec = AdjPositiveDef(spec, W)
  return(spec)
}
