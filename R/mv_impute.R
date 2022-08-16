# Function to carry out the MvLSWimpute method
# p is the number of terms in the clipped predictor, index is the set of indices where datapoints are missing

mv_impute <- function(data, p = 2, type = "forward", index = NULL){
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
  if(is.null(type)==TRUE){
    stop("Imputation type must be either forward or forward-backward.")
  }

  P <- dim(data)[2]
  T <- dim(data)[1]
  J <- log2(T)

  # Finds locations of missing datapoints in the original signal if these are not provided
  if(is.null(index)){
    missing.index <- which(apply(data, 1, function(x){any(is.na(x))}))
  }else{
    missing.index <- index
  }
  missing.index <- sort(missing.index)

  # Estimates the LWS spectrum of the original series
  spec <- spec_estimation(data)

  if(type=="forward"){
    # For each time where we have missing values, use one step ahead prediction to obtain an estimate
    len.missing = length(missing.index)
    data.forecast <- data
    data.update <- data
    for (l in 1:len.missing){
      lacv.forward <- form_lacv_forward(spec$spectrum, missing.index[l], p.len = p)
      pred.forward <- pred_eq_forward(lacv.array = lacv.forward, p = p, index = missing.index[l])
      b.forward <- solve(pred.forward$B, pred.forward$RHS)
      pmean.forward = t(b.forward) %*% c(t(data.forecast[((missing.index[l]-1):(missing.index[l]-p)), ]))
      update.index.forward <- which(is.na(data[missing.index[l], ]))
      data.forecast[missing.index[l], update.index.forward] <- pmean.forward[update.index.forward]
    }
    data.update <- data.forecast
  }

  if(type=="forward-backward"){
  # For each time where we have missing values, use one step ahead prediction to obtain an estimate
  len.missing = length(missing.index)
  data.forecast <- data
  data.backcast <- data
  data.update <- data
  for (l in 1:len.missing){
    lacv.forward <- form_lacv_forward(spec$spectrum, missing.index[l], p.len=p)
    lacv.backward <- form_lacv_backward(spec$spectrum, missing.index[len.missing+1-l], p.len = p)

    pred.forward <- pred_eq_forward(lacv.array = lacv.forward, p = p, index = missing.index[l])
    pred.backward <- pred_eq_backward(lacv.array = lacv.backward, p = p, index = missing.index[len.missing+1-l])

    b.forward <- solve(pred.forward$B, pred.forward$RHS)
    b.backward <- solve(pred.backward$B, pred.backward$RHS)

    pmean.forward = t(b.forward) %*% c(t(data.forecast[((missing.index[l]-1):(missing.index[l]-p)), ]))
    pmean.backward = t(b.backward) %*% c(t(data.backcast[((missing.index[len.missing+1-l]+1):(missing.index[len.missing+1-l]+p)), ]))

    update.index.forward <- which(is.na(data[missing.index[l], ]))
    update.index.backward <- which(is.na(data[missing.index[len.missing+1-l], ]))

    data.forecast[missing.index[l], update.index.forward] <- pmean.forward[update.index.forward]
    data.backcast[missing.index[len.missing+1-l], update.index.backward] <- pmean.backward[update.index.backward]
  }
  data.update <- as.ts(matrix(as.numeric(unlist((0.5*(data.forecast+data.backcast)))), nrow = T, ncol = P))
  }

  return(list(ImputedData = data.update, missing.index = missing.index))
}
