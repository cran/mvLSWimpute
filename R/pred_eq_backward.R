# Function to form B matrix and RHS matrix for the one step ahead prediction (backcasting step)

pred_eq_backward <- function(lacv.array, p = 2, index){
  k <- sqrt(dim(lacv.array)[1])
  len <- dim(lacv.array)[2]-1
  L <- (dim(lacv.array)[3]-1)/2

  B <- matrix(0, k*p, k*p)
  RHS <- matrix(0, k*p, k)

  for(m in 1:p){
    for(n in 1:p){
      B[((m-1)*k+1):(k*m), ((n-1)*k+1):(k*n)] <- matrix(lacv.array[(1:(k*k)), len-m+1, L+n-m+1], nrow = k, ncol = k)
    }
  }

  for (m in 1:p){
    RHS[(((m-1)*k)+1):(m*k), 1:k] <- matrix(lacv.array[1:(k*k), len+1, L+m+1], nrow = k, ncol = k)
  }
  return(list(B = B, RHS = RHS))
}
