infer_Z = function(X, seq_depth = "TS") {

  d = dim(X)[2] # number of rows
  n = dim(X)[1] # number of columns

  X_clr = matrix(0, n, d) # makes a zero matrix
  nz = X != 0 # n by d logical matrix of whether entry is zero

  if( seq_depth == "TS"){
    X_n = X / apply(X, 1, sum) # divide sum row totals
    for(i in 1:n) {
      X_clr[i,] = log(X_n[i,]) - mean(log(X_n[i,X_n[i,] != 0])) # log standardizing?
    }
  }
  if(seq_depth == "unif"){
    for(k in 1:d){
      nz_X_k = X[nz[,k],k]
      emp = ecdf(c(nz_X_k, max(nz_X_k) + 1))
      X_clr[,k] = qnorm(emp(X[,k]))
    }
    X_clr[!nz] = -Inf
  }

  log_ratios = X_clr
  a = double(n)

  for(i in 1:d){
    a[i] = min(log_ratios[nz[,i],i]) # what is this doing? Looking at the minimum value of non-zero
  }

  Z = matrix(0, n, d)
  Z[] = log_ratios

  # here's the meat of the function!!
  for(i in 1:d){
    if(sum(nz[,i]) > 1){ # only will do this for variables which have at least more than one nonzero val
      mu = mean(log_ratios[is.finite(log_ratios)[,i],i]) # mean of non-zero vals
     # general purpose optimization algorithm
       o = optim(c(0, 1), log_likelihood, log_x = X_clr[,i], a  =  a[i]) # this doesn't seem to work when there are no zeros in the input
      par = o$par
      mu = par[1] # why this second mu?
      sigma = par[2]
      exp_Z_0 = integrate(z_d, lower = -Inf, upper = a[i], mean = mu, sd =  sigma)
      Z_0 =  integrate(dnorm, lower = -Inf, upper = a[i], mean = mu, sd =  sigma)

      Z[X[,i] == 0, i] = exp_Z_0$value / Z_0$value # here for the values in the matrix that are 0, uses this estimate instead
      # Handle small sigma
      if(sigma < 1e-5){
        Z[X[,i] == 0, i] = -10
      }
    }
    else if(sum(nz[,i]) == 1){
      Z[X[,i] == 0, i] = -10 # why this -10 value if there is just one value in the column?
    }
    else{
      Z[,i] = 0
    }
  }
  return(Z)
}
