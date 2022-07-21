
log_likelihood = function(par, a = a,  log_x = log_x) {
  mu = par[1]
  sigma = par[2]
  #log_x = log(x)
  n = length(log_x)
  nz_log_x = log_x[is.finite(log_x)]
  nz = sum(is.finite(log_x))
  ll = (n - nz) * log(pnorm(q = a, mean = mu, sd = sigma)) - 0.5 * sum((nz_log_x - mu)^2 / sigma^2) - nz * log(sigma)
  return(-ll)
}

z_d = function(z, mean = 0, sd = 1){
  return(z * dnorm(z,  mean, sd))
}

precision_recall<- function (path, theta, verbose = TRUE, plot = FALSE, flip = T) {
  gcinfo(verbose = FALSE)
  ROC = list()
  theta = as.matrix(theta)
  d = ncol(theta)
  pos.total = sum(theta != 0)
  neg.total = d * (d - 1) - pos.total
  if (verbose)
    message("Computing F1 scores, false positive rates and true positive rates....", appendLF=FALSE)
  ROC$prec = rep(0, length(path))
  ROC$rec  = rep(0, length(path))
  ROC$F1 = rep(0, length(path))
  for (r in 1:length(path)) {
    tmp = as.matrix(path[[r]])
    tp.all = (theta != 0) * (tmp != 0)
    diag(tp.all) = 0
    ROC$tp[r] <- sum(tp.all != 0)/pos.total
    fp.all = (theta == 0) * (tmp != 0)
    diag(fp.all) = 0
    ROC$fp[r] <- sum(fp.all != 0)/neg.total
    fn = 1 - ROC$tp[r]
    precision = ROC$tp[r]/(ROC$tp[r] + ROC$fp[r])
    recall = ROC$tp[r]/(ROC$tp[r] + fn)

    # Correction
    tp = sum(tp.all)
    fp = sum(fp.all)
    p = sum(theta != 0)

    denominateur = tp + fp
    denominateur = max(denominateur, 1e-6)

    precision = tp / (denominateur)
    recall = tp / p

    ROC$prec[r] <- precision
    ROC$rec[r]  <- recall

    ROC$F1[r] = 2 * precision * recall/(precision + recall)
    if (is.na(ROC$F1[r]))
      ROC$F1[r] = 0
  }
  if (verbose)
    message("done.")
  rm(precision, recall, tp.all, fp.all, path, theta, fn)
  gc()
  if(flip) {
    ord.p = order(ROC$rec, na.last=NA)
    tmp2 = ROC$prec[ord.p]
    tmp1 = ROC$rec[ord.p]
    tmp2 = c(1, tmp2)
    tmp1 = c(0, tmp1)
  }
  else{
    ord.p = order(ROC$prec, na.last=NA)
    tmp1 = ROC$prec[ord.p]
    tmp2 = ROC$rec[ord.p]
  }

  if (plot) {
    par(mfrow = c(1, 1))
    plot(tmp1, tmp2, type = "b", main = "PR Curve", xlab = "Precision",
         ylab = "Recall", ylim = c(0, 1))
  }
  ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
  rm(ord.p, tmp1, tmp2)
  gc()
  class(ROC) = "roc"
  return(ROC)
}


rmvzinegbin_new <- function(n, mu, Sigma, munbs, ks, ps, ...) {
  # Generate an NxD matrix of Zero-inflated poisson data,
  # with counts approximately correlated according to Sigma
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(munbs) || missing(ps) || missing(ks)) {
    if (length(mu) != length(SDs)) stop("Sigma and mu dimensions don't match")
    munbs <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getLam(mu[i], SDs[i])))
    ps   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getP(mu[i], munbs[i])))
    ks   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getK(mu[i], SDs[i], munbs[i])))
  }
  if (length(munbs) != length(SDs)) stop("Sigma and mu dimensions don't match")
  d   <- length(munbs)


  normd  <- SpiecEasi::rmvnorm(n, rep(0, d), Sigma=Cor)
  unif   <- pnorm(normd)

  data <- t(matrix(VGAM::qzinegbin(t(unif), munb=munbs, size=ks, pstr0=ps, ...), d, n))
  return(data)
}

max_off_diagonal_value = function(S) {
  S_diag_off = S
  diag(S_diag_off) = 0
  max_value = max(S_diag_off)
  return(max_value)
}

lamda_path = function(lamda_max = 1, lamda_ratio = 0.1, nlambda=12){
  pen = seq(0,nlambda - 1, 1) / (nlambda - 1)
  lamda_min = lamda_max * lamda_ratio
  pen = rep(lamda_max, nlambda) * (lamda_min)^pen
  return(pen)
}


# from CRAN mirror on github since package is no longer available on CRAN
QUIC <- function(S, rho, path = NULL, tol = 1.0e-4, msg = 1, maxIter = 1000,
                 X.init = NULL, W.init = NULL) {
  ### $Id: QUIC.R,v 1.7 2012-05-01 02:12:19 sustik Exp $

  n <- nrow(S)
  if (is.null(path))
    npath <- 1
  else
    npath <- length(path)
  if (!is.matrix(rho) && length(rho) != 1 && length(rho) != n) {
    stop("Wrong number of elements in rho")
  }
  if (is.vector(rho)){
    rho <- matrix(sqrt(rho))%*%sqrt(rho)
  }
  if (length(rho) == 1){
    rho <- matrix(rho, ncol = n, nrow = n)
  }
  if (is.null(path)) {
    if (is.null(X.init)) {
      X <- diag(n)
      W <- diag(n)
    } else {
      X <- X.init
      W <- W.init
    }
  } else {
    if (is.null(X.init)) {
      X <- array(diag(n), c(n, n, npath))
      W <- array(diag(n), c(n, n, npath))
    } else {
      X <- array(0, c(n, n, npath))
      W <- array(0, c(n, n, npath))
      X[, , 1] <- X.init
      W[, , 1] <- W.init
    }
  }
  opt <- matrix(0, ncol = npath, nrow = 1)
  cputime <- matrix(0, ncol = npath, nrow = 1)
  iter <- matrix(0, ncol = npath, nrow = 1)
  dGap <- matrix(0, ncol = npath, nrow = 1)
  if (is.null(path))
    job <- "d"
  else
    job <- "p"

  storage.mode(job) <- "character"
  storage.mode(S) <- "double"
  storage.mode(rho) <- "double"
  storage.mode(npath) <- "integer"
  storage.mode(path) <- "double"
  storage.mode(tol) <- "double"
  storage.mode(msg) <- "integer"
  storage.mode(maxIter) <- "integer"
  storage.mode(X) <- "double"
  storage.mode(W) <- "double"
  storage.mode(opt) <- "double"
  storage.mode(cputime) <- "double"
  storage.mode(iter) <- "integer"
  storage.mode(dGap) <- "double"
  # tmp<-.C("QUICR",
  #         job, n, S, rho, npath, path, tol, msg, maxIter,
  #          X = X, W = W, opt = opt, cputime = cputime, iter = iter,
  #          dGap = dGap)
  return (list(X = X, W = W, opt = opt, cputime = cputime,
               iter = iter, regloglik = -(n/2)*opt,
               dGap = dGap))
}
