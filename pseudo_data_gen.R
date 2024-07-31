pseudo_data_gen <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (EISPACK) 
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(runif(p * n), n) #changed from rnorm to runif
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- apply(X, 2, function(column) ifelse(is.na(column), 0, column))
    X <- X %*% svd(X, nu = 0, nv = p)$v
    X <- scale(X, FALSE, TRUE)
    X <- apply(X, 2, function(column) ifelse(is.na(column), 0, column))
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) drop(X) else t(X)
}
