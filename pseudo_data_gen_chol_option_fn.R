pseudo_data_gen_chol <- function(grp_name){
  #Generate numbers randomly from uniform distribution
  n <- summary_stats[[grp_name]][1,'n']
  p <- length(all.vars(formula))
  desired_mean <- summary_stats[[grp_name]][summary_stats[[grp_name]][,'variable'] %in% all.vars(formula), 'mean']
  desired_cov_mat <- var_cov_mat[[grp_name]][all.vars(formula), all.vars(formula)]
  r <- matrix(rnorm(n * p), ncol = p)
  #Compute mean vector and covariance matrix of generated numbers
  r_mean <- apply(r, 2, mean)
  r_cov <- cov(r)
  #Compute centered values of r
  r_star <- r - matrix(rep(r_mean, rep(n,p)), ncol = p)
  #Compute scaled value of r
  r_0 <- solve(t(chol(r_cov))) %*% t(r_star)
  #Generate dataset with desired mean vector and covariance matrix
  ps <- t(chol(desired_cov_mat)) %*% r_0 + t(matrix(rep(desired_mean, rep(n,p)), ncol = p))
  pseudo_ <- as.data.frame(t(ps))
  names(pseudo_) <- all.vars(formula)
  pseudo_$clinic_name <- grp_name
  pseudo_
}