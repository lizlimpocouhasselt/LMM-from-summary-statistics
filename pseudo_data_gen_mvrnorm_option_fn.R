pseudo_data_gen_mvrnorm <- function(grp_name){
  pseudo_ <- as.data.frame(
    mvrnorm(summary_stats[[grp_name]][1,'n'],
            mu = summary_stats[[grp_name]][match(all.vars(formula), summary_stats[[grp_name]]$variable), 'mean'],
            Sigma = var_cov_mat[[grp_name]][all.vars(formula), all.vars(formula)],
            empirical = TRUE
    )
  )
  pseudo_$clinic_name <- grp_name
  pseudo_
}