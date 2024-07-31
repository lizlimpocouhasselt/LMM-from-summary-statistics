pseudo_data_ls_fn <- function(grp_name){
  if(summary_stats[[grp_name]][1,'n'] == 1){
    pseudo_ <- as.data.frame(matrix(summary_stats[[grp_name]][summary_stats[[grp_name]][,'variable'] %in% all.vars(formula), 'mean'], nrow = 1))
    names(pseudo_) <- all.vars(formula)
  } else{
    pseudo_ <- as.data.frame(
      pseudo_data_gen(summary_stats[[grp_name]][1,'n'],
                      mu = summary_stats[[grp_name]][match(all.vars(formula), summary_stats[[grp_name]]$variable), 'mean'],
                      Sigma = var_cov_mat[[grp_name]][all.vars(formula), all.vars(formula)],
                      empirical = TRUE
      )
    )
  }
  pseudo_$clinic_name <- grp_name
  pseudo_
}