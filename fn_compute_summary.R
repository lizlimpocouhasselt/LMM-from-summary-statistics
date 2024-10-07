fn_compute_summary <- function(n, Xy.df){
  # Prepare functions to compute moments
  mv_moment_3_4_bypair <- function(a,b,n){
    a2b <- sum(a^2 * b)/n
    ab2 <- sum(a * b^2)/n
    ab3 <- sum(a * b^3)/n
    a2b2 <- sum(a^2 * b^2)/n
    a3b <- sum(a^3 * b)/n
    return(list(a2b=a2b, ab2=ab2, ab3=ab3, a2b2=a2b2, a3b=a3b))
  }
  mv_moment_3_4_by3 <- function(a,b,c,n){
    abc <- sum(a * b * c)/n
    a2bc <- sum(a^2 * b * c)/n
    ab2c <- sum(a * b^2 * c)/n
    abc2 <- sum(a * b * c^2)/n
    return(list(abc=abc, a2bc=a2bc, ab2c=ab2c, abc2=abc2))
  }
  mv_moment_4 <- function(a,b,c,d,n){
    abcd <- sum(a * b * c * d)/n
    abcd
  }
  
 
  p <- length(names(Xy.df)) #y and p-1 preds
  
  xy <- combn(names(Xy.df), 2)
  # Identify numeric variables
  numeric_var_names <- names(Xy.df)
  
  # Construct design matrix (as dataframe) from current dataset
  data_design_df <- Xy.df
  
  # Compute summary statistics (w/ moments)
  summary_stats <- data.frame(variable = names(Xy.df),
                              #type = c('bin', 'num', 'bin', 'num', 'bin', 'num', 'num', 'num', 'num', 'num', 'num', 'num'),
                              #type = rep('num', p),
                              n = n,
                              mean = apply(Xy.df, 2, mean, na.rm = TRUE),
                              variance = apply(Xy.df, 2, var),
                              target_3_moment = apply(Xy.df, 2, function(x) sum(x^3))/n,
                              target_4_moment = apply(Xy.df, 2, function(x) sum(x^4))/n,
                              row.names = NULL)
  
  var_cov_mat <- cov(Xy.df)
  
  mv_moment_3_4_bypair_df <- data.frame(vars = apply(xy, 2, function(x) paste0(t(x)[, 1], "_", t(x)[, 2])),
                                        a2b = sapply(1:dim(xy)[2], function(x) mv_moment_3_4_bypair(Xy.df[, xy[1,x]], Xy.df[, xy[2,x]], n)$a2b),
                                        ab2 = sapply(1:dim(xy)[2], function(x) mv_moment_3_4_bypair(Xy.df[, xy[1,x]], Xy.df[, xy[2,x]], n)$ab2),
                                        a3b = sapply(1:dim(xy)[2], function(x) mv_moment_3_4_bypair(Xy.df[, xy[1,x]], Xy.df[, xy[2,x]], n)$a3b),
                                        a2b2 = sapply(1:dim(xy)[2], function(x) mv_moment_3_4_bypair(Xy.df[, xy[1,x]], Xy.df[, xy[2,x]], n)$a2b2),
                                        ab3 = sapply(1:dim(xy)[2], function(x) mv_moment_3_4_bypair(Xy.df[, xy[1,x]], Xy.df[, xy[2,x]], n)$ab3))
  
  # Compute central moments of standardized numeric variables
  # summary_stats <- cbind(summary_stats,
  #                        std_mean = apply(summary_stats[,c('type','mean')], 1, function(x) ifelse(x[1] == 'num', 0, as.numeric(x[2]))),
  #                        std_var = apply(summary_stats[,c('type','variance')], 1, function(x) ifelse(x[1] == 'num', 1, as.numeric(x[2]))),
  #                        std_target_3_moment = apply(summary_stats[,c('type','variance','target_3_moment')], 1, function(x) ifelse(x[1] == 'num', as.numeric(x[3])/sqrt(as.numeric(x[2]))^3, as.numeric(x[3]))),
  #                        std_target_4_moment = apply(summary_stats[,c('type','variance','target_4_moment')], 1, function(x) ifelse(x[1] == 'num', as.numeric(x[3])/sqrt(as.numeric(x[2]))^4, as.numeric(x[3]))))
  
  # for(numeric_var_name in numeric_var_names){
  #   sd_orig <- sqrt(var_cov_mat[numeric_var_name, numeric_var_name])
  #   var_cov_mat[numeric_var_name, ] <- var_cov_mat[numeric_var_name, ]/sd_orig
  #   var_cov_mat[, numeric_var_name] <- var_cov_mat[, numeric_var_name]/sd_orig
  #   var_cov_mat[numeric_var_name, numeric_var_name] <- 1
  # }
  
  # for(numeric_var_name in numeric_var_names){
  #   var_orig <- summary_stats[summary_stats[,'variable'] == numeric_var_name, 'variance']
  #   mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'a2b'] <- mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'a2b']/var_orig
  #   mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'a2b'] <- mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'a2b']/sqrt(var_orig)
  #   mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'ab2'] <- mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'ab2']/sqrt(var_orig)
  #   mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'ab2'] <- mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'ab2']/var_orig
  #   mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'a3b'] <- mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'a3b']/sqrt(var_orig)^3
  #   mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'a3b'] <- mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'a3b']/sqrt(var_orig)
  #   mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'ab3'] <- mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'ab3']/sqrt(var_orig)
  #   mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'ab3'] <- mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'ab3']/sqrt(var_orig)^3
  #   mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'a2b2'] <- mv_moment_3_4_bypair_df[xy[1, ] == numeric_var_name, 'a2b2']/var_orig
  #   mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'a2b2'] <- mv_moment_3_4_bypair_df[xy[2, ] == numeric_var_name, 'a2b2']/var_orig
  # }
  
  if(p >= 3){
    xyz <- combn(names(Xy.df), 3)
    
    mv_moment_3_4_by3_df <- data.frame(vars = apply(xyz, 2, function(x) paste0(t(x)[, 1], "_", t(x)[, 2],"_", t(x)[, 3])),
                                       abc = sapply(1:dim(xyz)[2], function(x) mv_moment_3_4_by3(Xy.df[, xyz[1,x]], Xy.df[, xyz[2,x]], Xy.df[, xyz[3,x]], n)$abc),
                                       a2bc = sapply(1:dim(xyz)[2], function(x) mv_moment_3_4_by3(Xy.df[, xyz[1,x]], Xy.df[, xyz[2,x]], Xy.df[, xyz[3,x]], n)$a2bc),
                                       ab2c = sapply(1:dim(xyz)[2], function(x) mv_moment_3_4_by3(Xy.df[, xyz[1,x]], Xy.df[, xyz[2,x]], Xy.df[, xyz[3,x]], n)$ab2c),
                                       abc2 = sapply(1:dim(xyz)[2], function(x) mv_moment_3_4_by3(Xy.df[, xyz[1,x]], Xy.df[, xyz[2,x]], Xy.df[, xyz[3,x]], n)$abc2))
    
    # Compute moments of standardized numeric variables
    # for(numeric_var_name in numeric_var_names){
    #   var_orig <- summary_stats[summary_stats[,'variable'] == numeric_var_name, 'variance']
    #   mv_moment_3_4_by3_df[grepl(numeric_var_name, mv_moment_3_4_by3_df[,'vars']), 'abc'] <- mv_moment_3_4_by3_df[grepl(numeric_var_name, mv_moment_3_4_by3_df[,'vars']), 'abc']/sqrt(var_orig)
    #   mv_moment_3_4_by3_df[xyz[1, ] == numeric_var_name, c('ab2c', 'abc2')] <- mv_moment_3_4_by3_df[xyz[1, ] == numeric_var_name, c('ab2c', 'abc2')]/sqrt(var_orig)
    #   mv_moment_3_4_by3_df[xyz[1, ] == numeric_var_name, 'a2bc'] <- mv_moment_3_4_by3_df[xyz[1, ] == numeric_var_name, 'a2bc']/var_orig
    #   mv_moment_3_4_by3_df[xyz[2, ] == numeric_var_name, c('a2bc', 'abc2')] <- mv_moment_3_4_by3_df[xyz[2, ] == numeric_var_name, c('a2bc', 'abc2')]/sqrt(var_orig)
    #   mv_moment_3_4_by3_df[xyz[2, ] == numeric_var_name, 'ab2c'] <- mv_moment_3_4_by3_df[xyz[2, ] == numeric_var_name, 'ab2c']/var_orig
    #   mv_moment_3_4_by3_df[xyz[3, ] == numeric_var_name, c('a2bc', 'ab2c')] <- mv_moment_3_4_by3_df[xyz[3, ] == numeric_var_name, c('a2bc', 'ab2c')]/sqrt(var_orig)
    #   mv_moment_3_4_by3_df[xyz[3, ] == numeric_var_name, 'abc2'] <- mv_moment_3_4_by3_df[xyz[3, ] == numeric_var_name, 'abc2']/var_orig
    # }
  }
  
  if(p >= 4){
    wxyz <- combn(names(Xy.df), 4)
    
    mv_moment_4_df <- data.frame(vars = apply(wxyz, 2, function(x) paste0(t(x)[, 1], "_", t(x)[, 2],"_", t(x)[, 3], "_", t(x)[, 4])),
                                 abcd = sapply(1:dim(wxyz)[2], function(x) mv_moment_4(Xy.df[, wxyz[1,x]], Xy.df[, wxyz[2,x]], Xy.df[, wxyz[3,x]], Xy.df[, wxyz[4,x]], n)))
    
    # Compute moments of standardized numeric variables
    # for(numeric_var_name in numeric_var_names){
    #   var_orig <- summary_stats[summary_stats[,'variable'] == numeric_var_name, 'variance']
    #   mv_moment_4_df[grepl(numeric_var_name, mv_moment_4_df[,'vars']), 'abcd'] <- mv_moment_4_df[grepl(numeric_var_name, mv_moment_4_df[,'vars']), 'abcd']/sqrt(var_orig) 
    # }
  }
  
   

                               
  
  
  
  return(list(summary_stats,
              var_cov_mat,
              mv_moment_3_4_bypair_df,
              mv_moment_3_4_by3_df,
              mv_moment_4_df
              ))
  
  
  
}