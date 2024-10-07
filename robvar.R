#------------------------------------------------------
# COMPUTING ROBUST STD ERRORS FROM SUMMARY STATISTICS
#------------------------------------------------------




### ------ PRELIMINARIES ------
# Clear memory
rm(list=ls(all=TRUE))


# Load packages / functions
# library(pracma)
# library(dplyr)
# library(MASS)
library(lmtest)
library(sandwich)
library(medicaldata)
library(dplyr)
source('fn_compute_summary.R')
source('pseudo_data_gen_fn.R') #modified mvrnorm



### ------ DATA PROVIDER TASK ------

# Load actual data
data("covid_testing")
data <- covid_testing %>%
  dplyr::select(gender, ct_result, demo_group, pan_day, age, drive_thru_ind, col_rec_tat, rec_ver_tat)
Xy.df <- as.data.frame(model.matrix(~ ., data))
n <- nrow(Xy.df)


summary.all <- fn_compute_summary(n, Xy.df[, -1])
summary_stats <- summary.all[[1]]
var_cov_mat <- summary.all[[2]]
mv_moment_3_4_bypair_df <- summary.all[[3]]
mv_moment_3_4_by3_df <- summary.all[[4]]
mv_moment_4_df <- summary.all[[5]]


### ------ DATA ANALYST TASK ------

# Generate pseudo-data --------------------------------------
formula <- ct_result ~ gendermale + pan_day + age + drive_thru_ind
set.seed(121314)
mu <- summary_stats[match(all.vars(formula), summary_stats$variable), 'mean']
Sigma <- var_cov_mat[all.vars(formula), all.vars(formula)]
pseudo_data <- as.data.frame(pseudo_data_gen(n, mu, Sigma, empirical = T))

# Estimate a linear model
lm_pseudo <- lm(formula, pseudo_data)
betas <- lm_pseudo$coefficients

# I. Construct XtX --> first set of summary statistics --------------
X.ps <- model.matrix(formula, pseudo_data)  #will be equal to the actual X
XtX <- t(X.ps) %*% X.ps

# II. Construct XtWX -- second set of summary statistics
p <- length(all.vars(formula)[-1]) #no. of predictors

# A. XtWX.1 --------------------------------------------
XtWX.1 <- matrix(0, nrow = p + 1, ncol = p + 1) #initial matrix
# Begin with the diagonal
y.row.ind <- summary_stats$variable == all.vars(formula)[1]
y2 <- summary_stats[y.row.ind, 'variance'] * (n-1) +
  summary_stats[y.row.ind, 'mean']^2 * n

y2x2_j <- sapply(1:p, function(j){
  vars <- mv_moment_3_4_bypair_df$vars
  var1 <- all.vars(formula)[1]
  var2 <- all.vars(formula)[1+j]
  row.ind <- grepl(var1, vars, fixed = T) &
    grepl(var2, vars, fixed = T) &
    !grepl(paste(c(":","X"), collapse = "|"), vars)
  mv_moment_3_4_bypair_df[row.ind, 'a2b2']})*n

diag(XtWX.1) <- c(y2, y2x2_j)


# Continue with the upper triangular part by transposing the lower triangular
y2x <- sapply(1:p, function(j){
  vars <- mv_moment_3_4_bypair_df$vars
  var1 <- all.vars(formula)[1]
  var2 <- all.vars(formula)[1+j]
  row.ind <- grepl(var1, vars, fixed = T) &
    grepl(var2, vars, fixed = T) &
    !grepl(paste(c(":","X"), collapse = "|"), vars)
  var.name <- vars[row.ind]
  if(unlist(gregexpr(var1, var.name)) == 1){
    mv_moment_3_4_bypair_df[row.ind, 'a2b']*n
    } else{ mv_moment_3_4_bypair_df[row.ind, 'ab2']*n}
})
y2x_jx_k <- unlist(sapply(1:(p-1), function(k){
  sapply((k+1):p, function(j){
    vars <- mv_moment_3_4_by3_df$vars
    var1 <- all.vars(formula)[1]
    var2 <- all.vars(formula)[1+k]
    var3 <- all.vars(formula)[1+j]
    row.ind <- grepl(var1, vars, fixed = T) &
      grepl(var2, vars, fixed = T) &
      grepl(var3, vars, fixed = T) &
      !grepl(paste(c(":","X"), collapse = "|"), vars)
    var.name <- vars[row.ind]
    if(unlist(gregexpr(var1, var.name)) == 1){
      mv_moment_3_4_by3_df[row.ind, 'a2bc']*n
    } else if(unlist(gregexpr(var1, var.name)) +
            nchar(var1) - 1 == nchar(var.name)){
      mv_moment_3_4_by3_df[row.ind, 'abc2']*n
    } else{mv_moment_3_4_by3_df[row.ind, 'ab2c']*n}
  })
})) 

XtWX.1[lower.tri(XtWX.1)] <- c(y2x, y2x_jx_k)
XtWX.1 <- t(XtWX.1)

# Complete the lower triangular part
XtWX.1[lower.tri(XtWX.1)] <- c(y2x, y2x_jx_k)


# B. XtWX.2 --------------------------------------------------------
XtWX.2 <- matrix(0, nrow = p + 1, ncol = p + 1) #initial matrix
# Begin with the diagonal
yx <- sapply(1:(p+1), function(j){
  var1 <- all.vars(formula)[1]
  vars <- summary_stats$variable
  if(j == 1){
    summary_stats[vars == var1, 'mean'] * n
  } else{
    var2 <- all.vars(formula)[j]
    var_cov_mat[var1, var2] * (n-1) +
      n * summary_stats[vars == var1, 'mean'] *
      summary_stats[vars == var2, 'mean']
  }}) 
yxTb <- t(yx * 2) %*% as.matrix(betas)

x2_jyx <- sapply(1:p, function(j){
  var1 <- all.vars(formula)[1]
  var2 <- all.vars(formula)[1+j]
  x2_jy_x_k <- sapply(1:(p+1), function(k){
    var3 <- all.vars(formula)[k]
    if(k == 1 | var2 == var3){
      vars <- mv_moment_3_4_bypair_df$vars
      row.ind <- grepl(var1, vars, fixed = T) &
        grepl(var2, vars, fixed = T) &
        !grepl(paste(c(":","X"), collapse = "|"), vars)
      var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
      if(k == 1){
        if(unlist(gregexpr(var1, var.name)) == 1){
          mv_moment_3_4_bypair_df[row.ind, 'ab2']*n
        } else{mv_moment_3_4_bypair_df[row.ind, 'a2b']*n}
      } else{
        if(unlist(gregexpr(var1, var.name)) == 1){
          mv_moment_3_4_bypair_df[row.ind, 'ab3']*n
        } else{mv_moment_3_4_bypair_df[row.ind, 'a3b']*n}
      }
    } else{
      vars <- mv_moment_3_4_by3_df$vars
      row.ind <- grepl(var1, vars, fixed = T) &
        grepl(var2, vars, fixed = T) &
        grepl(var3, vars, fixed = T) &
        !grepl(paste(c(":","X"), collapse = "|"), vars)
      var.name <- mv_moment_3_4_by3_df$vars[row.ind]
      if(unlist(gregexpr(var2, var.name)) == 1){
        mv_moment_3_4_by3_df[row.ind, 'a2bc']*n
      } else if(unlist(gregexpr(var2, var.name)) +
                nchar(var2) - 1 == nchar(var.name)){
        mv_moment_3_4_by3_df[row.ind, 'abc2']*n
      } else{mv_moment_3_4_by3_df[row.ind, 'ab2c']*n}
    }
  })
  return(t(x2_jy_x_k*2) %*% as.matrix(betas))
})

diag(XtWX.2) <- c(yxTb, x2_jyx)

# Continue with the upper triangular part by transposing the lower triangular
x_jyx <- sapply(1:p, function(j){
  var1 <- all.vars(formula)[1]
  var2 <- all.vars(formula)[1+j]
  x_jy_x_k <- sapply(1:(p+1), function(k){
    var3 <- all.vars(formula)[k]
    if(k == 1){
      vars <- summary_stats$variable
      var_cov_mat[var1, var2] * (n-1) +
        n * summary_stats[vars == var1, 'mean'] *
        summary_stats[vars == var2, 'mean'] 
    } else if(var2 == var3){
      vars <- mv_moment_3_4_bypair_df$vars
      row.ind <- grepl(var1, vars, fixed = T) &
        grepl(var2, vars, fixed = T) &
        !grepl(paste(c(":","X"), collapse = "|"), vars)
      var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
      if(unlist(gregexpr(var1, var.name)) == 1){
        mv_moment_3_4_bypair_df[row.ind, 'ab2']*n
      } else{mv_moment_3_4_bypair_df[row.ind, 'a2b']*n}
    } else{
      vars <- mv_moment_3_4_by3_df$vars
      row.ind <- grepl(var1, vars, fixed = T) &
        grepl(var2, vars, fixed = T) &
        grepl(var3, vars, fixed = T) &
        !grepl(paste(c(":","X"), collapse = "|"), vars)
      var.name <- mv_moment_3_4_by3_df$vars[row.ind]
      mv_moment_3_4_by3_df[row.ind, 'abc']*n
    }
  })
  return(t(x_jy_x_k*2) %*% as.matrix(betas))
})

x_jx_kyx <- unlist(sapply(1:(p-1), function(j){
  sapply((j+1):p, function(k){
    var1 <- all.vars(formula)[1]
    var2 <- all.vars(formula)[1+j]
    var3 <- all.vars(formula)[1+k]
    x_jx_ky_x_l <- sapply(1:(p+1), function(l){
      var4 <- all.vars(formula)[l]
      if(l == 1){
        vars <- mv_moment_3_4_by3_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          grepl(var3, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        mv_moment_3_4_by3_df[row.ind, 'abc']*n
      } else if(var2 == var4){
        vars <- mv_moment_3_4_by3_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          grepl(var3, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        var.name <- mv_moment_3_4_by3_df$vars[row.ind]
        if(unlist(gregexpr(var2, var.name)) == 1){
          mv_moment_3_4_by3_df[row.ind, 'a2bc']*n
        } else if(unlist(gregexpr(var2, var.name)) +
                  nchar(var2) - 1 == nchar(var.name)){
          mv_moment_3_4_by3_df[row.ind, 'abc2']*n
        } else{mv_moment_3_4_by3_df[row.ind, 'ab2c']*n}
      } else if(var3 == var4){
        vars <- mv_moment_3_4_by3_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          grepl(var3, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        var.name <- mv_moment_3_4_by3_df$vars[row.ind]
        if(unlist(gregexpr(var3, var.name)) == 1){
          mv_moment_3_4_by3_df[row.ind, 'a2bc']*n
        } else if(unlist(gregexpr(var3, var.name)) +
                  nchar(var3) - 1 == nchar(var.name)){
          mv_moment_3_4_by3_df[row.ind, 'abc2']*n
        } else{mv_moment_3_4_by3_df[row.ind, 'ab2c']*n}
      } else{
        vars <- mv_moment_4_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          grepl(var3, vars, fixed = T) &
          grepl(var4, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        mv_moment_4_df[row.ind, 'abcd']*n
      }
    })
    return(t(x_jx_ky_x_l*2) %*% as.matrix(betas))
  })
}))

XtWX.2[lower.tri(XtWX.2)] <- c(x_jyx, x_jx_kyx)

XtWX.2 <- t(XtWX.2)

# Complete the lower triangular part
XtWX.2[lower.tri(XtWX.2)] <- c(x_jyx, x_jx_kyx)





# C. XtWX.3 --------------------------------------------------------
XtWX.3 <- matrix(0, nrow = p + 1, ncol = p + 1) #initial matrix
# Begin with the diagonal
bTxxTb <- t(betas)%*% XtX %*% betas
bTx2_jxxTb <- sapply(1:p, function(j){
  var1 <- all.vars(formula)[1+j]
  x2_jXtX <- matrix(0, nrow = p + 1, ncol = p + 1) #initial matrix
  diag(x2_jXtX) <- sapply(1:(p+1), function(k){
    if(k == 1){
      vars <- summary_stats$variable
      var_cov_mat[var1, var1]*(n-1) + 
        summary_stats[vars == var1, 'mean']^2*n
    } else{
      var2 <- all.vars(formula)[k]
      if(var1 == var2){
        vars <- summary_stats$variable
        summary_stats[vars == var1, 'target_4_moment']*n
      } else{
        vars <- mv_moment_3_4_bypair_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        mv_moment_3_4_bypair_df[row.ind, 'a2b2']*n
      }
    }
  })

  x2_jx_k <- sapply(1:p, function(k){
    var2 <- all.vars(formula)[1+k]
    if(var1 == var2){
      vars <- summary_stats$variable
      summary_stats[vars == var1, 'target_3_moment'] * n
    } else{
      vars <- mv_moment_3_4_bypair_df$vars
      row.ind <- grepl(var1, vars, fixed = T) &
        grepl(var2, vars, fixed = T) &
        !grepl(paste(c(":","X"), collapse = "|"), vars)
      var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
      if(unlist(gregexpr(var1, var.name)) == 1){
        mv_moment_3_4_bypair_df[row.ind, 'a2b']*n
      } else{mv_moment_3_4_bypair_df[row.ind, 'ab2']*n}
    }
  })
  x2_jx_kx_l <- unlist(sapply(1:(p-1), function(k){
    var2 <- all.vars(formula)[k+1]
    sapply((k+1):p, function(l){
      var3 <- all.vars(formula)[l+1]
      var1to3 <- c(var1,var2,var3)
      unique.vars <- unique(var1to3)
      var.dup <- var1to3[duplicated(var1to3)]
      if(length(unique.vars) == 2){
        vars <- mv_moment_3_4_bypair_df$vars
        row.ind <- grepl(unique.vars[1], vars, fixed = T) &
          grepl(unique.vars[2], vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
        if(unlist(gregexpr(var1, var.name)) == 1){
          mv_moment_3_4_bypair_df[row.ind, 'a3b']*n
        } else{mv_moment_3_4_bypair_df[row.ind, 'ab3']*n}
      } else{
        vars <- mv_moment_3_4_by3_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          grepl(var3, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        var.name <- mv_moment_3_4_by3_df$vars[row.ind]
        if(unlist(gregexpr(var1, var.name)) == 1){
          mv_moment_3_4_by3_df[row.ind, 'a2bc']*n
        } else if(unlist(gregexpr(var1, var.name)) +
                  nchar(var1) - 1 == nchar(var.name)){
          mv_moment_3_4_by3_df[row.ind, 'abc2']*n
        } else{mv_moment_3_4_by3_df[row.ind, 'ab2c']*n}
      }
    })
  }))
  
  x2_jXtX[lower.tri(x2_jXtX)] <- c(x2_jx_k, x2_jx_kx_l)
  x2_jXtX <- t(x2_jXtX)
  x2_jXtX[lower.tri(x2_jXtX)] <- c(x2_jx_k, x2_jx_kx_l)
  t(betas)%*% x2_jXtX %*% betas
})

diag(XtWX.3) <- c(bTxxTb, bTx2_jxxTb)

# Continue with the upper triangular part by transposing the lower triangular
bTx_jxxTb <- sapply(1:p, function(j){
  var1 <- all.vars(formula)[1+j]
  x_jXtX <- matrix(0, nrow = p + 1, ncol = p + 1) #initial matrix
  diag(x_jXtX) <- sapply(1:(p+1), function(k){
    if(k == 1){
      vars <- summary_stats$variable
      summary_stats[vars == var1, 'mean'] * n
    } else{
      var2 <- all.vars(formula)[k]
      if(var1 == var2){
        vars <- summary_stats$variable
        summary_stats[vars == var1, 'target_3_moment']*n
      } else{
        vars <- mv_moment_3_4_bypair_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
        if(unlist(gregexpr(var1, var.name)) == 1){
          mv_moment_3_4_bypair_df[row.ind, 'ab2']*n
        } else{mv_moment_3_4_bypair_df[row.ind, 'a2b']*n}
      }
    }
  })
  
  x_jx_k <- sapply(1:p, function(k){
    var2 <- all.vars(formula)[1+k]
    vars <- summary_stats$variable
    if(var1 == var2){
      summary_stats[vars == var1, 'variance'] * (n-1) +
        summary_stats[vars == var1, 'mean']^2 * n
    } else{
      var_cov_mat[var1, var2] * (n-1) +
        summary_stats[vars == var1, 'mean'] *
        summary_stats[vars == var2, 'mean'] * n
    }
  })
  x_jx_kx_l <- unlist(sapply(1:(p-1), function(k){
    var2 <- all.vars(formula)[k+1]
    sapply((k+1):p, function(l){
      var3 <- all.vars(formula)[l+1]
      var1to3 <- c(var1,var2,var3)
      unique.vars <- unique(var1to3)
      var.dup <- var1to3[duplicated(var1to3)]
      if(length(unique.vars) == 2){
        vars <- mv_moment_3_4_bypair_df$vars
        row.ind <- grepl(unique.vars[1], vars, fixed = T) &
          grepl(unique.vars[2], vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
        if(unlist(gregexpr(var1, var.name)) == 1){
          mv_moment_3_4_bypair_df[row.ind, 'a2b']*n
        } else{mv_moment_3_4_bypair_df[row.ind, 'ab2']*n}
      } else{
        vars <- mv_moment_3_4_by3_df$vars
        row.ind <- grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          grepl(var3, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        mv_moment_3_4_by3_df[row.ind, 'abc']*n
      }
    })
  }))
  
  x_jXtX[lower.tri(x_jXtX)] <- c(x_jx_k, x_jx_kx_l)
  x_jXtX <- t(x_jXtX)
  x_jXtX[lower.tri(x_jXtX)] <- c(x_jx_k, x_jx_kx_l)
  t(betas)%*% x_jXtX %*% betas
})


bTx_jx_kxxTb <- unlist(sapply(1:(p-1), function(m){
  var0 <- all.vars(formula)[1+m]
  sapply((m+1):p, function(j){
    var1 <- all.vars(formula)[1+j]
    x_jx_kXtX <- matrix(0, nrow = p + 1, ncol = p + 1) #initial matrix
    diag(x_jx_kXtX) <- sapply(1:(p+1), function(k){
      if(k == 1){
        vars <- summary_stats$variable
        var_cov_mat[var0, var1] * (n-1) +
          summary_stats[vars == var0, 'mean'] *
          summary_stats[vars == var1, 'mean'] * n
      } else{
        var2 <- all.vars(formula)[k]
        var0to2 <- c(var0,var1,var2)
        unique.vars <- unique(var0to2)
        var.dup <- var0to2[duplicated(var0to2)]
        if(length(unique.vars) == 2){
          vars <- mv_moment_3_4_bypair_df$vars
          row.ind <- grepl(unique.vars[1], vars, fixed = T) &
            grepl(unique.vars[2], vars, fixed = T) &
            !grepl(paste(c(":","X"), collapse = "|"), vars)
          var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
          if(unlist(gregexpr(var.dup, var.name)) == 1){
            mv_moment_3_4_bypair_df[row.ind, 'a3b'] * n
          } else{mv_moment_3_4_bypair_df[row.ind, 'ab3'] * n}
        } else{
          vars <- mv_moment_3_4_by3_df$vars
          row.ind <- grepl(var0, vars, fixed = T) &
            grepl(var1, vars, fixed = T) &
            grepl(var2, vars, fixed = T) &
            !grepl(paste(c(":","X"), collapse = "|"), vars)
          var.name <- mv_moment_3_4_by3_df$vars[row.ind]
          if(unlist(gregexpr(var2, var.name)) == 1){
            mv_moment_3_4_by3_df[row.ind, 'a2bc']*n
          } else if(unlist(gregexpr(var2, var.name)) +
                    nchar(var2) - 1 == nchar(var.name)){
            mv_moment_3_4_by3_df[row.ind, 'abc2']*n
          } else{mv_moment_3_4_by3_df[row.ind, 'ab2c']*n}
        }
      }
    })
    
    x_jx_kx_l <- sapply(1:p, function(k){
      var2 <- all.vars(formula)[k+1]
      var0to2 <- c(var0,var1,var2)
      unique.vars <- unique(var0to2)
      var.dup <- var0to2[duplicated(var0to2)]
      if(length(unique.vars) == 2){
        vars <- mv_moment_3_4_bypair_df$vars
        row.ind <- grepl(unique.vars[1], vars, fixed = T) &
          grepl(unique.vars[2], vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        var.name <- mv_moment_3_4_bypair_df$vars[row.ind]
        if(unlist(gregexpr(var2, var.name)) == 1){
          mv_moment_3_4_bypair_df[row.ind, 'a2b']*n
        } else{mv_moment_3_4_bypair_df[row.ind, 'ab2']*n}
      } else{
        vars <- mv_moment_3_4_by3_df$vars
        row.ind <- grepl(var0, vars, fixed = T) &
          grepl(var1, vars, fixed = T) &
          grepl(var2, vars, fixed = T) &
          !grepl(paste(c(":","X"), collapse = "|"), vars)
        mv_moment_3_4_by3_df[row.ind, 'abc']*n
      }
    })
    x_jx_kx_lx_m <- unlist(sapply(1:(p-1), function(k){
      var2 <- all.vars(formula)[k+1]
      sapply((k+1):p, function(l){
        var3 <- all.vars(formula)[l+1]
        var0to3 <- c(var0,var1,var2,var3)
        unique.vars <- unique(var0to3)
        var.dup <- var0to3[duplicated(var0to3)]
        if(length(unique.vars) == 3){
          vars <- mv_moment_3_4_by3_df$vars
          row.ind <- grepl(unique.vars[1], vars, fixed = T) &
            grepl(unique.vars[2], vars, fixed = T) &
            grepl(unique.vars[3], vars, fixed = T) &
            !grepl(paste(c(":","X"), collapse = "|"), vars)
          var.name <- mv_moment_3_4_by3_df$vars[row.ind]
          if(unlist(gregexpr(var.dup, var.name)) == 1){
            mv_moment_3_4_by3_df[row.ind, 'a2bc']*n
          } else if(unlist(gregexpr(var.dup, var.name)) +
                    nchar(var.dup) - 1 == nchar(var.name)){
            mv_moment_3_4_by3_df[row.ind, 'abc2']*n
          } else{mv_moment_3_4_by3_df[row.ind, 'ab2c']*n}
        } else if(length(unique.vars) == 2){
          vars <- mv_moment_3_4_bypair_df$vars
          row.ind <- grepl(unique.vars[1], vars, fixed = T) &
            grepl(unique.vars[2], vars, fixed = T) &
            !grepl(paste(c(":","X"), collapse = "|"), vars)
          mv_moment_3_4_bypair_df[row.ind, 'a2b2'] * n
        } else{
          vars <- mv_moment_4_df$vars
          row.ind <- grepl(var0, vars, fixed = T) &
            grepl(var1, vars, fixed = T) &
            grepl(var2, vars, fixed = T) &
            grepl(var3, vars, fixed = T) &
            !grepl(paste(c(":","X"), collapse = "|"), vars)
          mv_moment_4_df[row.ind, 'abcd']*n
        }
      })
    }))
    
    x_jx_kXtX[lower.tri(x_jx_kXtX)] <- c(x_jx_kx_l, x_jx_kx_lx_m)
    x_jx_kXtX <- t(x_jx_kXtX)
    x_jx_kXtX[lower.tri(x_jx_kXtX)] <- c(x_jx_kx_l, x_jx_kx_lx_m)
    t(betas)%*% x_jx_kXtX %*% betas
  })
}))



XtWX.3[lower.tri(XtWX.3)] <- c(bTx_jxxTb, bTx_jx_kxxTb)
XtWX.3 <- t(XtWX.3)
XtWX.3[lower.tri(XtWX.3)] <- c(bTx_jxxTb, bTx_jx_kxxTb)

# D. Add the three matrices ----------------------------
XtWX <- XtWX.1 - XtWX.2 + XtWX.3

# III. Compute robust variance matrix ----------------------------
robvarmat <- solve(XtX) %*% (XtWX) %*% solve(XtX)
sqrt(diag(robvarmat))



# IV. Compare with estimates from individual-level data ------------
lm_actual <- lm(formula, data = Xy.df)
coeftest(lm_actual, vcov = vcovHC(lm_actual, type = 'HC0'))
