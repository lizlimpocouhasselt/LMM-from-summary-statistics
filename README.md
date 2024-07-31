# LMM of federated data from summary statistics

## Description
This README demonstrates (via R) the concepts explained in our paper entitled [Linear mixed modelling of federated data when only the mean, covariance, and sample size are available](https://arxiv.org/abs/2407.20796), although these can be implemented in any other statistical softwares.

We present how a data analyst, upon receiving summary data from a data provider (e.g. clinic), can estimate a linear mixed model (LMM) without having access to individual observations for privacy reasons. For possible data preparation steps on the side of the data provider, code is also available in lmm_chop.R.

## Outline of demo code
I. Load packages / functions / [summary] data

II. Generate pseudo-data

III. Estimate LMM using pseudo-data

## I. Load packages / functions / [summary] data

The package we will use for estimating an LMM is 
```r
library(lme4)
```
while the following functions 
Here, we use the publicly available data from the Children's Hospital of Pennsylvania (CHOP) which can be accessed from the R package medicaldata.
```r
source('pseudo_data_gen_fn.R') #modified mvrnorm
source('pseudo_data_ls_fn.R') 
#other options are source('pseudo_data_gen_mvrnorm_option_fn.R') or source('pseudo_data_gen_chol_option_fn.R')
summary_stats <- readRDS("summary_stats_chop.R")
var_cov_mat <- readRDS("var_cov_mat_chop.R")
# Preview of summary statistics from 1 data provider
vars <- c('log_ct_result', 'gendermale', 'age', 'drive_thru_ind_factor1', 'gendermaleXage')
summary_stats[[26]][summary_stats[[26]][,'variable'] %in% vars, ]
round(var_cov_mat[[26]][vars,vars], 3)
``` 
