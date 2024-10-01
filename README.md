# LMM of federated data from summary statistics

## Description
This README demonstrates (via R) the concepts explained in our paper entitled [Linear mixed modelling of federated data when only the mean, covariance, and sample size are available](https://arxiv.org/abs/2407.20796), although these can be implemented in any other statistical softwares.

We present how a data analyst, upon receiving summary data from a data provider (e.g. clinic), can estimate a linear mixed model (LMM) without having access to individual observations for privacy reasons. For possible data preparation steps on the side of the data provider, code is also available in lmm_chop.R.

## Strategy
I. Generate pseudo-data

II. Estimate LMM using pseudo-data

## Demo code
To generate pseudo-data, we utilize the following functions:
```r
source('pseudo_data_gen_fn.R') #modified mvrnorm
source('pseudo_data_ls_fn.R') 
``` 
We also need the following R packages: 
```r
library(lme4)
library(dplyr)
```
For this demo, we used the publicly available data from the Children's Hospital of Pennsylvania (CHOP) which can be accessed from the R package [medicaldata](https://cran.r-project.org/web/packages/medicaldata/medicaldata.pdf). It contains deidentified patient information and COVID-19 test results from 88 clinics, although we only utilized 70 clinics for this demo after filtering out incomplete observations, observations with invalid values, and clinics with only 1 patient record. From these, we computed summary statistics such as the mean vector, covariance matrix, and sample size per clinic (see DATA PROVIDER TASK portion of lmm_chop.R). To mimic a real-world scenario in this demo, we assume that the data analyst has received only these summary statistics: 
```r
# Load summary data
id_summary_stats <- "1dw6y5CN1d1Kuh5PMyZ7Nh6_5kbyrpN4L"
summary_stats_one <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id_summary_stats))
summary_stats <- summary_stats_one %>%        #breaks up data into 1 dataframe per clinic
  dplyr::select(-1) %>% 
  split(f = as.factor(.$clinic_name)) %>%
  lapply(function(df){
    df[-ncol(df)]
  })

id_var_cov_mat <- "1dwuaHNbUAtnzV8FtBiH7H_isRqkHjT1w"
var_cov_mat_one <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id_var_cov_mat))
var_cov_mat <- var_cov_mat_one %>%        #breaks up data into 1 dataframe per clinic 
  dplyr::select(-1) %>%
  split(f = as.factor(.$clinic_name)) %>%
  lapply(function(df) {
      rownames(df) <- colnames(df)[-ncol(df)]
    return(as.matrix(df[-ncol(df)]))
  })

# Show preview of summary statistics from 1 data provider
vars <- c('log_ct_result', 'gendermale', 'age', 'drive_thru_ind_factor1', 'gendermaleXage')
summary_stats[[26]][summary_stats[[26]][,'variable'] %in% vars, ]
round(var_cov_mat[[26]][vars,vars], 3)
``` 

Next, we specify the model formula we wish to use:
```r
# A. Provide model specification ------
formula <- log_ct_result ~ gendermale + std_age + drive_thru_ind_factor1 + gendermaleXstd_age 
```

Then we generate pseudo-data for the specified variables and combine them into 1 dataframe:
```r
# B. Generate pseudo-data ------
pseudo_data <- lapply(names(summary_stats), pseudo_data_ls_fn)
names(pseudo_data) <- names(summary_stats)

# C. Combine [list] into one df ------
pooled_pseudo_data <- do.call(rbind, pseudo_data)

```

With the pseudo-data, we can estimate an LMM with random-intercept only or with both random intercept and random slope and compare the output to that produced from actual data:

#### *Random intercept only*
```r
# A. With random-intercept only ------
lmm_pseudo <- lmer(log_ct_result ~ gendermale + std_age + drive_thru_ind_factor1 + gendermaleXstd_age + (1|clinic_name), data = pooled_pseudo_data)
summary(lmm_pseudo)

# Compare with LMM using actual data
id_data <- "1eSLj-vAwAOqlM1PV8WSdYl-fNcouWdme"
pooled_actual_data <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id_data))
lmm_actual <- lmer(log_ct_result ~ gender * std_age + drive_thru_ind_factor + (1|clinic_name), data = pooled_actual_data)
summary(lmm_actual)

# Compute AIC and confidence intervals (random intercept only)
AIC(lmm_pseudo)
AIC(lmm_actual)
BIC(lmm_pseudo)
BIC(lmm_actual)
confint(lmm_pseudo)
confint(lmm_actual)

```

#### *Random intercept + random slope*
```r
# B. With random intercept+slope) ------
lmm_pseudo <- lmer(log_ct_result ~ gendermale + std_age + drive_thru_ind_factor1 + gendermaleXstd_age + (1+std_age|clinic_name), data = pooled_pseudo_data)
summary(lmm_pseudo)

# Compare with LMM using actual data (random intercept+slope)
lmm_actual <- lmer(log_ct_result ~ gender * std_age + drive_thru_ind_factor + (1+std_age|clinic_name), data = pooled_actual_data)
summary(lmm_actual)


# Compute AIC and confidence intervals (random intercept+slope)
AIC(lmm_pseudo)
AIC(lmm_actual)
BIC(lmm_pseudo)
BIC(lmm_actual)
confint(lmm_pseudo)
confint(lmm_actual)

```

These output should be exactly the same except for the residuals, which would require the actual observations that are unavailable in practice.

