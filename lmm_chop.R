#---------------------------------------------------
# FEDERATED DATA ANALYSIS FOR LINEAR [MIXED] MODELS
# USING CHOP DATA
#--------------------------------------------------- 



### ------ PRELIMINARIES ------
# Clear memory
rm(list=ls(all=TRUE))

# Set working directory
#setwd("G:/My Drive/PhD/FDA")


### ------ DATA PROVIDER TASK ------

# Install packages
#install.packages('medicaldata')
#install.packages('skimr')
#install.packages('purrr')
#install.packages('dplyr')

## I. Load packages and data ------
library(medicaldata)
library(skimr)
library(purrr)
library(dplyr)
data("covid_testing")
skim(covid_testing) #overview of variables

## II. Preprocess data ------
# 1. include only valid values for ct_result
# 2. select variables for the analysis
# 3. convert from numeric to factor
data <- covid_testing %>%
  filter(ct_result > 0, result %in% c('positive', 'negative')) %>%
  dplyr::select(clinic_name, result, ct_result, gender, age, drive_thru_ind) %>%
  mutate(drive_thru_ind_factor = as.factor(drive_thru_ind), drive_thru_ind = NULL)

## III. Construct summary data table PER GROUP ------
# Include also summaries for interactions, log transformations, and quadratic terms

# A. Standardize numeric variablesm ------
data <- data %>%
  mutate(across(where(is.numeric), ~ scale(.)[, 1], .names = "std_{.col}"))

# B. Log-transform numeric variables ------
data <- data %>%
  mutate(across(where(is.numeric) &
                  !starts_with("std_"),
                ~ log(.x), 
                .names = "log_{.col}"))

# C. Compute quadratic terms for numeric variables ------
data <- data %>%
  mutate(across(where(is.numeric) &
                  !starts_with("std_") &
                  !starts_with("log_"),
                ~ (.x)^2, 
                .names = "{.col}_sq"))

# D. Compute interaction columns ------
# 1. Construct design matrix (as dataframe) from current dataset
data_design_df <- as.data.frame(model.matrix(~ ., data %>% dplyr::select(-c(clinic_name))))
data_design_df$clinic_name <- data$clinic_name
# 2. Break up data per group
group_data_design_df <- data_design_df %>%
  group_split(clinic_name) %>%
  keep(~ nrow(.x) > 1)
# 3. Extract unique group names
unique_group_names <- do.call(rbind, group_data_design_df) %>%
  pull(clinic_name) %>%
  unique()
# 4. Assign names to the list elements
names(group_data_design_df) <- unique_group_names
# 5. Create variable names for interaction terms
var_names_for_interaction <- combn(names(group_data_design_df[[1]] %>%
                                           dplyr::select(-starts_with('log') & -ends_with('_sq') & -'(Intercept)' & -starts_with('clinic'))), 2)
var_names_of_interaction <- apply(var_names_for_interaction, 2, 
                                  function(vars){paste0(vars[1],'X',vars[2])})
# 6. Include interaction columns in dataframe
group_data_design_df <- lapply(group_data_design_df,
                               function(grp){
                                 interaction_columns <- do.call(cbind,
                                                                apply(var_names_for_interaction, 2, function(var_name) grp[,var_name[1]]*grp[,var_name[2]]))
                                 names(interaction_columns) <- var_names_of_interaction
                                 cbind(grp, interaction_columns)
                               })

# E. Compute sufficient statistics per hospital ------
summary_stats <- lapply(group_data_design_df,
                        function(grp){
                          data.frame(
                            variable = names(grp[, c(-1)] %>% 
                                               dplyr::select(-(starts_with('clinic')))),
                            n = rep(nrow(grp), length(names(grp[, c(-1)] %>%
                                                                   dplyr::select(-(starts_with('clinic')))))),
                            mean = apply(grp[, c(-1)] %>%
                                           dplyr::select(-(starts_with('clinic'))),
                                         2, mean, na.rm = TRUE),
                            variance = ifelse(is.na(apply(grp[, c(-1)] %>%
                                               dplyr::select(-(starts_with('clinic'))),
                                             2, var)), 0, apply(grp[, c(-1)] %>%
                                                                  dplyr::select(-(starts_with('clinic'))),
                                                                2, var)),
                            row.names = NULL)
                          })
var_cov_mat <- lapply(group_data_design_df,
                      function(grp) cov(grp[, c(-1)] %>%
                                          dplyr::select(-(starts_with('clinic')))))
var_cov_mat <- lapply(var_cov_mat, function(grp) apply(grp, 2, function(column) ifelse(is.na(column),0,column)))

# Merge all dataframes into one
summary_stats_one <- do.call(rbind, summary_stats)
row.names(summary_stats_one) <- NULL
summary_stats_one$clinic_name <- rep(names(summary_stats), 
                                     each = nrow(summary_stats[[1]]))
var_cov_mat_one <- as.data.frame(do.call(rbind, var_cov_mat))
var_cov_mat_one$clinic_name <- rep(names(var_cov_mat), 
                                   each = nrow(var_cov_mat[[1]]))

# Save summary statistics and preprocessed data
data <- data %>% filter(clinic_name %in% names(group_data_design_df))
write.csv(summary_stats_one, file = "summary_stats_chop.csv")
write.csv(var_cov_mat_one, file = "var_cov_mat_chop.csv")
write.csv(data, file = "actual_data_preprocessed_chop.csv")



### ------ DATA ANALYST TASK ------
## I. Load packages / functions / [summary] data ------
#install.packages('lme4')
library(lme4)
source('pseudo_data_gen_fn.R') #modified mvrnorm
source('pseudo_data_ls_fn.R') 
#other options are source('pseudo_data_gen_mvrnorm_option_fn.R') or source('pseudo_data_gen_chol_option_fn.R')

id_summary_stats <- "1dw6y5CN1d1Kuh5PMyZ7Nh6_5kbyrpN4L"
summary_stats_one <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id_summary_stats))
summary_stats <- summary_stats_one %>% 
  dplyr::select(-1) %>% 
  split(f = as.factor(.$clinic_name)) %>%
  lapply(function(df){
    df[-ncol(df)]
  })

id_var_cov_mat <- "1dwuaHNbUAtnzV8FtBiH7H_isRqkHjT1w"
var_cov_mat_one <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id_var_cov_mat))
var_cov_mat <- var_cov_mat_one %>% 
  dplyr::select(-1) %>%
  split(f = as.factor(.$clinic_name)) %>%
  lapply(function(df) {
      rownames(df) <- colnames(df)[-ncol(df)]
    return(as.matrix(df[-ncol(df)]))
  })

# Preview of summary statistics from 1 data provider
vars <- c('log_ct_result', 'gendermale', 'age', 'drive_thru_ind_factor1', 'gendermaleXage')
summary_stats[[26]][summary_stats[[26]][,'variable'] %in% vars, ]
round(var_cov_mat[[26]][vars,vars], 3)

## II. Generate pseudo-data ------
# A. Provide model specification ------
formula <- log_ct_result ~ gendermale + std_age + drive_thru_ind_factor1 + gendermaleXstd_age 




# B. Generate pseudo-data ------
# (initial random numbers come from a std uniform distribution)
set.seed(121314)
pseudo_data <- lapply(names(summary_stats), pseudo_data_ls_fn)
names(pseudo_data) <- names(summary_stats)

# C. Combine [list] into one df ------
pooled_pseudo_data <- do.call(rbind, pseudo_data)

## III. Estimate LMM using pseudo-data
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



# ------
# Perform posthoc procedures via lmerTest
# lmm_pseudo <- lmerTest::lmer(log_ct_result ~ gendermale + age + drive_thru_ind_factor1 + gendermaleXage + (1|clinic_name), data = pooled_pseudo_data)
# summary(lmm_pseudo)
# lmm_actual <- lmerTest::lmer(log_ct_result ~ gender * std_age + drive_thru_ind_factor + (1|clinic_name), data = pooled_actual_data)
# summary(lmm_actual)
# 
# anova(lmm_pseudo, type = 2, ddf = "Kenward-Roger")
# anova(lmm_actual, type = 2, ddf = "Kenward-Roger")