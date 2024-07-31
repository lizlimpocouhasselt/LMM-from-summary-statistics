#---------------------------------------------------
# FEDERATED DATA ANALYSIS FOR LINEAR [MIXED] MODELS
# USING CHOP DATA
#--------------------------------------------------- 



# ------ PRELIMINARIES ------
# Clear memory
rm(list=ls(all=TRUE))

# Set working directory
setwd("G:/My Drive/PhD/FDA")

# Install packages
#install.packages('dplyr')
#install.packages('purrr')
#install.packages('lubridate')
#install.packages('skimr')
#install.packages('MASS')
#install.packages('lme4')
#install.packages('medicaldata')
#install.packages('lmerTest')


# Load packages
library(medicaldata)
library(skimr)
library(purrr)
library(dplyr)


library(lubridate)
library(MASS)
library(lme4)
library(lmerTest)





# ------ DATA PROVIDER TASK ------

# Load data
data("covid_testing")
skim(covid_testing) #overview of variables



# Preprocess data
# 1. include only valid values for ct_result
# 2. select variables for the analysis
# 3. convert from numeric to factor
data <- covid_testing %>%
  filter(ct_result > 0, result %in% c('positive', 'negative')) %>%
  dplyr::select(clinic_name, result, ct_result, gender, age, drive_thru_ind) %>%
  mutate(drive_thru_ind_factor = as.factor(drive_thru_ind), drive_thru_ind = NULL)



# Construct summary data table PER GROUP
# Include also summaries for interactions, log transformations, and quadratic terms


# Standardize numeric variables
data <- data %>%
  mutate(across(where(is.numeric), ~ scale(.)[, 1], .names = "std_{.col}"))


# Log-transform numeric variables
data <- data %>%
  mutate(across(where(is.numeric) &
                  !starts_with("std_"),
                ~ log(.x), 
                .names = "log_{.col}"))


# Compute quadratic terms for numeric variables
data <- data %>%
  mutate(across(where(is.numeric) &
                  !starts_with("std_") &
                  !starts_with("log_"),
                ~ (.x)^2, 
                .names = "{.col}_sq"))


# Construct design matrix (as dataframe) from current dataset
data_design_df <- as.data.frame(model.matrix(~ ., data %>% dplyr::select(-c(clinic_name))))
data_design_df$clinic_name <- data$clinic_name

# Break up data per group
group_data_design_df <- data_design_df %>%
  group_split(clinic_name) %>%
  keep(~ nrow(.x) > 1)

# Extract unique group names
unique_group_names <- do.call(rbind, group_data_design_df) %>%
  pull(clinic_name) %>%
  unique()

# Assign names to the list elements
names(group_data_design_df) <- unique_group_names



# Compute interactions per group
var_names_for_interaction <- combn(names(group_data_design_df[[1]] %>%
                                           dplyr::select(-starts_with('log') & -ends_with('_sq') & -'(Intercept)' & -starts_with('clinic')
                                                         )), 2)
var_names_of_interaction <- apply(var_names_for_interaction, 2, 
                                  function(vars){paste0(vars[1],'X',vars[2])})
group_data_design_df <- lapply(group_data_design_df,
                               function(grp){
                                 interaction_columns <- do.call(cbind,
                                                                apply(var_names_for_interaction, 2, function(var_name) grp[,var_name[1]]*grp[,var_name[2]]))
                                 names(interaction_columns) <- var_names_of_interaction
                                 cbind(grp, interaction_columns)
                               })



# Compute sufficient statistics per hospital
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



# ------ DATA ANALYST TASK ------
# Given the following summary statistics from each data provider
vars <- c('log_ct_result', 'gendermale', 'age', 'drive_thru_ind_factor1', 'gendermaleXage')
summary_stats[[26]][summary_stats[[26]][,'variable'] %in% vars, ]
round(var_cov_mat[[26]][vars,vars], 3)

#---------------------------------------------
# Generate pseudo-data

# Provide formula
formula <- log_ct_result ~ gendermale + std_age + drive_thru_ind_factor1 + gendermaleXstd_age 



# Use mvrnorm to generate pseudo-data (initial random numbers come from a std normal distribution)
#set.seed(121314)
#pseudo_data <- lapply(names(summary_stats),
#                      function(grp_name){
#                        pseudo_ <- as.data.frame(
#                          mvrnorm(summary_stats[[grp_name]][1,'n'],
#                                  mu = summary_stats[[grp_name]][match(all.vars(formula), summary_stats[[grp_name]]$variable), 'mean'],
#                                  Sigma = var_cov_mat[[grp_name]][all.vars(formula), all.vars(formula)],
#                                  empirical = TRUE
#                          )
#                        )
#                        pseudo_$clinic_name <- grp_name
#                        pseudo_
#                      })





# Use modified mvrnorm to generate pseudo-data (initial random numbers come from a std uniform distribution)
source('pseudo_data_gen.R')
set.seed(121314)
pseudo_data <- lapply(names(summary_stats),
                      function(grp_name){
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
                      })
names(pseudo_data) <- names(summary_stats)





# Use manual method via Cholesky decomposition
#set.seed(121314)
#pseudo_data <- lapply(names(summary_stats),
#                      function(grp_name){
                        # Generate numbers randomly from uniform distribution
#                        n <- summary_stats[[grp_name]][1,'n']
#                        p <- length(all.vars(formula))
#                        desired_mean <- summary_stats[[grp_name]][summary_stats[[grp_name]][,'variable'] %in% all.vars(formula), 'mean']
#                        desired_cov_mat <- var_cov_mat[[grp_name]][all.vars(formula), all.vars(formula)]
#                        r <- matrix(rnorm(n * p), ncol = p)
                        # Compute mean vector and covariance matrix of generated numbers
#                        r_mean <- apply(r, 2, mean)
#                        r_cov <- cov(r)
                        # Compute centered values of r 
#                        r_star <- r - matrix(rep(r_mean, rep(n,p)), ncol = p)
                        # Compute scaled value of r
#                        r_0 <- solve(t(chol(r_cov))) %*% t(r_star)
                        # Generate dataset with desired mean vector and covariance matrix
#                        ps <- t(chol(desired_cov_mat)) %*% r_0 + t(matrix(rep(desired_mean, rep(n,p)), ncol = p))
#                        pseudo_ <- as.data.frame(t(ps))
#                        names(pseudo_) <- all.vars(formula)
#                        pseudo_$clinic_name <- grp_name
#                        pseudo_
#                      })


# Combine list into one df
pooled_pseudo_data <- pseudo_data[[1]]
for(h in 2:length(summary_stats)){
  pooled_pseudo_data <- rbind(pooled_pseudo_data, pseudo_data[[h]])
}

# Estimate LMM using pseudo-data (random intercept only)
lmm_pseudo <- lmer(log_ct_result ~ gendermale + std_age + drive_thru_ind_factor1 + gendermaleXstd_age + (1|clinic_name), data = pooled_pseudo_data)
summary(lmm_pseudo)

# Compare with LMM using actual data (random intercept only)
pooled_actual_data <- data %>% filter(clinic_name %in% names(group_data_design_df))
lmm_actual <- lmer(log_ct_result ~ gender * std_age + drive_thru_ind_factor + (1|clinic_name), data = pooled_actual_data)
summary(lmm_actual)


# Compute AIC and confidence intervals (random intercept only)
AIC(lmm_pseudo)
AIC(lmm_actual)
BIC(lmm_pseudo)
BIC(lmm_actual)
confint(lmm_pseudo)
confint(lmm_actual)


# Perform posthoc procedures via lmerTest
# lmm_pseudo <- lmerTest::lmer(log_ct_result ~ gendermale + age + drive_thru_ind_factor1 + gendermaleXage + (1|clinic_name), data = pooled_pseudo_data)
# summary(lmm_pseudo)
# lmm_actual <- lmerTest::lmer(log_ct_result ~ gender * std_age + drive_thru_ind_factor + (1|clinic_name), data = pooled_actual_data)
# summary(lmm_actual)
# 
# anova(lmm_pseudo, type = 2, ddf = "Kenward-Roger")
# anova(lmm_actual, type = 2, ddf = "Kenward-Roger")


# Estimate LMM using pseudo-data (random intercept+slope)
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


