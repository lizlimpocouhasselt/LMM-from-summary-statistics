# LMM of federated data from summary statistics

## Description
This README demonstrates (via R) the concepts explained in our paper entitled [Linear mixed modelling of federated data when only the mean, covariance, and sample size are available](https://arxiv.org/abs/2407.20796), although these can be implemented in any other statistical softwares.

We present how a data analyst, upon receiving summary data from a data provider (e.g. clinic), can estimate a linear mixed model (LMM) without having access to individual observations for privacy reasons. For possible data preparation steps on the side of the data provider, code is also available in lmm_chop.R.

## Outline of demo code
I. Load packages / functions / [summary] data
II. Generate pseudo-data
III. Estimate LMM using pseudo-data

Run the following commands in R

```r
# Load packages
library(medicaldata)
library(skimr)
library(purrr)
library(dplyr)
``` 
