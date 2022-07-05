# computes bootstraps for standard error estimation of
# cross-validated risk curves

source('utils.R')

library(survival)
library(readr)
library(dplyr)


# read in data
data.orig = readRDS('data.rds')

# get initial values
init_vals_as = readRDS('results/init_vals_as.rds')

# specify covariates
covs = paste0('X', 1:5)

# set folds and repetitions
M = 1  # number of repetitions
K = 5  # number of folds

# set number of bootstraps
B = 200

# create matrices to store results
boot_estimates_0 = matrix(NA, nrow=B, ncol=60)
boot_estimates_1 = matrix(NA, nrow=B, ncol=60)
boot_estimates_dtr = matrix(NA, nrow=B, ncol=60)

for (b in 1:B) {
  
  # take bootstrap resample of data
  boot_ind = sample(1:nrow(data.orig), size=nrow(data.orig), replace=TRUE)
  data = data.orig %>%
    slice(boot_ind) %>%
    mutate(id = 1:n()) %>%
    arrange(obs_times)
  
  # estimate propensity scores
  prop_frmla = as.formula(paste0('treat ~', paste0(covs, collapse='+')))
  prop_model = glm(prop_frmla, data = data, family = binomial)
  prop_scores_boot = predict(prop_model, type = 'response')
  
  # censoring weights
  frmla = paste0('Surv(obs_times, ltfu) ~ treat +', paste0(covs, collapse='+'))
  censor_frmla = as.formula(frmla)
  cens_model = coxph(censor_frmla, data = data)
  cens_weights_boot = survfit(cens_model, newdata = data)$surv %>% t()
  
  # non-dtr curves
  curv_0_boot = 1-ipwe_km(data$obs_times, data$treat, data$delta, prop_scores_boot, cens_weights_boot, 0)
  curv_1_boot = 1-ipwe_km(data$obs_times, data$treat, data$delta, prop_scores_boot, cens_weights_boot, 1)
  
  # estimate DTR
  genetic_as = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M, 
                            prop_scores=prop_scores_boot, censor=cens_weights_boot,
                            extra_args=list(init_vals=init_vals_as, smooth=TRUE, augmented=TRUE, max.generations=50))
  
  # get risk estimate under chosen treatment rule
  curv_d_boot = extract_risk(genetic_as$treat_rule, data, prop_scores_boot, cens_weights_boot)
  
  # some bootstraps may not have all 60 unique times
  # but want risk estimates for all 60 unique times
  # solve by carrying forward previous estimates
  if (length(unique(data$obs_times)) < length(unique(data.orig$obs_times))) {
    which_miss = setdiff(unique(data.orig$obs_times), unique(data$obs_times))
    for (miss in which_miss) {
      if (miss==1) {
        curv_0_boot = c(0, curv_0_boot)
        curv_1_boot = c(0, curv_1_boot)
        curv_d_boot = c(0, curv_d_boot)
      }
      else {
        curv_0_boot = append(curv_0_boot, curv_0_boot[miss-1], after=miss-1)
        curv_1_boot = append(curv_1_boot, curv_1_boot[miss-1], after=miss-1)
        curv_d_boot = append(curv_d_boot, curv_d_boot[miss-1], after=miss-1)
      }
    }
  }
  
  # store estimates
  boot_estimates_0[b,] = curv_0_boot
  boot_estimates_1[b,] = curv_1_boot
  boot_estimates_dtr[b,] = curv_d_boot
  
}

# get estimate standard errors
boot_sd_0 = apply(boot_estimates_0, 2, sd)
boot_sd_1 = apply(boot_estimates_1, 2, sd)
boot_sd_dtr = apply(boot_estimates_dtr, 2, sd)


# Save results ------------------------------------------------------------

saveRDS(mget(c('boot_sd_0', 'boot_sd_1', 'boot_sd_dtr')),
        'results/se_estimates.rds')