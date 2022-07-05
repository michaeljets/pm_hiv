source('utils.R')

library(survival)
library(readr)
library(dplyr)
library(grf, lib.loc='../RLibs')
require(WeightSVM, lib.loc='../RLibs')
require(modelObj, lib.loc='../RLibs')
require(DynTxRegime, lib.loc='../RLibs')


# read in data
data = readRDS('data.rds')

# covariate matrix
covs = paste0('X', 1:5)
X = data %>% select_at(covs) %>% as.matrix()

# set folds and repetitions
M = 1  # number of repetitions
K = 5  # number of folds


# Estimate propensity score and censoring model ---------------------------

# estimate propensity scores
prop_frmla = as.formula(paste0('treat ~', paste0(covs, collapse='+')))
prop_model = glm(prop_frmla, data = data, family = binomial)
prop_scores = predict(prop_model, type = 'response')

# censoring weights
frmla = paste0('Surv(obs_times, ltfu) ~ treat +', paste0(covs, collapse='+'))
censor_frmla = as.formula(frmla)
cens_model = coxph(censor_frmla, data = data)
cens_weights = survfit(cens_model, newdata = data)$surv %>% t()


# Ridge -------------------------------------------------------------------

ridge = cv_estimator(data, covs=covs, method='cox', nfolds=K, reps=M, extra_args='ridge')


# Lasso -------------------------------------------------------------------

lasso = cv_estimator(data, covs=covs, method='cox', nfolds=K, reps=M, extra_args='lasso')


# Elastic net -------------------------------------------------------------

elastic = cv_estimator(data, covs=covs, method='cox', nfolds=K, reps=M, extra_args='elastic')


# CSF ---------------------------------------------------------------------

csf = cv_estimator(data, covs=covs, method='csf', nfolds=K, reps=M, prop_scores=prop_scores, censor=cens_weights)


# Genetic algorithm -------------------------------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# run on whole data
whole_genetic = get_genetic_dtr(data=dat, smooth=FALSE, augmented=FALSE, covs=covs,
                                prop_scores=prop_scores, censor=cens_weights,
                                max.generations=50)

# set initial values
init_vals = whole_genetic$eta %>% round(0)

# estimate DTR
genetic = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M, 
                       prop_scores=prop_scores, censor=cens_weights,
                       extra_args=list(init_vals=init_vals, max.generations=50))



# Genetic algorithm, augmented --------------------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# run on whole data
whole_genetic_a = get_genetic_dtr(data=dat, smooth=FALSE, augmented=TRUE, covs=covs,
                                  prop_scores=prop_scores, censor=cens_weights,
                                  max.generations=50)

# set initial values
init_vals_a = whole_genetic_a$eta %>% round(0)

# estimate DTR
genetic_a = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M, 
                         prop_scores=prop_scores, censor=cens_weights,
                         extra_args=list(init_vals=init_vals_a, augmented=TRUE, max.generations=50))



# Genetic algorithm, smooth -----------------------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# run on whole data
whole_genetic_s = get_genetic_dtr(data=dat, smooth=TRUE, augmented=FALSE, covs=covs,
                                  prop_scores=prop_scores, censor=cens_weights,
                                  max.generations=50)

# set initial values
init_vals_s = whole_genetic_s$eta %>% round(0)

# estimate DTR
genetic_s = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M, 
                         prop_scores=prop_scores, censor=cens_weights,
                         extra_args=list(init_vals=init_vals_s, smooth=TRUE, max.generations=50))


# Genetic algorithm, augmented + smooth -----------------------------------

# format data
dat = cbind(treat=data$treat, X, delta=data$delta, obs_times=data$obs_times) %>% as.data.frame()

# run on whole data
whole_genetic_as = get_genetic_dtr(data=dat, smooth=TRUE, augmented=TRUE, covs=covs,
                                   prop_scores=prop_scores, censor=cens_weights, 
                                   max.generations=50)

# set initial values
init_vals_as = whole_genetic_as$eta %>% round(0)

# estimate DTR
genetic_as = cv_estimator(data, covs=covs, method='genetic', nfolds=K, reps=M, 
                          prop_scores=prop_scores, censor=cens_weights,
                          extra_args=list(init_vals=init_vals_as, smooth=TRUE, augmented=TRUE, max.generations=50))


# RIST --------------------------------------------------------------------

# impute censored values using RIST
source('RISTfunctions.r')

# set up data
frmla = as.formula(paste0('~-1+treat*(', paste0(covs, collapse='+'), ')') )
rist_data = model.matrix(frmla, data = data)
rist_data = as.data.frame(cbind(rist_data, delta=data$r_delta, obs_times=data$obs_times))

# parameters
P = ncol(rist_data)-2   # number of dimension for X
K2 = 3  # number of covariate considered per spilt
nmin = 10   # minimum number of observed data in each node
M2 = 250     # number of trees in each fold
L = 2    # number of folds
tao = 60    # length of study

# split dataset into censored
rist_data_c = rist_data[data$r_delta==0, ]

# RIST
R_Muti_ERT_build = Muti_ERT_fit(as.matrix(rist_data), M2, K2, L, nmin, SupLogRank=1, tao=tao, impute="random")
R_Muti_ERT_predict= Muti_ERT_Predict(as.matrix(rist_data_c), R_Muti_ERT_build$Forest_seq[[L]], R_Muti_ERT_build$SurvMat_seq[[L]], R_Muti_ERT_build$time_intrest)

# predictions
times_pred_c = R_Muti_ERT_predict$TMT_predict
times_pred = data$obs_times
times_pred[data$r_delta==0] = times_pred_c

# add imputed survival times (i.e. now complete case) to dataset
data$rist = times_pred


# OWL ---------------------------------------------------------------------

# estimate DTR
owl_rist = cv_estimator(data, method='owl', nfolds=K, reps=M, covs=covs,
                        prop_scores=prop_scores, extra_args=list('rist', kernel='linear'))


# RWL ---------------------------------------------------------------------

rwl = cv_estimator(data, method='rwl', covs=covs, nfolds=K, reps=M)


# Save results ------------------------------------------------------------

# save treatment rules
saveRDS(mget(c('ridge', 'lasso', 'elastic', 'csf', 
               'genetic', 'genetic_a', 'genetic_s', 'genetic_as',
               'owl_rist', 'rwl')),
        file = 'results/cv_treatrules.rds')

# save imputed survival times from RIST
saveRDS(times_pred, 'results/rist_imputation.rds')

# save initial values for genetic, aug + smooth
saveRDS(init_vals_as, 'results/init_vals_as.rds')
