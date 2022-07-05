

# Simulation functions ----------------------------------------------------

# functions for simulating data


GenSimTimes <- function(linear.part, lambda, gamma, censor=NULL, coarsen=TRUE, end=60, seed=NULL) {
  
  # simulates survival times using a linear function of data
  # generates data using a gompertz distribution with user-given parameters
  
  # set seed
  if (!is.null(seed)) (set.seed(seed))
  
  # generate true survival times
  true_times = (1/gamma) * log(((-gamma*log(runif(length(linear.part),0,1)))/(lambda*exp(linear.part))) + 1)
  
  # observed times
  obs_times = pmin(true_times, censor)
  
  # coarsen
  if (coarsen) {
    delta = as.numeric(ceiling(true_times) <= ceiling(censor))
    obs_times = ceiling(obs_times)
  }
  else {
    delta = as.numeric(true_times <= censor)
  }
  
  # return
  return(list(true_times = true_times, delta = delta, obs_times = obs_times))
}


GenSimCens <- function(linear.part, lambda, gamma, coarsen=TRUE, end=60, seed=NULL) {
  
  # simulates censoring times using a linear function of data
  # generates data using a gompertz distribution with user-given parameters
  
  # set seed
  if (!is.null(seed)) (set.seed(seed))
  
  # generate true survival times
  censor_times = (1/gamma) * log(((-gamma*log(runif(length(linear.part),0,1)))/(lambda*exp(linear.part))) + 1)
  
  # distinguish between LTFU and administrative censor
  if (coarsen) (ltfu = as.numeric(ceiling(censor_times) <= end))
  else (ltfu = as.numeric(censor_times <= end))
  
  # observed censored times
  obs_censor = pmin(censor_times, end)
  
  return(list(obs_censor=obs_censor, ltfu=ltfu, true_censor=censor_times))
}


GenSurvCurv <- function(linear.part, lambda, gamma, tpt=1:60) {
  
  # generate true survival probabilities
  sapply(tpt, function(x) exp( (-lambda*(exp(gamma*x)-1) / gamma) * exp(linear.part) ) )
}


# Dynamic treatment rules -------------------------------------------------

# functions in this section take in data and output treatment rules
# the output is a list where `treat_rule` is the treatment recommendations
# and `eta` are the coefficients to the treatment rule, if applicable


get_cox_dtr <- function(data, frmla, lambda='1se', penalty='lasso', standardize=TRUE, tpt=60, test_data=NULL) {
  
  require(survival)
  require(glmnet)
  require(dplyr)
  
  # set alpha by penalty
  if (penalty=='lasso') (alpha=1)
  else if (penalty=='elastic') (alpha=0.5)
  else if (penalty=='ridge') (alpha=0)
  else (stop('invalid penalty'))
  
  # covariate matrix
  X = model.matrix(frmla, data = data)
  X = X[,-1]  # remove intercept
  
  # Cox model
  y = Surv(data$obs_times, data$delta)
  cox.fit.reg = glmnet(X, y, family = 'cox', alpha = alpha, standardize = standardize)
  
  # cv to tune lambda
  cox.cvfit.reg = cv.glmnet(X, y, family = 'cox', type.measure = 'C', alpha = alpha, standardize = standardize)
  if (lambda=='1se') (lambda = cox.cvfit.reg$lambda.1se)
  else if (lambda=='min') (lambda = cox.cvfit.reg$lambda.min)
  else (stop('invalid lambda'))
  
  # get coefficients
  coefs = as.numeric(coef(cox.fit.reg, s=lambda))
  names(coefs) = row.names(coef(cox.fit.reg, s=lambda))
  
  # if input test data then get treatment recommendations on test data
  if (!is.null(test_data)) (data = test_data)
  
  # covariate matrices for ZOM
  X0 = model.matrix(frmla, data = data %>% mutate(treat=0))
  X0 = X0[,-1]  # remove intercept
  
  X1 = model.matrix(frmla, data = data %>% mutate(treat=1))
  X1 = X1[,-1]  # remove intercept
  
  # survival curves under ZOMs
  tmp0 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X0)
  tmp1 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X1)
  cox.surv0 = summary(tmp0, times = 1:tpt, extend = TRUE)$surv
  cox.surv1 = summary(tmp1, times = 1:tpt, extend = TRUE)$surv
  
  # get treatment rule
  if (is.null(nrow(cox.surv0))) {
    treat_rule = as.numeric(cox.surv1[tpt] >= cox.surv0[tpt])
  }
  else {
    treat_rule = as.numeric(cox.surv1[tpt,] >= cox.surv0[tpt,])
  }
  
  # return treatment rule and coefficients
  return(list(treat_rule=treat_rule, eta=coefs))
  
}


get_csf_dtr <- function(data, X, prop_scores, censor, tpt=1:60, test_data=NULL) {
  
  require(grf)
  
  # propensity score
  prop_scores = data$treat*prop_scores + (1-data$treat)*(1-prop_scores)
  
  # train CSF
  cs.forest = causal_survival_forest(X, data$obs_times, data$treat, data$r_delta,
                                     W.hat = prop_scores, C.hat = censor,
                                     mtry = floor(ncol(X)/3), num.trees = 5000, 
                                     failure.times = tpt, honesty=TRUE)
  
  # if input test data then get treatment recommendations on test data
  if (!is.null(test_data)) (X = test_data)
  
  # predict
  cs.pred = predict(cs.forest, newdata = X, estimate.variance = FALSE)
  
  # get treatment assignment
  treat_rule = as.numeric(cs.pred$predictions >= 0)
  
  return(list(treat_rule=treat_rule, eta=NA))
  
}


get_owl_dtr <- function(data, covs, prop_scores, kernel='linear', 
                        costs=NULL, gamma=NULL, tpt=60, test_data=NULL) {
  
  require(survival)
  require(WeightSVM)
  require(modelObj)
  require(DynTxRegime)
  require(dplyr)
  
  # svm formula
  svm_frmla = as.formula(paste0('~ ', paste0(covs, collapse='+')))
  
  # costs to consider
  if (is.null(costs)) (costs = 2^(seq(-8, 2, by=4)))
  else (costs = costs)
  
  # gamma to consider
  if (is.null(gamma) & kernel=='radial') (gamma = c(1/nrow(data), 2^(-5), 2^(-10)))
  else (gamma = gamma)
  
  # propensity score
  prop_scores = data$treat*prop_scores + (1-data$treat)*(1-prop_scores)
  
  # weights
  times_pred = data$rist
  w = times_pred / prop_scores
  
  # prep work
  P = ncol(data)-2
  Y = 2*data$treat - 1
  X = model.matrix(svm_frmla, data = data)
  X = scale(X[,-1])
  
  # svm
  if (length(costs)>1 | length(gamma)>1) {
        
    if (!is.null(gamma)) (search.grid = list(cost = costs, gamma = gamma))
    else (search.grid = list(cost = costs))
    
    svm_model = best.tune_wsvm(train.x=X, train.y=as.factor(data$treat), weight=w, ranges=search.grid,
                               tunecontrol=tune.control(sampling='cross', cross=3), 
                               kernel=kernel, scale=FALSE)
  }
  
  else {
    svm_model = wsvm(x=X, y=as.factor(data$treat), weight=w,
                     kernel=kernel, scale=FALSE, cost=costs, gamma=gamma)
  }
  
  # treatment rule
  treat_rule = predict(svm_model, newdata=X) %>% as.character() %>% as.numeric()
  
  # get coefficients of decision rule
  beta = drop(t(svm_model$coefs) %*% svm_model$SV)
  beta0 = -svm_model$rho
  eta = c(beta0, beta)
  
  # apply to test set
  if (!is.null(test_data)) {
    X_test = model.matrix(svm_frmla, data = test_data)
    X_test = scale(X_test[,-1])
    treat_rule = predict(svm_model, newdata=X_test) %>% as.character() %>% as.numeric()
    return(list(treat_rule=treat_rule, eta=eta))
  }
  
  return(list(treat_rule=treat_rule, eta=eta, svm_model=svm_model))
  
}


get_rwl_dtr <- function(data, covs, costs=NULL, tpt=60, test_data=NULL) {
  
  require(survival)
  require(WeightSVM)
  require(modelObj)
  require(DynTxRegime)
  require(dplyr)
  
  # formula
  frmla = as.formula(paste0('~ ', paste0(covs, collapse='+')))
  
  # costs to consider
  if (is.null(costs)) (costs = 2^(seq(-8, 2, by=2)))
  else (costs = costs)
  
  # reward
  times_pred = data$rist
  
  # define regression models
  moPropen <- modelObj::buildModelObj(model = frmla, 
                                      solver.method = 'glm', 
                                      solver.args = list(family='binomial'), 
                                      predict.method = 'predict.glm', 
                                      predict.args = list(type='response'))
  
  moMain <- modelObj::buildModelObj(model = frmla,
                                    solver.method = 'lm')
  
  # specify class of treatment regimes (linear)
  kernel = 'linear'
  
  # scale data
  data = data %>%
    mutate_at(covs, function(x) as.numeric(scale(x)))
  
  data = as.data.frame(data)
  
  if (length(costs)>1) {
    svm_model <- rwl(moPropen = moPropen,
                     moMain = moMain,
                     data = data,
                     response = times_pred,
                     txName = 'treat',
                     regime = frmla,
                     lambdas = costs,
                     cvFolds = 5L,
                     kernel = kernel,
                     responseType = 'continuous',
                     verbose = 1L)
  }
  
  else {
    svm_model <- rwl(moPropen = moPropen,
                     moMain = moMain,
                     data = data,
                     response = times_pred,
                     txName = 'treat',
                     regime = frmla,
                     lambdas = costs,
                     cvFolds = 0L,
                     kernel = kernel,
                     responseType = 'continuous',
                     verbose = 1L)
  }
  
  # apply to test set
  if (!is.null(test_data)) {
    eta = regimeCoef(object = svm_model)
    X_test = model.matrix(frmla, data = test_data)
    X_test = scale(X_test[,-1])
    treat_rule = as.numeric(cbind(1, X_test) %*% eta >= 0)
    return(list(treat_rule=treat_rule, eta=eta))
  }
  
  # return treatment rule
  treat_rule = optTx(x=svm_model, newdata=data)$optimalTx
  
  # get coefficients
  eta = regimeCoef(object = svm_model)
  
  return(list(treat_rule=treat_rule, eta=eta, svm_model=svm_model))
}


get_genetic_dtr <- function(data, prop_scores, censor, covs, smooth=FALSE, augmented=FALSE, 
                            standardize=TRUE, max.generations=100, pop.size=1000, 
                            tpt=60, init_vals=NULL, test_data=NULL) {
  
  # data must be ordered X (treat first column), delta, times
  
  library(compiler)
  library(rgenoud)
  library(survival)
  library(glmnet)
  require(dplyr)


  # make sure data is sorted
  stopifnot(!is.unsorted(data$obs_times))
  
  # X2 is treatment regime covariates
  X2 = as.matrix(data[,2:(ncol(data)-2)])
  
  
  # augmented
  if (augmented==TRUE) {
    
    # failure model
    cox_frmla = as.formula(paste0('~ treat * (', paste0(covs, collapse='+'), ')'))
    
    X = model.matrix(cox_frmla, data = data)
    X = X[,-1]  # remove intercept
    
    # penalized Cox model
    y = Surv(data$obs_times, data$delta)
    cox.fit.reg = glmnet(X, y, family = 'cox', alpha = 0.5)

    # cv to tune lambda
    cox.cvfit.reg = cv.glmnet(X, y, family = 'cox', type.measure = 'C', alpha = 0.5)
    lambda = cox.cvfit.reg$lambda.min

    # covariate matrices
    X0 = model.matrix(cox_frmla, data = data %>% mutate(treat=0))
    X0 = X0[,-1]  # remove intercept

    X1 = model.matrix(cox_frmla, data = data %>% mutate(treat=1))
    X1 = X1[,-1]  # remove intercept

    # surv fit
    failure.fit.0 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X0)
    failure.surv.0 = summary(failure.fit.0, times = 1:tpt, extend = TRUE)$surv %>% t()

    failure.fit.1 = survfit(cox.fit.reg, s = lambda, x = X, y = y, newx = X1)
    failure.surv.1 = summary(failure.fit.1, times = 1:tpt, extend = TRUE)$surv %>% t()
    
    
    # failure model hazard
    failure.cumhaz = summary(failure.fit.0, times = 1:tpt, extend = TRUE)$cumhaz %>% t()
    failure.haz0 = failure.cumhaz[,1]
    failure.haz = apply(failure.cumhaz, 1, diff) %>% t()
    failure.haz.0 = cbind(failure.haz0, failure.haz)
    
    failure.cumhaz = summary(failure.fit.1, times = 1:tpt, extend = TRUE)$cumhaz %>% t()
    failure.haz0 = failure.cumhaz[,1]
    failure.haz = apply(failure.cumhaz, 1, diff) %>% t()
    failure.haz.1 = cbind(failure.haz0, failure.haz)
    
    aug_terms = list(failure.surv.0=failure.surv.0,
                     failure.surv.1=failure.surv.1,
                     lambda.surv.0=failure.haz.0,
                     lambda.surv.1=failure.haz.1)
    
    # use algorithm to get eta
    eta = Genetic.IPWE(fn=aipwe_km_opt, X=X2, u=data$obs_times, z=data$treat, delta=data$delta,
                       prop_scores=prop_scores, censor=censor, smooth=smooth, 
                       nvars=ncol(X2)+1, max.generations=max.generations, pop.size=pop.size,
                       aug_terms=aug_terms, init_vals=init_vals, standardize=standardize)
    val_hat = eta[length(eta)]
    eta = eta[1:(length(eta)-1)]
    names(eta) = c('intercept', colnames(X2))
    
  }
  
  else {
    
    # use algorithm to get eta
    eta = Genetic.IPWE(fn=ipwe_km_opt, X=X2, u=data$obs_times, z=data$treat, delta=data$delta,
                       prop_scores=prop_scores, censor=censor, smooth=smooth, 
                       nvars=ncol(X2)+1, max.generations=max.generations, pop.size=pop.size,
                       init_vals=init_vals, standardize=standardize)
    val_hat = eta[length(eta)]
    eta = eta[1:(length(eta)-1)]
    names(eta) = c('intercept', colnames(X2))
    
  }
  
  # treatment rule
  if (standardize) (X2 = scale(X2))
  treat_rule = as.numeric(cbind(1,X2) %*% eta >= 0)
  
  if (!is.null(test_data)) {
    test_data = test_data[, !grepl('drugstart', colnames(test_data))]
    if (standardize) (test_data = scale(test_data))
    treat_rule = as.numeric(cbind(1,test_data) %*% eta >= 0)
    return(list(treat_rule=treat_rule, eta=eta))
  }
  
  return(list(treat_rule=treat_rule, eta=eta, val_hat=val_hat))
}


# CV estimator ------------------------------------------------------------


cv_estimator <- function(data, method, covs, prop_scores=NULL, censor=NULL, nfolds=NULL, reps=NULL, extra_args=NULL) {
  
  # computes CV treatment rules
  # see `get_treat_rules.R` for examples of usage 
  
  # Inputs:
  #   - data: data frame with survival times (`obs_times`), event indicator (`delta`), 
  #           treatment indicator (`treat`), and covariates (specified in `covs` argument)
  #   - method: character vector of length 1, which precision medicine method? options:
  #         - "cox": Cox regression
  #         - "genetic": genetic algorithm
  #         - "owl": outcome weighted learning, requires `rist` in `data`, imputed survival times
  #         - "rwl": residual weighted learning, requires `rist` in `data`, imputed survival times
  #         - "csf": causal survival forests
  #   - covs: vector of covariates, assumed to be same for all of propensity score model,
  #           probability of being uncensored model, and treatment rule inputs
  #   - prop_scores: some methods require user to pass in propensity scores
  #   - censor: matrix (length of data x number of unique times), some methods require user to pass 
  #             in matrix of censoring survival probabilities
  #   - nfolds: number of folds in cross-validation
  #   - reps: number of repetitions for cross-validation
  #   - extra_args: named list, some methods require extra arguments
  
  N = nrow(data)
  
  # make sure sorted by obs times
  stopifnot(!is.unsorted(data$obs_times))
  
  # default to LOOCV
  if (is.null(nfolds)) (nfolds = N)
  if (is.null(reps)) (reps = 1)
  
  # for each repetition
  treat_rules_stsh = list()
  eta = list()
  
  for (i in 1:reps) {
    
    # CV indices
    if (nfolds==N) (folds = sample(1:N, size=N, replace=FALSE))
    else (folds = sample(1:nfolds, size = N, replace = TRUE))
    
    # treatment rules within each rep
    treat_rules_stsh[[i]] = list()
    eta[[i]] = list()
    
    for (j in 1:nfolds) {
      
      # set up data
      data.train = data[folds != j, ]
      data.test = data[folds == j, ]
      
      data.test.order = order(data.test$obs_times)
      data.test = data.test[data.test.order, ]
      
      # compute treatment rule using training data
      if (method=='cox') {
        
        # type of penalization
        penalties = c('lasso', 'elastic', 'ridge')
        if (any(penalties %in% extra_args)) (penalty=penalties[penalties %in% extra_args])
        else (penalty = 'lasso')
        
        # formatting / prep work
        cox_frmla = as.formula(paste0('~ treat * (', paste0(covs, collapse='+'), ')'))
        
        # get treatment rule
        cox_res = get_cox_dtr(data.train, frmla=cox_frmla, lambda='min', penalty=penalty, tpt=60, test_data=data.test)
        treat_rule = cox_res$treat_rule
        eta[[i]][[j]] = cox_res$eta
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)

      }
      
      if (method=='genetic') {
        
        # requires treatment and censoring weights
        stopifnot(!is.null(prop_scores))
        stopifnot(!is.null(censor))
        
        # smooth or not smooth
        if ('smooth' %in% extra_args) (smooth=TRUE)
        else (smooth=FALSE)
        
        # augmented or not augmented
        if ('augmented' %in% extra_args) (augmented=TRUE)
        else (augmented=FALSE)
        
        # initial values to genetic algorithm
        if ('init_vals' %in% names(extra_args)) (init_vals = extra_args$init_vals)
        else (init_vals=NULL)
        
        # max generations to genetic algorithm
        if ('max.generations' %in% names(extra_args)) (max.generations = extra_args$max.generations)
        else (max.generations=100)
        
        # population size to genetic algorithm
        if ('pop.size' %in% names(extra_args)) (pop.size = extra_args$pop.size)
        else (pop.size=1000)
        
        # formatting / prep work
        X = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                         data = data)
        X = X[,-1]
        
        X.train = X[folds !=j, ]
        X.test = X[folds == j, ]
        
        dat = cbind(treat=data.train$treat, X.train, delta=data.train$delta, obs_times=data.train$obs_times) %>% as.data.frame()
        
        # get treatment rule
        genetic_res = get_genetic_dtr(data=dat, smooth=smooth, augmented=augmented, prop_scores=prop_scores[folds!=j], censor=censor[folds!=j,],
                                      covs=covs, max.generations=max.generations, pop.size=pop.size, test_data=X.test, init_vals=init_vals)
        treat_rule = genetic_res$treat_rule
        eta[[i]][[j]] = genetic_res$eta
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
      }
      
      
      if (method=='owl') {
        
        # requires treatment weights
        stopifnot(!is.null(prop_scores))
        
        # kernel
        if ('kernel' %in% names(extra_args)) (kernel = extra_args$kernel)
        else (kernel='linear')
        
        # cost
        if ('costs' %in% names(extra_args)) (costs = extra_args$costs)
        else (costs=NULL)
        
        # run owl w/ rist
        owl = get_owl_dtr(data=data.train, kernel=kernel, covs=covs, prop_scores=prop_scores[folds!=j],
                          costs=costs, tpt=60, test_data=data.test)
        
        # return treatment rule
        treat_rule = owl$treat_rule
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
        # get coefficients
        eta[[i]][[j]] = as.vector(owl$eta)

      }
      
      
      if (method=='rwl') {
        
        # cost
        if ('costs' %in% names(extra_args)) (costs = extra_args$costs)
        else (costs=NULL)
        
        # run rwl w/ rist
        rwl = get_rwl_dtr(data=data.train, covs=covs, costs=costs, tpt=60, test_data=data.test)
        
        # return treatment rule
        treat_rule = rwl$treat_rule
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
        # get coefficients
        eta[[i]][[j]] = as.vector(rwl$eta)
        
      }
      
      
      if (method=='csf') {
        
        # requires treatment and censoring weights
        stopifnot(!is.null(prop_scores))
        stopifnot(!is.null(censor))
        
        # covariates for training
        X = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                         data = data.train)
        X = X[,-1]
        
        # covariates for test
        X.test = model.matrix(as.formula(paste0('~ ', paste0(covs, collapse='+'))),
                              data = data.test)
        X.test = X.test[,-1]
        
        # run csf
        csf = get_csf_dtr(data.train, X, prop_scores=prop_scores[folds!=j], censor=censor[folds!=j,], test_data=X.test)
        
        # get treatment rule
        treat_rule = csf$treat_rule
        treat_rules_stsh[[i]][[j]] = data.frame(id=data.test$id, treat_rule=treat_rule)
        
        # leave coefficients blank
        eta[[i]][[j]] = NA
        
      }
      
    }
    
  }
  
  return(list(treat_rule = treat_rules_stsh,
              folds = folds,
              eta = eta))
  
}



# Value estimators --------------------------------------------------------


ipwe_km <- function(u, z, delta, prop_scores, censor, treat_rule) {
  
  # K-M value estimator from jiang et al.
  
  # Inputs:
  #   - u: event times
  #   - z: treatment indicator
  #   - delta: event indicator
  #   - prop_scores: propensity scores
  #   - censor: matrix, probability of being uncensored
  #   - treat_rule: treatment recommendations

  # Outputs:
  #   - vector of survival probabilities with length
  #     equal to number of unique event times
  
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(u))
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  # compute quantity inside \prod
  unique_times = 1:ncol(censor)
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in unique_times) {
    index = (u==unique_times[i])
    num[i] = sum(w[index & delta, i])
    den[i] = sum(w[u >= unique_times[i], i])
    if (num[i]==0 & den[i]==0) (den[i]=1)
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0
  
  return(curv)
}



aipwe_km <- function(u, z, delta, prop_scores, censor, failure.surv, lambda.surv, treat_rule) {
  
  # same as ipwe_km but with augmentation
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(u))
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  # compute quantity inside \prod
  Ntimes = ncol(censor)
  unique_times = 1:Ntimes
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in 1:Ntimes) {
    index = (u==unique_times[i])
    num[i] = sum(w[index, i] * delta[index]) + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i] * lambda.surv[index,i])
    den[i] = sum(w[u>unique_times[i] | (index & delta==0), i]) + num[i] + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i])
    if (num[i]==0 & den[i]==0) (den[i]=1)
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0
  
  return(curv)
}



ipwe_km_opt <- function(eta, X, obs_times, z, delta, prop_scores, censor, smooth=FALSE) {
  
  # K-M value estimator from jiang et al. - to be passed into Genetic.IPWE
  
  # same as ipwe_km except takes as a function eta instead of treatment rule
  # eta are coefficients to a linear treatment rule with covariates X
  # outputs estimated survival probability at last time point only
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(obs_times))
  
  len_eta = length(eta)
  
  # smooth treatment rule?
  if (smooth) {
    sd.etaX <- sd(c(cbind(1,X) %*% eta))
    if (!is.finite(sd.etaX)) return(-1000)
    if (sd.etaX > 0) eta <- eta/sd.etaX else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, len_eta-1))
    treat_rule = pnorm(c(cbind(1,X) %*% eta)/((length(z)/4)^(-1/3)))
  }
  
  else {
    treat_rule = as.numeric((cbind(1,X) %*% eta) >= 0)
  }
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  # compute quantity inside \prod
  Ntimes = ncol(censor)
  unique_times = 1:Ntimes
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in 1:Ntimes) {
    index = (obs_times==unique_times[i])
    num[i] = sum(w[index, i] * delta[index])
    den[i] = sum(w[obs_times>unique_times[i] | (index & delta==0), i]) + num[i]
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0
  
  # output only last time point
  return(curv[length(curv)])
}


aipwe_km_opt <- function(eta, X, obs_times, z, delta, prop_scores, censor, failure.surv.0, failure.surv.1, lambda.surv.0, lambda.surv.1, smooth=FALSE) {
  
  # K-M value estimator from jiang et al. - to be passed into Genetic.IPWE
  
  # same as aipwe_km except takes as a function eta instead of treatment rule
  # eta are coefficients to a linear treatment rule with covariates X
  # outputs estimated survival probability at last time point only
  
  # assume all vectors sorted by u
  stopifnot(!is.unsorted(obs_times))
  
  len_eta = length(eta)
  
  # smooth treatment rule?
  if (smooth) {
    sd.etaX <- sd(c(cbind(1,X) %*% eta))
    if (!is.finite(sd.etaX)) return(-1000)
    if (sd.etaX > 0) eta <- eta/sd.etaX else eta <- c(ifelse(eta[1] >= 0, 1, -1), rep(0, len_eta-1))
    treat_rule = pnorm(c(cbind(1,X) %*% eta)/((length(z)/4)^(-1/3)))
  }
  
  else {
    treat_rule = as.numeric((cbind(1,X) %*% eta) >= 0)
  }
  
  # inputed propensity scores should be P(Z=1|X)
  # change to P(Z=z|X)
  prop_scores2 = prop_scores*(z) + (1-prop_scores)*(1-z)
  follow_treat = z*treat_rule + (1-z)*(1-treat_rule)
  
  failure.surv = failure.surv.1*treat_rule + failure.surv.0*(1-treat_rule)
  lambda.surv = lambda.surv.1*treat_rule + lambda.surv.0*(1-treat_rule)
  
  # compute quantity inside \prod
  Ntimes = ncol(censor)
  unique_times = 1:Ntimes
  w = follow_treat / (prop_scores2 * censor)  # weights
  num = c()
  den = c()
  for (i in 1:Ntimes) {
    index = (obs_times==unique_times[i])
    num[i] = sum(w[index, i] * delta[index]) + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i] * lambda.surv[index,i])
    den[i] = sum(w[obs_times>unique_times[i] | (index & delta==0), i]) + num[i] + sum((1-w[index, i]) * failure.surv[index,i] * censor[index,i])
    if (num[i]==0 & den[i]==0) (den[i]=1)
  }
  
  # compute K-M curve
  curv = cumprod(1-num/den)
  curv[is.na(curv)] = 0
  
  # output only last time point
  return(curv[length(curv)])
}



# Eval functions ----------------------------------------------------------


extract_risk <- function(treat_list, data, prop_scores, censor) {
  
  # extract risk estimates (do this for each imputation)
  
  # Inputs:
  # treat_list: treat_rule output from cv_estimator function
  # data: data, sorted by obs time
  # prop_scores: propensity scores
  # censor: censoring weights
  # Outputs:
  # single risk curve
  
  
  # index to arrange data in original order (i.e. by id) -> sorted obs_times order
  u.order = rank(data$id, ties='first')
  
  # combine treatment recommendations for each M and arrange
  # i.e. get list of df where each df is treat rule for all data
  treat_list_flat = lapply(treat_list, function(x) bind_rows(lapply(x, function(y) as.data.frame(y))))
  treat_list_flat = lapply(treat_list_flat, function(x) x %>% arrange(id) %>% slice(u.order))
  
  # for each M, compute risk estimate
  risk_M = lapply(treat_list_flat, function(x) 1-ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, censor, x$treat_rule))
  
  # format
  risk_M_df = risk_M %>% as.data.frame() %>% t()
  
  # return average across the M
  return(colMeans(risk_M_df))
}


extract_treat <- function(treat_list, data) {
  
  # extract risk estimates (do this for each imputation)
  
  # Inputs:
  # treat_list: treat_rule output from cv_estimator function
  # data: data, sorted by obs time
  # Outputs:
  # vector of treatment rules, in same order as `data`
  
  
  # only works for M=1
  stopifnot('M must be 1' = length(treat_list)==1)
  
  # index to arrange data in original order (i.e. by id) -> sorted obs_times order
  u.order = rank(data$id, ties='first')
  
  # combine treatment recommendations for each M and arrange
  # i.e. get list of df where each df is treat rule for all data
  treat_list_flat = lapply(treat_list, function(x) bind_rows(lapply(x, function(y) as.data.frame(y))))
  treat_list_flat = lapply(treat_list_flat, function(x) x %>% arrange(id) %>% slice(u.order))
  
  
  # treatment rule as named vector
  treat_rule = treat_list_flat[[1]]$treat_rule
  names(treat_rule) = treat_list_flat[[1]]$id
  return(treat_rule)
}


# Other -------------------------------------------------------------------


do_rist <- function(data, nmin=6, M=50, L=2, tao=60) {
  
  # data: data must be ordered X, delta, times
  # nmin: minimum number of observed data in each node
  # M: number of trees in each fold
  # L: number of folds
  # tao: length of study
  
  # below file from: https://sites.google.com/site/teazrq/software
  source('RISTfunctions.r')
  
  # set other hyperparameters
  P = ncol(rist_data)-2   # number of dimension for X
  K = round(sqrt(P))
  
  # split dataset into not censored
  data_c = data[data$delta==0, ]
  
  # RIST
  R_Muti_ERT_build = Muti_ERT_fit(as.matrix(data[, 1:P]), M, K, L, nmin, SupLogRank=1, tao=tao, impute="random")
  R_Muti_ERT_predict = Muti_ERT_Predict(as.matrix(data_c), R_Muti_ERT_build$Forest_seq[[L]], 
                                        R_Muti_ERT_build$SurvMat_seq[[L]], R_Muti_ERT_build$time_intrest)
  
  # predictions
  times_pred_c = R_Muti_ERT_predict$TMT_predict
  times_pred = data$obs_times
  times_pred[data$delta==0] = times_pred_c
  
  # return observed times with censored times imputed
  return(times_pred)
  
}


Genetic.IPWE <- function(fn, X, u, z, delta, prop_scores, censor, smooth, nvars, standardize,
                         max.generations, pop.size, aug_terms=NULL, init_vals=NULL) {  
  
  # set initial values to 0 if not specified
  if (is.null(init_vals)) (init_vals = rep(0,nvars))
  
  # make sure number of initial values equals number of variables
  stopifnot(length(init_vals) == nvars)
  
  # standardize X?
  if (standardize) {
    X = scale(X)
  }
  
  # augmented version
  if (!is.null(aug_terms)) {
    
    failure.surv.0 = aug_terms$failure.surv.0
    failure.surv.1 = aug_terms$failure.surv.1
    lambda.surv.0 = aug_terms$lambda.surv.0
    lambda.surv.1 = aug_terms$lambda.surv.1
    
    temp <- genoud(fn=fn, X=X, obs_times=u, z=z, delta=delta, prop_scores=prop_scores, censor=censor, 
                   failure.surv.0=failure.surv.0, failure.surv.1=failure.surv.1, 
                   lambda.surv.0=lambda.surv.0, lambda.surv.1=lambda.surv.1, smooth=smooth,
                   nvars=nvars, 
                   Domains=cbind(rep(-5,nvars),rep(5,nvars)),
                   starting.values=init_vals,
                   max=TRUE,
                   pop.size = pop.size,
                   
                   print.level=1,
                   BFGS=FALSE, 
                   optim.method="Nelder-Mead",
                   P9=0, 
                   
                   max.generations=max.generations,
                   hard.generation.limit=TRUE)
  }
  
  # non-augmented version
  else {
    temp <- genoud(fn=fn, X=X, obs_times=u, z=z, delta=delta, prop_scores=prop_scores, censor=censor, smooth=smooth,
                   nvars=nvars, 
                   Domains=cbind(rep(-5,nvars),rep(5,nvars)),
                   starting.values=init_vals,
                   max=TRUE,
                   pop.size = pop.size,
                   
                   print.level=1,
                   BFGS=FALSE, 
                   optim.method="Nelder-Mead",
                   P9=0, 
                   
                   max.generations=max.generations,
                   hard.generation.limit=TRUE)
  }
  
  # extract estimates
  eta.est <- temp$par
  if (sum(eta.est^2)==0) eta.est <- c(1,rep(0,nvars-1))
  valhat.etahat <- temp$value
  
  return(c(eta.est, valhat.etahat))
}


# function to center and scale a vector by 2 sd
scale_2sd <- function(vec) {
  vec = vec - mean(vec)
  sdx = sd(vec)
  vec = vec / (2*sdx)
  return(vec)
}
