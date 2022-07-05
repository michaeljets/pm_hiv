# Creates a simulated dataset 

source('utils.R')

library(survival)
library(tidyverse)


# set seed for reproducibility
set.seed(1162819)

# set size of data
n = 1000

# simulate covariates
data = data.frame(X1 = runif(n, -0.5, 1),
                  X2 = rbinom(n, size=1, prob=0.5),
                  X3 = rbinom(n, size=1, prob=0.5),
                  X4 = rbinom(n, size=1, prob=0.5),
                  X5 = rbinom(n, size=1, prob=0.5))

# create treatment variable (treat), depending on covariates
data$treat = rbinom(n, size=1, prob = plogis(2*data$X1))

# parameters for censor distribution, Gompertz
censor.linear.part = 0.2 + 0.8*data$X2
censor.lambda = 0.008  # scale
censor.gamma = 0.02  # shape

# plot true censor curves
true_cens = GenSurvCurv(censor.linear.part, censor.lambda, censor.gamma)
plot(1:60, apply(true_cens, 2, mean), type='l', xlab='time', ylab='survival', main='Censoring survival curve', ylim=c(0,1))

# generate censoring times from Gompertz distribution
# with administrative censoring at time 60
censor_data = GenSimCens(censor.linear.part, censor.lambda, censor.gamma, end=60)
data$censor_times = censor_data$obs_censor


# parameters for survival event distribution, Gompertz
event.lambda = 0.001  # scale
event.gamma = 0.02  # shape
coarsen = TRUE

# linear part under different treatment rules
linear.part.0 = 1 + 1.5*data$X1
linear.part.1 = 1 + 1.5*data$X1 + 1 - 2*data$X2
linear.part.d = pmin(linear.part.0, linear.part.1)

# plot true survival curves
true_surv_0 = GenSurvCurv(linear.part.0, event.lambda, event.gamma)
true_surv_1 = GenSurvCurv(linear.part.1, event.lambda, event.gamma)
true_surv_d = GenSurvCurv(linear.part.d, event.lambda, event.gamma)

plot(1:60, 1-apply(true_surv_0, 2, mean), col='red', type='l', xlab='time', ylab='risk', main='True risk curves', ylim=c(0,0.5))
lines(1:60, 1-apply(true_surv_1, 2, mean), col='blue', type='l')
lines(1:60, 1-apply(true_surv_d, 2, mean), col='purple', type='l')
legend('bottomright', legend = c('treat 0 only', 'treat 1 only', 'dtr'), col = c('red', 'blue', 'purple'), lty = 1)

# risk difference at 60 months
apply(true_surv_d, 2, mean)[60] - max(apply(true_surv_0, 2, mean)[60], apply(true_surv_1, 2, mean)[60])

# generate event data from Gompertz distribution
event.linear.part = 1 + 1.5*data$X1 + data$treat*(1 - 2*data$X2)
times = GenSimTimes(event.linear.part, event.lambda, event.gamma,
                    data$censor_times, coarsen=coarsen, end=end)
data$true_times = times$true_times
data$delta = times$delta
data$obs_times = times$obs_times
data$ltfu = ifelse(data$delta==1, 0, censor_data$ltfu)

# add id column
data$id = 1:nrow(data)

# sort by observed time
data = data %>% arrange(obs_times, id)

# for RMST, last time point is considered not censored
data = data %>%
  mutate(r_delta = ifelse(obs_times==60, 1, delta))

# save data
saveRDS(data, 'data.rds')
