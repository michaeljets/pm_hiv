# computes final results

source('utils.R')

library(survival)
library(readr)
library(ggrepel)
library(tidyverse)


# set ggplot theme
theme_set(theme_bw())

# read in data
data = readRDS('data.rds')

# read in CV treatment rules
cv_treat = readRDS('results/cv_treatrules.rds')

# read in standard errors
boot_sd = readRDS('results/se_estimates.rds')

# covariates
covs = paste0('X', 1:5)

# estimate propensity scores
prop_frmla = as.formula(paste0('treat ~', paste0(covs, collapse='+')))
prop_model = glm(prop_frmla, data = data, family = binomial)
prop_scores = predict(prop_model, type = 'response')

# censoring weights
frmla = paste0('Surv(obs_times, ltfu) ~ treat +', paste0(covs, collapse='+'))
censor_frmla = as.formula(frmla)
cens_model = coxph(censor_frmla, data = data)
cens_weights = survfit(cens_model, newdata = data)$surv %>% t()

# non-dtr curves
curv_0 = 1-ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights, 0)
curv_1 = 1-ipwe_km(data$obs_times, data$treat, data$delta, prop_scores, cens_weights, 1)

# dtr curve, using genetic algorithm with augmentation and smoothing
curv_d = extract_risk(cv_treat$genetic_as$treat_rule, data, prop_scores, cens_weights)

# true curve

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

curv_0_true = 1 - colMeans(true_surv_0)
curv_1_true = 1 - colMeans(true_surv_1)
curv_d_true = 1 - colMeans(true_surv_d)



# Plot --------------------------------------------------------------------


# wrangle data
plot.data.point = data.frame(time = 1:60,
                             treat_0 = curv_0,
                             treat_1 = curv_1,
                             dtr = curv_d,
                             true_0 = curv_0_true,
                             true_1 = curv_1_true,
                             true_d = curv_d_true) %>%
  pivot_longer(cols = c('treat_0', 'treat_1', 'dtr', 'true_0', 'true_1', 'true_d'), names_to = 'treat_rule')

plot.data.lower = data.frame(time = 1:60,
                             treat_0 = curv_0 - 1.96*boot_sd$boot_sd_0,
                             treat_1 = curv_1 - 1.96*boot_sd$boot_sd_1,
                             dtr = curv_d - 1.96*boot_sd$boot_sd_dtr) %>%
  pivot_longer(cols = c('treat_0', 'treat_1', 'dtr'), names_to = 'treat_rule', values_to = 'lower')

plot.data.upper = data.frame(time = 1:60,
                             treat_0 = curv_0 + 1.96*boot_sd$boot_sd_0,
                             treat_1 = curv_1 + 1.96*boot_sd$boot_sd_1,
                             dtr = curv_d + 1.96*boot_sd$boot_sd_dtr) %>%
  pivot_longer(cols = c('treat_0', 'treat_1', 'dtr'), names_to = 'treat_rule', values_to = 'upper')

plot.data = left_join(plot.data.point, plot.data.lower) %>% left_join(plot.data.upper) %>%
  mutate(label = ifelse(time!=max(time), NA, treat_rule))

# add time 0 to plots
dat0 = tibble(time = 0,
              treat_rule = c('treat_0', 'treat_1', 'dtr', 'true_0', 'true_1', 'true_d'),
              value = 0,
              lower = NA,
              upper = NA,
              label = NA)

plot.data = bind_rows(dat0, plot.data)

# put on percent scale
plot.data = plot.data %>%
  mutate(value = value*100,
         lower = lower*100,
         upper = upper*100)

# add color indicator and treatment group variable
plot.data = plot.data %>%
  mutate(true = ifelse(grepl('true', treat_rule), 'true', 'estimated'),
         treat_group = gsub('.*_', '', treat_rule),
         treat_group = ifelse(treat_group=='d', 'dtr', treat_group))

# plots

# without confidence intervals
ggplot(data=plot.data %>% filter(true=='estimated'), 
       aes(x=time, y=value, group=treat_rule, label=label)) +
  geom_step(aes(linetype=treat_rule)) +
  geom_label_repel(nudge_x=2, force=1, min.segment.length = Inf, na.rm=TRUE) +
  ylim(0, 55) +
  labs(x = 'Time (years)',
       y = 'Risk (%)') +
  scale_x_continuous(breaks = seq(0, 60, by=12),
                     labels = 0:5,
                     limits = c(0, 65)) + 
  scale_color_manual(values = c('black', 'red')) +
  theme(legend.position = 'none',
        legend.title = element_blank())

ggsave('results/point_estimates.png', width=5, height=5, units='in')

# confidence intervals
# warning "In max(ids, na.rm = TRUE) : no non-missing arguments to max; returning -Inf" ok
ggplot(data=plot.data %>% filter(true=='estimated') %>% mutate(lower=ifelse(lower<0, 0, lower)), 
       aes(x=time, y=value, label=label, group=treat_rule)) +
  geom_step() +
  pammtools::geom_stepribbon(aes(ymin=lower, ymax=upper), alpha=0.5, fill='grey') +
  facet_wrap(vars(treat_group)) +
  ylim(0, 55) +
  labs(x = 'Time (years)',
       y = 'Risk (%)') +
  scale_x_continuous(breaks = seq(0, 60, by=12),
                     labels = 0:5,
                     limits = c(0, 65)) + 
  theme(legend.position = c(0.9, 0.8),
        legend.title = element_blank())

ggsave('results/ci_estimates.png', width=5, height=5, units='in')



# below plots are same as above but include true risk curves

# # without confidence intervals
# ggplot(data=plot.data, aes(x=time, y=value, group=treat_rule, label=label, color=true)) +
#   geom_step(aes(linetype=treat_rule)) +
#   geom_label_repel(nudge_x=2, force=1, min.segment.length = Inf, na.rm=TRUE) +
#   ylim(0, 55) +
#   labs(x = 'Time (years)',
#        y = 'Risk (%)') +
#   scale_x_continuous(breaks = seq(0, 60, by=12),
#                      labels = 0:5,
#                      limits = c(0, 65)) +
#   scale_color_manual(values = c('black', 'red')) +
#   theme(legend.position = 'none',
#         legend.title = element_blank())
# 
# 
# # confidence intervals
# # warning "In max(ids, na.rm = TRUE) : no non-missing arguments to max; returning -Inf" ok
# ggplot(data=plot.data %>% mutate(lower=ifelse(lower<0, 0, lower)),
#        aes(x=time, y=value, label=label, group=treat_rule, color=true)) +
#   geom_step() +
#   pammtools::geom_stepribbon(aes(ymin=lower, ymax=upper), alpha=0.5, fill='grey') +
#   facet_wrap(vars(treat_group)) +
#   ylim(0, 55) +
#   labs(x = 'Time (years)',
#        y = 'Risk (%)') +
#   scale_x_continuous(breaks = seq(0, 60, by=12),
#                      labels = 0:5,
#                      limits = c(0, 65)) +
#   theme(legend.position = c(0.9, 0.8),
#         legend.title = element_blank())
