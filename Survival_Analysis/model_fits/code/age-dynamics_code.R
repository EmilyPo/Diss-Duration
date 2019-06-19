# code for age explorations 

library(tidyverse)
library(here)
library(flexsurv)

#### Load data ####
dat <- readRDS(here("model_fits", "survdat.rds"))

########################################################################
# try using flexsurv package for regression
# similar to 'survival' but allows for left truncated data in parametric regression
# NOTE: can't do mixture models 
#######################################################################

#### baseline exp & weib & gamma w/o covariates ####

expFlex <- flexsurvreg(Surv(t_o, t_c, censored) ~ 1, 
                       data = dat, dist = "exp", weights = e.weight)

weibFlex <- flexsurvreg(Surv(t_o, t_c, censored) ~ 1, 
                        data = dat, dist = "weibull", weights = e.weight)

gammaFlex <- flexsurvreg(Surv(t_o, t_c, censored) ~ 1, 
                        data = dat, dist = "gamma", weights = e.weight)

saveRDS(expFlex, "model_fits/nocovs/exp.rds")
saveRDS(weibFlex, "model_fits/nocovs/weib.rds")
saveRDS(gammaFlex, "model_fits/nocovs/gamma.rds")

######## continuous covariates ##########
# ignoring gamma model from here on out 

## current ego age 
ea <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age, 
                      data = dat, dist = "exp", weights = e.weight) 

eaWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age, 
                          data = dat, dist = "weibull", weights = e.weight) 

saveRDS(ea, "model_fits/age/expEgoAgeCurrent.rds")
saveRDS(eaWeib, "model_fits/age/weibEgoAgeCurrent.rds")

## initial ego age 
ea.init <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age.initial, 
                           data = dat, dist = "exp", weights = e.weight) 

ea.initWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age.initial, 
                               data = dat, dist = "weibull", weights = e.weight) 

saveRDS(ea.init, "model_fits/age/expEgoAgeInit.rds")
saveRDS(ea.initWeib, "model_fits/age/weibEgoAgeInit.rds")


## diff sqrt age
ageDiff <- flexsurvreg(Surv(t_o, t_c, censored) ~ diff.sqrt.age, 
                       data = dat, dist = "exp", weights = e.weight)

ageDiffWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ diff.sqrt.age, 
                           data = dat, dist = "weibull", weights = e.weight)

saveRDS(ageDiff, "model_fits/age/expAgeDiff.rds")
saveRDS(ageDiffWeib, "model_fits/age/weibAgeDiff.rds")

# current age + diff sqrt age
AgeDiffCur <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age + diff.sqrt.age, 
                           data = dat, dist = "exp", weights = e.weight)

AgeDiffCurWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age + diff.sqrt.age, 
                               data = dat, dist = "weibull", weights = e.weight)

saveRDS(AgeDiffCur, "model_fits/age/expAgeDiff_CurAge.rds")
saveRDS(AgeDiffCurWeib, "model_fits/age/weibAgeDiff_CurAge.rds")


##### categorical covariates ######
# age cat current
eaCatFlex <- flexsurvreg(Surv(t_o, t_c, censored) ~ as.factor(e.agecat), 
                         data = dat, dist = "exp", weights = e.weight) 

eaCatFlexWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ as.factor(e.agecat), 
                             data = dat, dist = "weibull", weights = e.weight) 

saveRDS(eaCatFlex, "model_fits/age/expAgeCatCurrent.rds")
saveRDS(eaCatFlexWeib, "model_fits/age/weibAgeCatCurrent.rds")

eaCatAgeDiff <- flexsurvreg(Surv(t_o, t_c, censored) ~ as.factor(e.agecat) + diff.sqrt.age, 
                            data = dat, dist = "exp", weights = e.weight)

eaCatAgeDiffWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ as.factor(e.agecat) + diff.sqrt.age, 
                                data = dat, dist = "weibull", weights = e.weight) 

saveRDS(eaCatAgeDiff, "model_fits/age/expAgeDiffAgeCatCurrent.rds")
saveRDS(eaCatAgeDiffWeib, "model_fits/age/weibAgeDiffAgeCatCurrent.rds")

