########################################################
## code for kaplan-meier models 
########################################################

##### libs ######
library(here)
library(survival)

#### Load data ####
dat <- readRDS(here("model_fits", "survdat.rds"))

#### Modifed Kaplan-Meier Models ####

## no covariates
km <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ 1,
  data = dat,
  error = "greenwood"
)

saveRDS(km, "model_fits//km/km-nsfg.rds")

## current ego age category
km_agecat_current <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat,
  data = dat,
  error = "greenwood"
)

saveRDS(km_agecat_current, "model_fits/km/km-agecat-current.rds")

## initial ego age category
km_agecat_initial <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat.initial,
  data = dat,
  error = "greenwood"
)

saveRDS(km_agecat_initial, "model_fits/km/km-agecat-initial.rds")

#### Weighted K-Ms ############################

## no covariates
km_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ 1,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_weighted, "model_fits/km/km-weighted.rds")

## current ego age category
km_agecat_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_agecat_weighted, "model_fits/km/km_ac_weighted.rds")

## initial ego age category
km_i_agecat_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat.initial,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_i_agecat_weighted, "model_fits/km/km_i_agecat_weighted.rds")


