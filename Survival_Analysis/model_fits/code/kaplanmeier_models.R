########################################################
## code for kaplan-meier models 
########################################################

##### libs ######
library(here)
library(survival)

#### Load data ####
dat <- readRDS("~/NSFG_DATA/Objects/altersegos_survdat.rds")

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

# no log-rank test for left-truncated data
survdiff(Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat,
         data = dat)
# lets try to fit a coz ph model to get a statistic


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

e.agecat.coxph <- coxph(Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat,
                        data = dat,
                        weights = e.weight,
                        robust = T)
saveRDS(e.agecat.coxph, "model_fits/km/agecat-coxph.rds")

## initial ego age category
km_i_agecat_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat.initial,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_i_agecat_weighted, "model_fits/km/km_i_agecat_weighted.rds")

## by reltype
km_reltype_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ reltype,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_reltype_weighted, "model_fits/km/km_reltype_weighted.rds")

reltype.coxph <- coxph(Surv(time = t_o, time2 = t_c, event = censored) ~ reltype,
                        data = dat,
                        weights = e.weight,
                        robust = T)
saveRDS(reltype.coxph, "model_fits/km/reltype-coxph.rds")


## by ego # parts in last year
km_eparts_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.partsyr3,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_eparts_weighted, "model_fits/km/km_eparts_weighted.rds")


partsyr.coxph <- coxph(Surv(time = t_o, time2 = t_c, event = censored) ~ e.partsyr3,
                       data = dat,
                       weights = e.weight,
                       robust = T)
saveRDS(partsyr.coxph, "model_fits/km/partsyr-coxph.rds")

## by ego race
km_erace_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.race,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_erace_weighted, "model_fits/km/km_erace_weighted.rds")

erace.coxph <- coxph(Surv(time = t_o, time2 = t_c, event = censored) ~ e.race,
                       data = dat,
                       weights = e.weight,
                       robust = T)
saveRDS(erace.coxph, "model_fits/km/erace-coxph.rds")


## by network (mar/cohab and other)
km_network_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ network1,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_network_weighted, "model_fits/km/km_network_weighted.rds")

network.coxph <- coxph(Surv(time = t_o, time2 = t_c, event = censored) ~ network1,
                       data = dat,
                       weights = e.weight,
                       robust = T)
saveRDS(network.coxph, "model_fits/km/network-coxph.rds")
