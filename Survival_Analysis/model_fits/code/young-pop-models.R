
#############################################################
# Restrict Data to egos & alters aged 15-29
#############################################################

library(tidyverse)
library(here)
library(flexsurv)
library(survival)
library(maxLik)

#### Load data ####
dat <- readRDS(here("model_fits", "survdat.rds"))
young.dat <- dat %>% filter(e.age < 30) %>% filter(age < 30)
x <- young.dat

#### kaplan-meier ######
km.young.unweighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ 1,
  data = young.dat,
  error = "greenwood"
)

km.young.weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ 1,
  data = young.dat, weights = e.weight,
  error = "greenwood"
)

## current ego age category
km.young_agecat_unweighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat,
  data = young.dat,
  error = "greenwood"
)

km.young_agecat_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat,
  weights = e.weight,
  data = young.dat,
  error = "greenwood"
)

saveRDS(km.young.unweighted, "model_fits/1529/km-young-unweighted.rds")
saveRDS(km.young.weighted, "model_fits/1529/km-young-weighted.rds")
saveRDS(km.young_agecat_unweighted, "model_fits/1529/km-young_agecat_unweighted.rds")
saveRDS(km.young_agecat_weighted, "model_fits/1529/km-young_agecat_weighted.rds")

######## flexsurv models ######################
# no covs for reference
exp.young <- flexsurvreg(Surv(t_o, t_c, censored) ~ 1, 
                         data = young.dat, dist = "exp", weights = e.weight)

weib.young <- flexsurvreg(Surv(t_o, t_c, censored) ~ 1, 
                          data = young.dat, dist = "weibull", weights = e.weight)

saveRDS(exp.young, "model_fits/1529/exp.young.rds")
saveRDS(weib.young, "model_fits/1529/weib.young.rds")

# ego current age 
expYoungCurrAge <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age, 
                               data = young.dat, dist = "exp", weights = e.weight)

weibYoungCurrAge <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.age, 
                                data = young.dat, dist = "weibull", weights = e.weight)

saveRDS(expYoungCurrAge, "model_fits/1529/expYoungCurrAge.rds")
saveRDS(weibYoungCurrAge, "model_fits/1529/weibYoungCurrAge.rds")

# sqrt age diff 
expYoungAgeDiff <- flexsurvreg(Surv(t_o, t_c, censored) ~ diff.sqrt.age, 
                               data = young.dat, dist = "exp", weights = e.weight)

weibYoungAgeDiff <- flexsurvreg(Surv(t_o, t_c, censored) ~ diff.sqrt.age, 
                                data = young.dat, dist = "weibull", weights = e.weight) 

saveRDS(expYoungAgeDiff, "model_fits/1529/expYoungAgeDiff.rds")
saveRDS(weibYoungAgeDiff, "model_fits/1529/weibYoungAgeDiff.rds")

# age cat current 
expYoungAgeCat <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat, 
                              data = young.dat, dist = "exp", weights = e.weight)

weibYoungAgeCat <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat, 
                               data = young.dat, dist = "weibull", weights = e.weight) 

saveRDS(expYoungAgeCat, "model_fits/1529/expYoungAgeCat.rds")
saveRDS(weibYoungAgeCat, "model_fits/1529/weibYoungAgeCat.rds")


######## latent mixture models ################

# here: latent model w/ weights for each category
logit = function(p){
  log(p/(1-p))
}

invLogit = function(phi){
  1/(1+exp(-phi))
}

#### NO COVS
mixML <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight1 = param[3]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  
  ifelse(censored == 1,
         yes = log(
           (dexp(closed, exp(logLambda1)) * invLogit(logitWeight1) +
              dexp(closed, exp(logLambda2)) * (1 - invLogit(logitWeight1))) / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
         ),
         no = log(
           (pexp(closed, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
              pexp(closed, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1))) / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
         )
  )
  
}

mixExpStartParam = c(
  logLambda1 = log(0.2),
  logLambda2 = log(0.01),
  logitWeight1 = logit(0.5)
)

mixExpFityoung = maxLik(
  logLik = mixML,
  start = mixExpStartParam,
  method = "BHHH"
)

summary(mixExpFityoung)
exp(coef(mixExpFityoung)[1:2])
invLogit(coef(mixExpFityoung)[3])
AIC(mixExpFityoung)

saveRDS(mixExpFityoung, "model_fits/1529/mixExpYoung.rds")


#### WITH AGE 
# for ind age cat models, can pull from "model_fits/mixed..."


#### mixed exp with weights by age cat #################

eaMixMLweights.young <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight15 = param[3]
  logitWeight20 = param[4]
  logitWeight25 = param[5]
  
  agecat15 <- x$e.15
  agecat20 <- x$e.20
  agecat25 <- x$e.25
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  weights <- logitWeight15*agecat15 + logitWeight20*agecat20 + logitWeight25*agecat25
  
  ifelse(censored == 1,
         yes = log(
           (dexp(closed, exp(logLambda1)) * invLogit(weights) +
              dexp(closed, exp(logLambda2)) * (1 - invLogit(weights)))  / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(weights) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(weights)))
         ),
         no = log(
           (pexp(closed, exp(logLambda1), lower.tail = FALSE) * invLogit(weights) +
              pexp(closed, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(weights))) / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(weights) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(weights)))
         )
  )
  
}

mixExpStartParam = c(
  logLambda1 = log(0.2),
  logLambda2 = log(0.01),
  logitWeight15 = logit(0.5),
  logitWeight20 = logit(0.5),
  logitWeight25 = logit(0.5)
)


weightsMixExp.young <-  maxLik(logLik = eaMixMLweights.young, start = mixExpStartParam, method = "BHHH")

summary(weightsMixExp.young)
saveRDS(weightsMixExp.young, "model_fits/1529/weightsMixExp-young.rds")
