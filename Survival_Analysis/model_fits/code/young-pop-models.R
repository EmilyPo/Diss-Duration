
#############################################################
# Restrict Data to egos & alters aged 15-29
#############################################################

library(tidyverse)
library(here)
library(flexsurv)
library(survival)
library(maxLik)

#### Load data ####
dat <- readRDS("~/NSFG_DATA/Objects/altersegos_survdat.rds")
dat <- dat %>%
  filter(e.age < 30) %>% 
  filter(age < 30) %>%
  mutate(e.race = as.factor(e.race),
         reltype = as.factor(reltype),
         e.partsyr3 = as.factor(e.partsyr3),
         e.osnpyr3 = as.factor(e.osnpyr3),
         e.maxospyr3 = as.factor(e.maxospyr3), 
         network1 = as.factor(network1), 
         e.agecat = as.factor(e.agecat),
         e.agecat.initial = as.factor(e.agecat.initial),
         e.deg.main = as.factor(e.deg.main))  

x <- dat

#### K-Ms ############################

## no covariates
km_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ 1,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_weighted, "model_fits/1529/km-weighted.rds")

km_unweighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ 1,
  data = dat,
  error = "greenwood"
)

saveRDS(km_unweighted, "model_fits/1529/km-unweighted.rds")

## current ego age category
km_agecat_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_agecat_weighted, "model_fits/1529/km_ac_weighted.rds")

## initial ego age category
km_i_agecat_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.agecat.initial,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_i_agecat_weighted, "model_fits/1529/km_i_agecat_weighted.rds")

## by reltype
km_reltype_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ reltype,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_reltype_weighted, "model_fits/1529/km_reltype_weighted.rds")


## by ego # parts in last year
km_eparts_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.partsyr3,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_eparts_weighted, "model_fits/1529/km_eparts_weighted.rds")

## by ego race
km_erace_weighted <- survfit(
  formula = Surv(time = t_o, time2 = t_c, event = censored) ~ e.race,
  weights = e.weight,
  data = dat,
  error = "greenwood"
)

saveRDS(km_erace_weighted, "model_fits/1529/km_erace_weighted.rds")


######## Flexsurv models ######################
# no covs for reference
expFlex <- flexsurvreg(Surv(t_o, t_c, censored) ~ 1, 
                       data = dat, dist = "exp", weights = e.weight)

weibFlex <- flexsurvreg(Surv(t_o, t_c, censored) ~ 1, 
                        data = dat, dist = "weibull", weights = e.weight)

saveRDS(expFlex, "model_fits/1529/expFlex.rds")
saveRDS(weibFlex, "model_fits/1529/weibFlex.rds")


#----------------------- AGE CAT ---------------------------------------------------
  
e.agecat <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat, 
                              data = dat, dist = "exp", weights = e.weight)

w.agecat <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat, 
                               data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.agecat, "model_fits/1529/eagecat.rds")
saveRDS(w.agecat, "model_fits/1529/wagecat.rds")

# ----------------------- RACE -------------------------------------------------------

e.race <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race, 
                      data = dat, dist = "exp", weights = e.weight) 

w.race <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race, 
                      data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.race, "model_fits/1529/erace.rds")
saveRDS(w.race, "model_fits/1529/wrace.rds")

# ----------------------- RELTYPE ----------------------------------------------------

e.reltype <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype, 
                         data = dat, dist = "exp", weights = e.weight) 

saveRDS(e.reltype, "model_fits/1529/ereltype.rds")

w.reltype <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype, 
                         data = dat, dist = "weibull", weights = e.weight) 

saveRDS(w.reltype, "model_fits/1529/wreltype.rds")

# ----------------------- NPARTS LAST YEAR -------------------------------------------

e.partsyr <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.partsyr3, 
                         data = dat, dist = "exp", weights = e.weight) 

saveRDS(e.partsyr, "model_fits/1529/epartsyr.rds")

w.partsyr <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.partsyr3, 
                         data = dat, dist = "weibull", weights = e.weight) 

saveRDS(w.partsyr, "model_fits/1529/wpartsyr.rds")



# ----------------------- MOMENTARY DEGREE -------------------------------------------

e.deg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main, 
                     data = dat, dist = "exp", weights = e.weight) 

w.deg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main, 
                     data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.deg, "model_fits/1529/edeg.rds")
saveRDS(w.deg, "model_fits/1529/wdeg.rds")


# ----------------------- network type -------------------------------------------

e.network <- flexsurvreg(Surv(t_o, t_c, censored) ~ network1, 
                         data = dat, dist = "exp", weights = e.weight) 


w.network <- flexsurvreg(Surv(t_o, t_c, censored) ~ network1, 
                         data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.network, "model_fits/1529/enetwork.rds")
saveRDS(w.network, "model_fits/1529/wnetwork.rds")

# ----------------------- MULTIVARIATE -------------------------------------------

e.relparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.partsyr3, 
                          data = dat, dist = "exp", weights = e.weight)

w.relparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.partsyr3, 
                          data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.relparts, "model_fits/1529/erelparts.rds")
saveRDS(w.relparts, "model_fits/1529/wrelparts.rds")

#---

e.ageparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.partsyr3, 
                          data = dat, dist = "exp", weights = e.weight) 

w.ageparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.partsyr3, 
                          data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.ageparts, "model_fits/1529/eageparts.rds")
saveRDS(w.ageparts, "model_fits/1529/wageparts.rds")

#---

e.raceparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.partsyr3, 
                           data = dat, dist = "exp", weights = e.weight)

w.raceparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.partsyr3, 
                           data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.raceparts, "model_fits/1529/eraceparts.rds")
saveRDS(w.raceparts, "model_fits/1529/wraceparts.rds")

#---

e.racerels <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.race, 
                          data = dat, dist = "exp", weights = e.weight)

w.racerels <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.race, 
                          data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.racerels, "model_fits/1529/eracerels.rds")
saveRDS(w.racerels, "model_fits/1529/wracerels.rds")

#---

e.raceage <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.race, 
                         data = dat, dist = "exp", weights = e.weight)

w.raceage <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.race, 
                         data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.raceage, "model_fits/1529/eraceage.rds")
saveRDS(w.raceage, "model_fits/1529/wraceage.rds")

#---

e.agerels <-  flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + reltype, 
                          data = dat, dist = "exp", weights = e.weight)

w.agerels <-  flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + reltype, 
                          data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.agerels, "model_fits/1529/eagerels.rds")
saveRDS(w.agerels, "model_fits/1529/wagerels.rds")

#---

e.racedeg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.deg.main, 
                         data = dat, dist = "exp", weights = e.weight)

w.racedeg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.deg.main, 
                         data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.racedeg, "model_fits/1529/eracedeg.rds")
saveRDS(w.racedeg, "model_fits/1529/wracedeg.rds")

#---

e.agedeg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main + e.agecat, 
                        data = dat, dist = "exp", weights = e.weight)

# here - model unfit. "system is computationally singular"
w.agedeg <-  flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main + e.agecat, 
                         data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.agedeg, "model_fits/1529/eagedeg.rds")
saveRDS(w.agedeg, "model_fits/1529/wagedeg.rds")



################### other age models ######################
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
