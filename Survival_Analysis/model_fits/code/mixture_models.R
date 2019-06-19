library(maxLik)
library(tidyverse)

#### Load data ####
dat <- readRDS(here("model_fits", "survdat.rds"))

# set data up for maxLik package
x <- dat


#### Mixed Exponential No Covs ####
# exponential mixture model of 2 classes of relationship length
# (similar setup to Model 4 of Kirk's dissertation chapter 2)

logit = function(p){
  log(p/(1-p))
}

invLogit = function(phi){
  1/(1+exp(-phi))
}

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

mixExpFit = maxLik(
  logLik = mixML,
  start = mixExpStartParam,
  method = "BHHH"
)

summary(mixExpFit)
exp(coef(mixExpFit)[1:2])
invLogit(coef(mixExpFit)[3])
AIC(mixExpFit)

saveRDS(mixExpFit, "model_fits/mixed/mix-exp-nsfg.rds")


#### ego age + latent components SEPARATE MODELS ####
# separate models split by age cat 
# ego initial age category 
# 2 latent components (long/short rels)

# First assumed that age effects the lambda for each component in the same way, but got ests w/ inf SE
# Below contains 2 parameters for age effects by short rels and long rels, same
# might just be that ego age really does not explain a lot -- certainly doesn't help straight exponential much
# returning to this later 

logit = function(p){
  log(p/(1-p))
}

invLogit = function(phi){
  1/(1+exp(-phi))
}


eaMixML <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight1 = param[3]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  ifelse(censored == 1,
         yes = log(
           (dexp(closed, exp(logLambda1)) * invLogit(logitWeight1) +
              dexp(closed, exp(logLambda2)) * (1 - invLogit(logitWeight1)))  / 
             
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

x <- dat %>% filter(e.agecat %in% "15-19")
mix1519 <-  maxLik(logLik = eaMixML, start = mixExpStartParam, method = "BHHH")

summary(mix1519)
exp(coef(mix1519)[c(1:2)])
invLogit(coef(mix1519)[3])
AIC(mix1519 )
saveRDS(mix1519, "model_fits/mixed/mix1519.rds")


x <- dat %>% filter(e.agecat %in% "20-24")
mix2024 <-  maxLik(logLik = eaMixML, start = mixExpStartParam, method = "BHHH")

summary(mix2024)
1/exp(coef(mix2024)[1]); 1/exp(coef(mix2024)[2])
invLogit(coef(mix2024)[3])
AIC(mix2024)
saveRDS(mix2024, "model_fits/mixed/mix2024.rds")

x <- dat %>% filter(e.agecat %in% "25-29")
mix2529 <-  maxLik(logLik = eaMixML, start = mixExpStartParam, method = "BHHH")

summary(mix2529)
1/exp(coef(mix2529)[1]); 1/exp(coef(mix2529)[2])
invLogit(coef(mix2529)[3])
AIC(mix2529)
saveRDS(mix2529, "model_fits/mixed/mix2529.rds")

x <- dat %>% filter(e.agecat %in% "30-34")
mix3034 <-  maxLik(logLik = eaMixML, start = mixExpStartParam, method = "BHHH")

summary(mix3034)
1/exp(coef(mix3034)[1]); 1/exp(coef(mix3034)[2])
invLogit(coef(mix3034)[3])
AIC(mix3034)
saveRDS(mix3034, "model_fits/mixed/mix3034.rds")

x <- dat %>% filter(e.agecat %in% "35-39")
mix3539 <-  maxLik(logLik = eaMixML, start = mixExpStartParam, method = "BHHH")

summary(mix3539)
1/exp(coef(mix3539)[1]); 1/exp(coef(mix3539)[2])
invLogit(coef(mix3539)[3])
AIC(mix3539)
saveRDS(mix3539, "model_fits/mixed/mix3539.rds")

x <- dat %>% filter(e.agecat %in% "40-44")
mix4044 <-  maxLik(logLik = eaMixML, start = mixExpStartParam, method = "BHHH")

summary(mix4044)
1/exp(coef(mix4044)[1]); 1/exp(coef(mix4044)[2])
invLogit(coef(mix4044)[3])
AIC(mix4044)
saveRDS(mix4044, "model_fits/mixed/mix4044.rds")



#### mixed exp with weights by age cat #################

eaMixMLweights <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight15 = param[3]
  logitWeight20 = param[4]
  logitWeight25 = param[5]
  logitWeight30 = param[6]
  logitWeight35 = param[7]
  logitWeight40 = param[8]
  
  agecat15 <- x$e.15
  agecat20 <- x$e.20
  agecat25 <- x$e.25
  agecat30 <- x$e.30
  agecat35 <- x$e.35
  agecat40 <- x$e.40
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  weights <- logitWeight15*agecat15 + logitWeight20*agecat20 + logitWeight25*agecat25 +
    logitWeight30*agecat30 + logitWeight35*agecat35 + logitWeight40*agecat40
  
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
  logitWeight25 = logit(0.5),
  logitWeight30 = logit(0.5),
  logitWeight35 = logit(0.5),
  logitWeight40 = logit(0.5)
)

x <- dat
weightsMixExp <-  maxLik(logLik = eaMixMLweights, start = mixExpStartParam, method = "BHHH")

summary(weightsMixExp)
saveRDS(weightsMixExp, "model_fits/mixed/weightsMixExp.rds")



