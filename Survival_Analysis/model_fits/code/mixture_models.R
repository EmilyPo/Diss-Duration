library(maxLik)
library(tidyverse)

#### Load data ####
dat <- readRDS("~/NSFG_DATA/Objects/altersegos_survdat.rds")

dat <- dat %>%
  mutate(e.race = as.factor(e.race),
         reltype = as.factor(reltype),
         e.partsyr3 = as.factor(e.partsyr3),
         e.osnpyr3 = as.factor(e.osnpyr3),
         e.maxospyr3 = as.factor(e.maxospyr3), 
         network1 = as.factor(network1), 
         e.agecat = as.factor(e.agecat),
         e.agecat.initial = as.factor(e.agecat.initial),
         e.deg.main = as.factor(e.deg.main)) 

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


######### WITH IMMUNE FRACTION ################

immuneMixML <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight1 = param[3]

  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  frac_immune <- 0.057
  
  ifelse(censored == 1,
         yes = log(
           (1-frac_immune) * 
           (dexp(closed, exp(logLambda1)) * invLogit(logitWeight1) +
              dexp(closed, exp(logLambda2)) * (1 - invLogit(logitWeight1)))  / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
         + frac_immune
         ),
         no = log(
           (1-frac_immune) * 
           (pexp(closed, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
              pexp(closed, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1))) / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
           + frac_immune
         )
  )
  
}

mixExpStartParam = c(
  logLambda1 = log(0.2),
  logLambda2 = log(0.01),
  logitWeight1 = logit(0.5)
)


immuneMix <-  maxLik(logLik = immuneMixML, start = mixExpStartParam, method = "BHHH")

summary(immuneMix)
exp(coef(mix1519)[c(1:2)])
invLogit(coef(mix1519)[3])
AIC(mix1519 )


#------------------------- WITH SURVEY WEIGHTS ----------------------------------------------

# weighted age dynamics
# immune fraction 

library(maxLik)
library(here)
library(tidyverse)

#---------------------------------------------------------
#### survey weights ####
#---------------------------------------------------------
#### Data ####
dat <- readRDS("~/NSFG_DATA/Objects/altersegos_survdat.rds")

dat <- dat %>%
  mutate(e.race = as.factor(e.race),
         reltype = as.factor(reltype),
         e.partsyr3 = as.factor(e.partsyr3),
         e.osnpyr3 = as.factor(e.osnpyr3),
         e.maxospyr3 = as.factor(e.maxospyr3), 
         network1 = as.factor(network1), 
         e.agecat = as.factor(e.agecat),
         e.agecat.initial = as.factor(e.agecat.initial),
         e.deg.main = as.factor(e.deg.main)) 

# set data up for maxLik package
x <- dat

logit = function(p){
  log(p/(1-p))
}

invLogit = function(phi){
  1/(1+exp(-phi))
}

#### tester - exponential, no mixture ####
# to make sure same results as flexsurv (aka using survey weights right)
# this works, except AIC is somehow different (but same coef estimate)

expMLweights <- function(param){
  lambda <- param
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  sw <- x$weight
  sum <- sum(x$weight)
  sw <- sw/sum
  
  ifelse(
    test = censored == 1, 
    yes = # if exact dur
      sw * 
      log(
        dexp(closed, lambda, log=F) /
          pexp(open, lambda, log=F, lower.tail = F)
      ), 
    # if right censored
    no = 
      sw* 
      log(
        pexp(closed, lambda, log=F, lower.tail = F)/
          pexp(open, lambda, log=F, lower.tail = F)
      )
  )
}

expFitweighted <- maxLik(logLik = expMLweights,  start=c(scale=.02), method="BHHH", control = list(iterlim = 1000))

summary(expFitweighted)
exp(coef(expFitweighted))
AIC(expFitweighted)


plot(expFit)
lines(x = durVec,
      y= pexp(q = durVec, coef(expFitweighted)[1], lower.tail = FALSE), col = "blue")

saveRDS(expFitweighted, here("model_fits", "nocovs", "exp-nsfg-weighted.rds"))


#### Mixed Exponential No Covs ####

mixMLweighted <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight1 = param[3]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  sw <- x$weight
  sum <- sum(x$weight)
  sw <- sw/sum
  
  ifelse(censored == 1,
         yes = sw * log(
           (dexp(closed, exp(logLambda1)) * invLogit(logitWeight1) +
              dexp(closed, exp(logLambda2)) * (1 - invLogit(logitWeight1))) / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
         ),
         no = sw * log(
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

mixExpFitweighted = maxLik(
  logLik = mixMLweighted,
  start = mixExpStartParam,
  method = "BHHH"
)

summary(mixExpFitweighted)
exp(coef(mixExpFitweighted)[1:2])
invLogit(coef(mixExpFitweighted)[3])
AIC(mixExpFitweighted)

saveRDS(mixExpFitweighted, here("model_fits", "mixed", "mix-exp-nsfg-weighted.rds"))


#### Mixed Exponential fit by age cat ####

mixMLweighted <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight1 = param[3]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  sw <- x$weight
  sum <- sum(x$weight)
  sw <- sw/sum
  
  ifelse(censored == 1,
         yes = sw * log(
           (dexp(closed, exp(logLambda1)) * invLogit(logitWeight1) +
              dexp(closed, exp(logLambda2)) * (1 - invLogit(logitWeight1))) / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
         ),
         no = sw * log(
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
mix1519weighted <-  maxLik(logLik = mixMLweighted, start = mixExpStartParam, method = "BHHH")

summary(mix1519weighted)
exp(coef(mix1519weighted)[c(1:2)])
invLogit(coef(mix1519weighted)[3])
saveRDS(mix1519weighted, "model_fits/mixed/mix1519weighted.rds")


x <- dat %>% filter(e.agecat %in% "20-24")
mix2024weighted <-  maxLik(logLik = mixMLweighted, start = mixExpStartParam, method = "BHHH")

summary(mix2024weighted)
1/exp(coef(mix2024weighted)[1]); 1/exp(coef(mix2024weighted)[2])
invLogit(coef(mix2024weighted)[3])
saveRDS(mix2024weighted, "model_fits/mixed/mix2024weighted.rds")

x <- dat %>% filter(e.agecat %in% "25-29")
mix2529weighted <-  maxLik(logLik = mixMLweighted, start = mixExpStartParam, method = "BHHH")

summary(mix2529weighted)
1/exp(coef(mix2529weighted)[1]); 1/exp(coef(mix2529weighted)[2])
invLogit(coef(mix2529weighted)[3])
saveRDS(mix2529weighted, "model_fits/mixed/mix2529weighted.rds")

x <- dat %>% filter(e.agecat %in% "30-34")
mix3034weighted <-  maxLik(logLik = mixMLweighted, start = mixExpStartParam, method = "BHHH")

summary(mix3034weighted)
1/exp(coef(mix3034weighted)[1]); 1/exp(coef(mix3034weighted)[2])
invLogit(coef(mix3034weighted)[3])
saveRDS(mix3034weighted, "model_fits/mixed/mix3034weighted.rds")

x <- dat %>% filter(e.agecat %in% "35-39")
mix3539weighted <-  maxLik(logLik = mixMLweighted, start = mixExpStartParam, method = "BHHH")

summary(mix3539weighted)
1/exp(coef(mix3539weighted)[1]); 1/exp(coef(mix3539weighted)[2])
invLogit(coef(mix3539weighted)[3])
saveRDS(mix3539weighted, "model_fits/mixed/mix3539weighted.rds")

x <- dat %>% filter(e.agecat %in% "40-44")
mix4044weighted <-  maxLik(logLik = mixMLweighted, start = mixExpStartParam, method = "BHHH")

summary(mix4044weighted)
1/exp(coef(mix4044weighted)[1]); 1/exp(coef(mix4044weighted)[2])
invLogit(coef(mix4044weighted)[3])
saveRDS(mix4044weighted, "model_fits/mixed/mix4044weighted.rds")



#### Mixed Exponential w/ 2 params & weights by age cat 

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
  
  sw <- x$weight
  sum <- sum(x$weight)
  sw <- sw/sum
  
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
         yes = sw * log(
           (dexp(closed, exp(logLambda1)) * invLogit(weights) +
              dexp(closed, exp(logLambda2)) * (1 - invLogit(weights)))  / 
             
             (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(weights) +
                pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(weights)))
         ),
         no = sw * log(
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
saveRDS(weightsMixExp, "model_fits/mixed/weightsMixExp_surveyweighted.rds")



#---------------------------------------------------------
#### immune fraction ####
#---------------------------------------------------------

############## exponential, no covs ##############

#### Exponential #### 

# NO LEFT TRUNCATION CORRECTION HERE)

expMLimmune <- function(param){
  lambda <- param[1]
  p <- param[2] 
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  ifelse(
    test = censored == 1, 
    
    yes = # if exact dur
      log(
        invLogit(p) * lambda * exp(-lambda*closed) 
      ), 
    # if right censored
    no = 
      log(
        1 - invLogit(p) + invLogit(p) * exp(-lambda*closed)
      )
  )
}


expImmuneStartParams <- c(
  lambda = 1,
  p = logit(0.7)
)

expFitImmune <- maxLik(
  logLik = expMLimmune,  
  start=expImmuneStartParams, 
  method="BHHH", control = list(iterlim = 1000))


summary(expFitImmune)
exp(coef(expFitImmune))
AIC(expFitImmune)

#saveRDS(expFit, "model_fits/exp-nsfg.rds")

#### try out using flexcure

fit <- curereg(Surv(t_o, t_c, censored)~1, cureformula = ~1, data = dat,
               timedist = "exp", ncausedist = "bernoulli")

fit2 <- curereg(Surv(t_o, t_c, censored) ~ reltype, cureformula = ~reltype, data = dat,
                timedist = "exp", ncausedist = "bernoulli")


################ immune with mixture (not working currently) ################
immuneMixML <- function(param){
  logLambda1 = param[1]
  logLambda2 = param[2]
  logitWeight1 = param[3]
  frac_immune = param[4]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  #frac_immune <- nrow(dat[which(dat$t_o > 0 & dat$censored==0),]) / nrow(dat)
  
  ifelse(censored == 1,
         yes = log(1-invLogit(frac_immune)) +
           log(
             (dexp(closed, exp(logLambda1)) * invLogit(logitWeight1) +
                dexp(closed, exp(logLambda2)) * (1 - invLogit(logitWeight1)))  / 
               
               (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                  pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
           ),
         no = invLogit(frac_immune) +
           log(1-invLogit(frac_immune)) + 
           log(
             (pexp(closed, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pexp(closed, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1))) / 
               
               (pexp(open, exp(logLambda1), lower.tail = FALSE) * invLogit(logitWeight1) +
                  pexp(open, exp(logLambda2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
           )
  )
  
}

mixImmuneStartParam = c(
  logLambda1 = log(0.2),
  logLambda2 = log(0.01),
  logitWeight1 = logit(0.5),
  logitFracImmune = logit(0.5)
)

immune_fracML = maxLik(
  logLik = immuneMixML,
  start = mixExpStartParam,
  method = "BHHH"
)

summary(immune_fracML)
exp(coef(immune_fracML)[1:2])
invLogit(coef(immune_fracML)[3])
AIC(immune_fracML)
