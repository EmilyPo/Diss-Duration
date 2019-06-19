
# hand coded models originally used in CSSS 544 final project, and beginnings of dissertation work 
# all except mixed exponentials replaced by flexsurv package (easier to incorporate weights & covariates)

library(maxLik)

########## No Covariates #########################
#### Exponential #### 
expML <- function(param){
  lambda <- param
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  ifelse(
    test = censored == 1, 
    yes = # if exact dur
      log(
        dexp(closed, lambda, log=F) /
          pexp(open, lambda, log=F, lower.tail = F)
      ), 
    # if right censored
    no = 
      log(
        pexp(closed, lambda, log=F, lower.tail = F)/
          pexp(open, lambda, log=F, lower.tail = F)
      )
  )
}

expFit<- maxLik(logLik = expML,  start=c(scale=1), method="BHHH", control = list(iterlim = 1000))

summary(expFit)
exp(coef(expFit))
AIC(expFit)

saveRDS(expFit, "model_fits/exp-nsfg.rds")

#### Mixed Exponential ####

# exponential mixture model of 2 classes of relationship length
# (similar setup to Model 4 of Kirk's dissertation chapter 2)

# NOTE THAT WE DO NOT USED LOG = TRUE HERE
# we need to add funcations up first and then log them (see parens)

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

saveRDS(mixExpFit, "model_fits/mix-exp-nsfg.rds")

#### Weibull ####

weibML <- function(param){
  shape <- param[1]
  scale <- param[2]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  ifelse(censored == 1, 
         # if exact dur
         yes = log(
           dweibull(closed, shape, scale, log=F) /
             pweibull(open, shape, scale, log.p=F, lower.tail = FALSE)
         ), 
         # if right censored
         no = log(
           pweibull(closed, shape, scale, log.p=F, lower.tail = FALSE) /
             pweibull(open, shape, scale, log.p=F, lower.tail = FALSE)
         )
  )
}

weibMaxLik <- maxLik(logLik = weibML,  start=c(shape=1, scale=1), method="BHHH")

summary(weibMaxLik)
exp(coef(weibMaxLik))
AIC(weibMaxLik)

saveRDS(weibMaxLik, "model_fits/weib-nsfg.rds")

#### Gamma ####

gammaML <- function(param){
  shape <- param[1]
  scale <- param[2]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  ifelse(censored == 1, 
         # if exact dur
         yes = log(
           dgamma(closed, shape=shape, scale=scale, log=F) /
             pgamma(open, shape=shape, scale=scale, log.p=F, lower.tail = F)
         ), 
         # if right censored
         no = log(
           pgamma(closed, shape=shape, scale=scale, log.p=F, lower.tail = F) /
             pgamma(open, shape=shape, scale=scale, log.p=F, lower.tail = F)
         )
  )
  
}

gammaMaxLik <- maxLik(logLik = gammaML,  start=c(shape=1, scale=1), method="BHHH")

summary(gammaMaxLik)
exp(coef(gammaMaxLik))
AIC(gammaMaxLik)

saveRDS(gammaMaxLik, "model_fits/gamma-nsfg.rds")

#### Mixed Models ####
mixMLweib <- function(param){
  logShape1 = param[1]
  logScale1 = param[2]
  logShape2 = param[3]
  logScale2 = param[4]
  logitWeight1 = param[5]
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  ifelse(censored == 1,
         yes = log(
           (dweibull(closed, exp(logShape1), exp(logScale1), log=F) * invLogit(logitWeight1) +
              dweibull(closed, exp(logShape2), exp(logScale2), log=F) * (1 - invLogit(logitWeight1))) / 
             
             (pweibull(open, exp(logShape1), exp(logScale1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pweibull(open, exp(logShape2), exp(logScale2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
         ),
         no = log(
           (pweibull(closed, exp(logShape1), exp(logScale1), lower.tail = FALSE) * invLogit(logitWeight1) +
              pweibull(closed, exp(logShape2), exp(logScale2), lower.tail = FALSE) * (1 - invLogit(logitWeight1))) / 
             
             (pweibull(open, exp(logShape1), exp(logScale1), lower.tail = FALSE) * invLogit(logitWeight1) +
                pweibull(open, exp(logShape2), exp(logScale2), lower.tail = FALSE) * (1 - invLogit(logitWeight1)))
         )
  )
}

mixWeibStartParam = c(
  logShape1 = log(0.2), 
  logScale1 = log(15),
  logShape2 = log(0.2),
  logScale2 = log(15),
  logitWeight1 = logit(0.7)
)

mixWeibFit = maxLik(
  logLik = mixMLweib,
  start = mixWeibStartParam,
  method = "BHHH"
)

summary(mixWeibFit)
#exp(coef(mixExpFit)[1:2])
#invLogit(coef(mixExpFit)[3])
#AIC(mixExpFit)






########## Age Dynamics #########################
#### ego age ####
# ego age at time of interview 
# no latent components
# parameter covariate model - regression

eaML <- function(param){
  lambda <- param[1]
  age_cov <- param[2]
  
  age <- x$e.age
  truncated <- x$truncated
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  cov_effects <- age_cov*age
  
  ifelse(censored == 1, 
         # if exact dur
         log(
           dexp(closed, lambda*exp(cov_effects), log=F) /
             pexp(open, lambda*exp(cov_effects), log=F, lower.tail = FALSE)
         ), 
         # if right censored
         log(
           pexp(closed, lambda*exp(cov_effects), log=F, lower.tail = FALSE) /
             pexp(open, lambda*exp(cov_effects), log=F, lower.tail = FALSE)
         )
  )
}


eaFit<- maxLik(logLik = eaML,  start=c(scale=1, age=0), method="BHHH", control = list(iterlim = 1000))

summary(eaFit)
exp(coef(eaFit)[2])
AIC(eaFit)

saveRDS(eaFit, "model_fits/egoAge-nsfg.rds")


#### #### ego age categories  ####
# ego age category at time of interview 
# no latent components
# several models here - comparing 1 age group to all others, attempting to fit all age cats results in INF std errors


ea.youngML <- function(param){
  lambda <- param[1]
  age_cov1529 <- param[2]
  
  cat15 <- x$e.15
  cat20 <- x$e.20
  cat25 <- x$e.25
  
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  effects15 <- age_cov1529*cat15
  effects20 <- age_cov1529*cat20
  effects25 <- age_cov1529*cat25
  
  covs <- lambda * exp(effects15 + effects20 + effects25)
  
  ifelse(censored == 1, 
         # if exact dur
         log(
           dexp(closed, covs, log=F) /
             pexp(open, covs, log.p=F, lower.tail = FALSE)
         ), 
         # if right censored
         log(
           pexp(closed, covs, log.p=F, lower.tail = FALSE) /
             pexp(open, covs, log.p=F, lower.tail = FALSE)
         )
  )
}

ea.youngStartParams <- c(
  lambda=1,
  ea1529=0
)

ea.youngFit <- maxLik(logLik = ea.youngML,  start=ea.youngStartParams, method="BHHH", 
                      control = list(iterlim = 1000))

summary(ea.youngFit)
AIC(ea.youngFit)

# not a great model fit, but better than exponential and suggests young age cat has roughly 1/2 of expected duration than the older set 

saveRDS(ea.youngFit, "model_fits/ea.youngFit-nsfg.rds")



#### ego age at beginning of rel ####
intial.ageML <- function(param){
  lambda <- param[1]
  age <- param[2]
  
  initial.age <- x$e.age.initial
  censored <- x$censored
  open <- x$t_o
  closed <- x$t_c
  
  age_effect <- age*initial.age
  
  ifelse(censored == 1, 
         # if exact dur
         log(
           dexp(closed, lambda*exp(age_effect), log=F) /
             pexp(open, lambda*exp(age_effect), log=F, lower.tail = FALSE)
         ), 
         # if right censored
         log(
           pexp(closed, lambda*exp(age_effect), log=F, lower.tail = FALSE) /
             pexp(open, lambda*exp(age_effect), log=F, lower.tail = FALSE)
         )
  )
  
} 

initial.eaFit <- maxLik(logLik = intial.ageML,  start=c(lambda=1, age=0), method="BHHH", control = list(iterlim = 1000))

summary(initial.eaFit)
exp(coef(initial.eaFit)[2])
AIC(initial.eaFit)

saveRDS(initial.eaFit, "model_fits/initialAge-nsfg.rds")

#### plot - old ####

# THIS DOES NOT CAPTURE CONTINUOUS VARS -- HOW TO VISUALALLY COMPARE TO K-M? Possble??
km <- readRDS("model_fits/km-nsfg.rds")

durVec <- 0:max(x$t_c)

plot(
  km, conf.int = FALSE,
  xlab = "t", ylab = "S(t)",
  main="Ego Age, Exponential Dists"
)
points(
  x = x$t_c,
  y = pexp(x$t_c, coef(eaFit)[1] * exp(coef(eaFit)[2]*x$e.age), lower.tail = FALSE),
  col = alpha("darkred", 0.1), lwd = 2
)

lines(
  x = durVec,
  y = pexp(durVec, coef(initial.eaFit)[1] * exp(coef(initial.eaFit)[2]), lower.tail = FALSE),
  col = "steelblue3", lwd = 2
)
legend(300,1, legend = c("Kaplan-Meier", "Current Age", "Initial Age"),
       col = c("black", "darkred", "steelblue3"), cex=0.7, lwd=1)


ea.surv <- sort((pexp(x$t_c, coef(eaFit)[1] * exp(coef(eaFit)[2]*x$e.age), lower.tail = FALSE))) 
ea.surv$t <- x$t_c
head(ea.surv)
plot(x=ea.surv[,2], y=ea.surv[,1])


