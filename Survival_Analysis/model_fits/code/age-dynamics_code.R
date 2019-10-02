# code for covariate models 

library(tidyverse)
library(here)
library(flexsurv)


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

# ----------------------- RACE -------------------------------------------------------

e.race <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race, 
                  data = dat, dist = "exp", weights = e.weight) 

w.race <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race, 
                      data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.race, "model_fits/covs/erace.rds")
saveRDS(w.race, "model_fits/covs/wrace.rds")

# ----------------------- RELTYPE ----------------------------------------------------

e.reltype <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype, 
                      data = dat, dist = "exp", weights = e.weight) 

saveRDS(e.reltype, "model_fits/covs/ereltype.rds")

w.reltype <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype, 
                         data = dat, dist = "weibull", weights = e.weight) 

saveRDS(w.reltype, "model_fits/covs/wreltype.rds")

# ----------------------- NPARTS LAST YEAR -------------------------------------------

e.partsyr <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.partsyr3, 
                         data = dat, dist = "exp", weights = e.weight) 

saveRDS(e.partsyr, "model_fits/covs/epartsyr.rds")

w.partsyr <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.partsyr3, 
                         data = dat, dist = "weibull", weights = e.weight) 

saveRDS(w.partsyr, "model_fits/covs/wpartsyr.rds")

# compare with osnpyr & maxospyr (all opposite sex partners, not just vaginal sex partners)
# very little differences, but AIC for e.partsyr smaller and I am modeling vaginal intercourse transmitted sex anyway...

e.osnpyr <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.osnpyr3, 
                        data = dat, dist = "exp", weights = e.weight) 

e.maxospyr <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.maxospyr3, 
                        data = dat, dist = "exp", weights = e.weight) 


# ----------------------- MOMENTARY DEGREE -------------------------------------------

e.deg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main, 
                         data = dat, dist = "exp", weights = e.weight) 

w.deg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main, 
                         data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.deg, "model_fits/covs/edeg.rds")
saveRDS(w.deg, "model_fits/covs/wdeg.rds")


# ----------------------- network type -------------------------------------------

e.network <- flexsurvreg(Surv(t_o, t_c, censored) ~ network1, 
                       data = dat, dist = "exp", weights = e.weight) 


w.network <- flexsurvreg(Surv(t_o, t_c, censored) ~ network1, 
                       data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.network, "model_fits/covs/enetwork.rds")
saveRDS(w.network, "model_fits/covs/wnetwork.rds")

# ----------------------- MULTIVARIATE -------------------------------------------

e.relparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.partsyr3, 
                          data = dat, dist = "exp", weights = e.weight)

w.relparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.partsyr3, 
                          data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.relparts, "model_fits/covs/erelparts.rds")
saveRDS(w.relparts, "model_fits/covs/wrelparts.rds")

#---

e.ageparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.partsyr3, 
                          data = dat, dist = "exp", weights = e.weight) 

w.ageparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.partsyr3, 
                          data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.ageparts, "model_fits/covs/eageparts.rds")
saveRDS(w.ageparts, "model_fits/covs/wageparts.rds")

#---

e.raceparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.partsyr3, 
                          data = dat, dist = "exp", weights = e.weight)

w.raceparts <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.partsyr3, 
                          data = dat, dist = "weibull", weights = e.weight) 

saveRDS(e.raceparts, "model_fits/covs/eraceparts.rds")
saveRDS(w.raceparts, "model_fits/covs/wraceparts.rds")

#---

e.racerels <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.race, 
                          data = dat, dist = "exp", weights = e.weight)

w.racerels <- flexsurvreg(Surv(t_o, t_c, censored) ~ reltype + e.race, 
                          data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.racerels, "model_fits/covs/eracerels.rds")
saveRDS(w.racerels, "model_fits/covs/wracerels.rds")

#---

e.raceage <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.race, 
                          data = dat, dist = "exp", weights = e.weight)

w.raceage <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + e.race, 
                         data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.raceage, "model_fits/covs/eraceage.rds")
saveRDS(w.raceage, "model_fits/covs/wraceage.rds")

#---

e.agerels <-  flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + reltype, 
                       data = dat, dist = "exp", weights = e.weight)

w.agerels <-  flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + reltype, 
                          data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.agerels, "model_fits/covs/eagerels.rds")
saveRDS(w.agerels, "model_fits/covs/wagerels.rds")

#---

e.racedeg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.deg.main, 
                           data = dat, dist = "exp", weights = e.weight)

w.racedeg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.race + e.deg.main, 
                         data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.racedeg, "model_fits/covs/eracedeg.rds")
saveRDS(w.racedeg, "model_fits/covs/wracedeg.rds")

#---

e.agedeg <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main + e.agecat, 
                        data = dat, dist = "exp", weights = e.weight)

w.agedeg <-  flexsurvreg(Surv(t_o, t_c, censored) ~ e.deg.main + e.agecat, 
                          data = dat, dist = "weibull", weights = e.weight)

saveRDS(e.agedeg, "model_fits/covs/eagedeg.rds")
saveRDS(w.agedeg, "model_fits/covs/wagedeg.rds")



# ----------------------- AGE --------------------------------------------------------
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
eaCatFlex <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat, 
                         data = dat, dist = "exp", weights = e.weight) 

eaCatFlexWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat, 
                             data = dat, dist = "weibull", weights = e.weight) 

saveRDS(eaCatFlex, "model_fits/age/expAgeCatCurrent.rds")
saveRDS(eaCatFlexWeib, "model_fits/age/weibAgeCatCurrent.rds")

eaCatAgeDiff <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + diff.sqrt.age, 
                            data = dat, dist = "exp", weights = e.weight)

eaCatAgeDiffWeib <- flexsurvreg(Surv(t_o, t_c, censored) ~ e.agecat + diff.sqrt.age, 
                                data = dat, dist = "weibull", weights = e.weight) 

saveRDS(eaCatAgeDiff, "model_fits/age/expAgeDiffAgeCatCurrent.rds")
saveRDS(eaCatAgeDiffWeib, "model_fits/age/weibAgeDiffAgeCatCurrent.rds")

