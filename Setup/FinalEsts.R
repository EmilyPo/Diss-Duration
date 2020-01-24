# script for final network estimation

# ---------- setup -------------
# libs
library(ergm.ego)
library(EpiModel)

# data
dat <- readRDS("~/NSFG_DATA/Objects/fullEgodata.rds")

# adjust attributes to what epimodelHIV will expect
# sex = male 1 or 0
dat$egos$male <- ifelse(dat$egos$sex %in% "F", 0, 1)
dat$altersOther$male <- ifelse(dat$altersOther$sex %in% "F", 0, 1)
dat$altersMarCoh$male <- ifelse(dat$altersMarCoh$sex %in% "F", 0, 1)
# debuted = HADSEX
dat$egos$debuted <- ifelse(dat$egos$HADSEX==1, 1, 0)

# agecat, 1:6
breaks <- c(15, 20, 25, 30, 35, 40, 45)
dat$egos$agecat <- cut(dat$egos$age, breaks, labels=F, right=F)
dat$altersOther$agecat <- cut(dat$altersOther$age, breaks, labels=F, right=F)
dat$altersMarCoh$agecat <- cut(dat$altersMarCoh$age, breaks, labels=F, right=F)

# make deg.other & deg.marcoh binary:
# 0 if 0 "other" parters, 1 if >= 1. 
# same for marcoh - there are 47 people who report active mar/cohabs > 1. 
# length(dat$egos[dat$egos$deg.marcoh > 1])
dat$egos$deg.other.binary <- ifelse(dat$egos$deg.other == 0, 0, 1)
dat$egos$deg.marcoh.binary <- ifelse(dat$egos$deg.marcoh == 0, 0, 1)

# reduce attributes down to what's actually used in the ergm.ego estimation
alter_attrs <- c("ego", "age", "weight","male", "agecat", "race", "sqrtage")
egos_attrs <- c("ego", "weight", "age", "male", "agecat", "sqrtage", "race", "debuted", "deg.marcoh.binary", "deg.other.binary")
dat$egos <- dat$egos[,egos_attrs]
dat$altersMarCoh <- dat$altersMarCoh[,alter_attrs]
dat$altersOther <- dat$altersOther[,alter_attrs]

#all alters
#egodat <- egodata(egos=dat$egos, alters=dat$altersAllActive, egoWt = dat$egos$weight, egoIDcol = "ego")

# marriage/cohab alters
egodat_marcoh <- egodata(egos=dat$egos, alters=dat$altersMarCoh, egoWt = dat$egos$weight, egoIDcol = "ego")

# other alters 
egodat_other <- egodata(egos=dat$egos, alters=dat$altersOther, egoWt = dat$egos$weight, egoIDcol = "ego")

# ------------ Estimation for Two-Network Scenario -------------------

ppopsize=15000

##### Marriage & Cohab Network ######
fit.marcoh <- ergm.ego(egodat_marcoh ~ edges + 
                                nodecov("age") + 
                                absdiff("sqrtage") + 
                                nodefactor("deg.other.binary") + 
                                offset(nodematch("male", diff = FALSE)) +
                                offset(nodefactor("debuted", levels=1)) + 
                                offset("concurrent"), 
                                offset.coef = c(-Inf, -Inf, -Inf),
                                control = control.ergm.ego(ppopsize = ppopsize, 
                                                          ergm.control = control.ergm(MCMLE.maxit = 100,
                                                                                      MCMC.interval = 1024*10)))

saveRDS(fit.marcoh, "~/Documents/Dissertation/R/Duration/Setup/fits/final/fit.marcoh.RDS")

#fit.marcoh <- readRDS("~/Documents/Dissertation/R/Duration/Setup/fits/final/fit.marcoh.RDS")

summary(fit.marcoh)

mcmc.diagnostics(fit.marcoh)

plot(gof(fit.marcoh, GOF = "model"))



##### Other Network ######


# what about degree in other network? 

fit.other <- ergm.ego(egodat_other ~ edges + 
                          nodecov("age") + 
                          absdiff("sqrtage") + 
                          nodefactor("deg.marcoh.binary") +
                          concurrent +
                          offset(nodematch("male", diff = FALSE)) +
                          offset(nodefactor("debuted", levels=1)),
                          offset.coef = c(-Inf, -Inf),
                          control = control.ergm.ego(ppopsize = ppopsize, 
                                                     ergm.control = control.ergm(MCMLE.maxit = 100,
                                                                                 MCMC.interval = 1024*10)))

saveRDS(fit.other, "~/Documents/Dissertation/R/Duration/Setup/fits/final/fit.other.RDS")

#fit.other.cov <- readRDS("~/Documents/Dissertation/R/Duration/Setup/fits/final/other-nodecov.RDS")

summary(fit.other)

mcmc.diagnostics(fit.other)

plot(gof(fit.other, GOF = "model"))

plot(gof(fit.other, GOF = "degree"))
  