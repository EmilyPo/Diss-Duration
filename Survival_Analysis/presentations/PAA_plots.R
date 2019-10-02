

##### Plots for PAA Extended Abstract #####
library(here)

# load models 
expFit <- readRDS(here("model_fits", "nocovs", "exp.rds")) 
km_weighted <- readRDS(here("model_fits", "km", "km-weighted.rds"))
weibFit <- readRDS(here("model_fits", "nocovs", "weib.rds"))
mixExpFitweighted <- readRDS(here("model_fits", "mixed", "mix-exp-nsfg-weighted.rds"))
e.agecat <- readRDS(here("model_fits", "age", "expAgeCatCurrent.rds"))
km_agecat_weighted <- readRDS(here("model_fits",  "km", "km_ac_weighted.rds"))
e.reltype <- readRDS(here("model_fits", "covs", "ereltype.rds")) 
e.network <- readRDS(here("model_fits", "covs", "enetwork.rds"))
gammaFit <- readRDS(here("model_fits", "nocovs", "gamma.rds"))

# functions
logit = function(p){
  log(p/(1-p))
}

invLogit = function(phi){
  1/(1+exp(-phi))
}

durVec <- c(1:max(km_weighted$time))

s <- summary(e.agecat)

par(mfrow=c(2,2))
# General Plot, no covariates
plot(km_weighted, 
     conf.int = FALSE,
     xlab = "t, Months", ylab = "S(t)",
     main="Survival of Relationships, No Covariates"
)
lines(expFit, col="blue", lwd=2)
lines(weibFit, col="red", lwd=2)
lines(gammaFit, col="darkgreen", lwd=2)
#lines(x = durVec,
#      y = invLogit(coef(mixExpFitweighted)[3]) *
#        pexp(q = durVec, exp(coef(mixExpFitweighted)[1]), lower.tail = FALSE) +
#       (1-invLogit(coef(mixExpFitweighted)[3])) *
#        pexp(q = durVec, exp(coef(mixExpFitweighted)[2]), lower.tail = FALSE),
#      col ="purple", lwd=2)

legend("topright", legend = c("K-M Reference", "Exponential", "Gamma", "Weibull"),
       col = c("black", "blue", "darkgreen", "red"), cex=0.9, lwd=2, lty = c(1,1,1,1))



# Age Category 
plot(km_agecat_weighted, lty=2, col = c(2:6,"orange"), lwd=2,
     ylab="S(t)", xlab="t, Months", 
     main="Current Age Category (Exponential Dist)")
lines(s$`e.agecat=15-19`$est, type = "lines", col = 2, lwd=1.5)
lines(s$`e.agecat=20-24`$est, type = "lines", col = 3, lwd=1.5)
lines(s$`e.agecat=25-29`$est, type = "lines", col = 4, lwd=1.5)
lines(s$`e.agecat=30-34`$est, type = "lines", col = 5, lwd=1.5)
lines(s$`e.agecat=35-39`$est, type = "lines", col = 6, lwd=1.5)
lines(s$`e.agecat=40-44`$est, type = "lines", col = "orange", lwd=1.5)
legend("topright", lty=c(2,1,1,1,1,1,1,1), lwd=c(rep(2,8)), col=c(1,1,2,3,4,5,6,"orange"), 
       legend = c("K-M", "Fitted Model", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44"), cex = 0.7)

# Relationships: Two Types (Exponential Dist)
plot(e.network, main = "Relationships: 2 Types (Exponential Dist)", 
     col=c("red", "blue"), xlab="t, Months", ylab = "S(t)")
legend("topright", legend = c("K-M Reference", "Marriage/Cohabitation", "Other"), 
       col=c("black", "red","blue"), lwd=2, lty=1, cex=.9)

# Relationships: Three Types (Exponential Dist)
plot(e.reltype, main = "Relationships: Three Types (Exponential Dist)", 
     col=c("red", "blue", "springgreen3"), xlab="t, Months", ylab = "S(t)")
legend("right", legend = c("K-M Reference", "Marriage", "Cohab", "Other"), 
       col=c("black", "red", "springgreen3", "blue"), lwd=2, lty=1, cex=.7)
