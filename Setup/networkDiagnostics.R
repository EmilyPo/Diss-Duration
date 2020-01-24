# network diagnostics before epimodel simulation 
library(EpiModel)
#### load fitted models 
# first the models without the deg.other/deg.cohab terms because 
# that's what I'm currently having trouble with

fit.marcoh <- readRDS("~/Documents/Dissertation/R/Duration/Setup/fits/final/fit.marcoh.RDS")
fit.other <- readRDS("~/Documents/Dissertation/R/Duration/Setup/fits/final/fit.other.RDS")


mcmc.diagnostics(fit.marcoh)


plot(gof(fit.marcoh)) # issues here, why does model only need 7 params? 
# or gof(ergmFitObject, GOF=~model)


# plot(gof(fit.marcoh, GOF = "duration"))

