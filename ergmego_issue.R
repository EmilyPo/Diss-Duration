# ergm.ego issue with nodefactor for ego-only attribute 

# uninstalled cran version, installed github dev version, same error: 

library(ergm.ego)
session_info()

# data
dat <- readRDS("~/NSFG_DATA/Objects/fullEgodata.rds")

#all alters
egodat <- egodata(egos=dat$egos, alters=dat$altersAllActive, egoWt = dat$egos$weight, egoIDcol = "ego")

# marriage/cohab alters
egodat_marcoh <- egodata(egos=dat$egos, alters=dat$altersMarCoh, egoWt = dat$egos$weight, egoIDcol = "ego")

fit.marcoh.nodecov <- ergm.ego(egodat_marcoh ~ edges + 
                                 nodecov("age") + 
                                 absdiff("sqrtage") + 
                                 nodefactor("deg.other") +
                                 offset(nodematch("sex", diff = FALSE)) +
                                 offset("concurrent"), 
                               offset.coef = c(-Inf, -Inf),
                               control = control.ergm.ego(ppopsize = 5000))


# "There appears to be a mismatch between estimated statistic and the sufficient statistic of the ERGM: 
# statistics ‘nodefactor.deg.other.2’ estimated from data are extraneous to the ERGM. A common cause of 
# this is that egos and alters do not have a consistent set of levels for one or more factors."

