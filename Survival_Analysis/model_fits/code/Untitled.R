# weighted age dynamics

library(survival)

#### Data ####
dat <- readRDS("~/nsfg_analysis/data/alters_egos.rds")

## remove valid skips, missings, once partners
dat <- dat[which(dat$active < 2),]
dat <- dat[which(dat$once==0),]

## new variables

# time
# t_o = relationship dur at 1 year prior to interview = edge_age_month-12, if edge_age_month-12 is negative, 0
# t_c = relationship month end (time = 1-12), if still active, 12
dat$t_o <- ifelse((dat$edge_age_month - 12 >= 0), dat$edge_age_month-12, 0)
dat$t_c <- dat$edge_age_month

# censoring variable (0 if censored, 1 if exact duration data) (active vs inactive)
dat$censored <- ifelse(dat$active==1, 0, 1)

# set up dummy vars for ego age cat
dat$e.15 <- ifelse(dat$e.agecat %in% "15-19", 1, 0)
dat$e.20 <- ifelse(dat$e.agecat %in% "20-24", 1, 0)
dat$e.25 <- ifelse(dat$e.agecat %in% "25-29", 1, 0)
dat$e.30 <- ifelse(dat$e.agecat %in% "30-34", 1, 0)
dat$e.35 <- ifelse(dat$e.agecat %in% "35-39", 1, 0)
dat$e.40 <- ifelse(dat$e.agecat %in% "40-44", 1, 0)

# set data up for maxLik package
x <- dat

# set up longer dat for weighted analysis in handcode MLE
datLonger <- dat[rep(seq_len(nrow(dat)), dat$e.weight), 1:ncol(dat)]

