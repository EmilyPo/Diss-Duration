# 
# this file takes the alters_egos rds from nsfg_analysis data folder, makes some adjustments
# and saves it in survival_analysis folder for ease of access
# this needs to be re-run if alters_egos gets updated
# last run: April 24 2019

# import
dat <- readRDS("~/nsfg_analysis/data/alters_egos.rds")

## remove valid skips, missings, once partners
dat <- dat[which(dat$active < 2),]
dat <- dat[which(dat$once==0),]
dat <- dat[which(dat$rel <= 3),]
dat$rel <- factor(dat$rel, 
                  labels = c("Spouse", "Cohab", "Other"))
dat$rel2 <- ifelse((dat$optype == 1 | dat$optype == 3), 1, 
            ifelse((dat$optype == 2 | dat$optype == 4), 2, 
            ifelse(dat$optype == 5, 3, NA)))
dat$rel2 <- factor(dat$rel2, 
                   labels = c("C/F Spouse", "C/F Cohab", "Other"))

## new variables

# time
# t_o = relationship dur at 1 year prior to interview = edge_age_month-12, if edge_age_month-12 is negative, 0
# t_c = relationship month end (time = 1-12), if still active, 12
dat$t_o <- ifelse((dat$edge_age_month - 12 >= 0), dat$edge_age_month-12, 0)
dat$t_c <- dat$edge_age_month

dat$t_o_years <- dat$t_o/12
dat$t_c_years <- dat$t_c/12

# censoring variable (0 if censored, 1 if exact duration data) (active vs inactive)
dat$censored <- ifelse(dat$active==1, 0, 1)

# set up dummy vars for ego age cat
dat$e.15 <- ifelse(dat$e.agecat %in% "15-19", 1, 0)
dat$e.20 <- ifelse(dat$e.agecat %in% "20-24", 1, 0)
dat$e.25 <- ifelse(dat$e.agecat %in% "25-29", 1, 0)
dat$e.30 <- ifelse(dat$e.agecat %in% "30-34", 1, 0)
dat$e.35 <- ifelse(dat$e.agecat %in% "35-39", 1, 0)
dat$e.40 <- ifelse(dat$e.agecat %in% "40-44", 1, 0)

# save out 
saveRDS(dat, "model_fits/survdat.rds")