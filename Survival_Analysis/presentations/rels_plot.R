library(tidyverse)

#### load nsfg -- not done here ###
#### set up data for plots #### 
# observed rels - first 75 listed in data 
examples <- dat[1:75,]
examples <- examples %>% 
              select(dfs, dls, active) %>% # select first date of sex, last date of sex, and if rel is active
              mutate(dfs=-dfs, dls=-dls) %>% # make dates negative so date of interview = time 0
              mutate(dfs_c = dls, dls_c=-dfs) # flag time rels may be right censored and make hypothetical "end" date 

examples$id <- c(1:nrow(examples))

# left truncated rels - based on next 75 observed in data
second <- dat[76:150,] %>% select(dfs, dls, active) %>% mutate(dfs=-dfs, dls=-dls)
# randomly assign offset to visually place rels further back in time 
# (we don't see these rels in the data because they ended at least 1 year prior to interview)
offset <- sample(12:30, 75, replace=T) 
second <- cbind(second,offset)
second <- second %>% mutate(dfs=dfs-offset, dls=dls-offset)
second <- second[1:37,]
second$id <- seq(1,74, by=2)

# right censored rels (rels active on day of interview)
# randomly assign how long the rel will continue
rights <- examples %>% filter(dls==0 & active==1)
rights$dls_c <- sample(1:80, nrow(rights), replace=T)

#### plot ####
plot(x=NA,
     xlim=c(-350, 60), xaxt='na', xlab=NA,
     ylim=c(0, 75), yaxt='n', ylab=NA,
     main = "Example: Known and Unknown Relationships in NSFG"
)
# observation window (time 0 to 12 months previous)
abline(v=-12, lty=2)
abline(v=0, lty=2)
# observed rels
segments(
  x0=examples$dls,
  x1=examples$dfs,
  y0=order(examples$id), 
  col="darkblue"
)
# truncated (unobserved) rels 
segments(
  x0=second$dls,
  x1=second$dfs,
  y0=second$id,
  col="darkgreen"
)
# right-censored examples
segments(
  x0=rights$dfs_c,
  x1=rights$dls_c,
  y0=rights$id, 
  col="darkred"
)
# legend
legend(-348, 17, legend = c("Known Durations", "Possible Right Censored Durations", "Possible Unobserved (Truncated) Durations"), 
       col = c("darkblue", "darkred", "darkgreen"), lwd=2, cex=0.6)