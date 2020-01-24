#' @title Sample Relationships with NSFG Design 
#'
#' @description This function takes the output from simulation cel.temp and cel.complete, 
#' removes relationships that ended more than 12 months from last time step,
#' adds exact duration flag and limits to 3 relationships per ID. 
#' 
#' this function currenly only works for 1 simulation?
#'              
#' @param x is the netsim output
#' 
#' @keywords 
#'
#' @export
#'

as_nsfg_rels <- function(x){

  nsteps <- x$control$nsteps
  time.step <- x$param$time.unit
  nsims <- x$control$nsims
  
  allRelsComplete <- NULL
  
  for (i in 1:nsims){
  # rels that ended
  rels <- x$cel.complete[[i]] 
  
  # rels that are ongoing
  relscur <- x$cel.temp[[i]]
  
  # full 
  allRels <- rbind(rels, relscur)
  
  # make all end == NA 100 and create censored flag (1 if exact, 0 if ongoing) 
  # (previous NA, newly 100 are all ongoing)
  allRels[,"exact"] <- ifelse(is.na(allRels[,"end"]), 0, 1) 
  allRels[,"end"] <- ifelse(is.na(allRels[,"end"]), nsteps, allRels[,"end"]) 
  allRels[,"len"] <- allRels[,"end"]-allRels[,"start"]
  
  # limit to those that ended in last 12 months or are ongoing
  cap <- nsteps - floor(365/time.step)
  
  current <- which(allRels[,"end"] >= cap)
  
  allRels <- allRels[current,]
  
    # limit that to 3 rels per UID
    tabP1 <- table(allRels$p1)
    tabP2 <- table(allRels$p2)
    
    relsp1 <- NULL
    relsp2 <- NULL
    
    if (length(which(tabP1>3)) > 0) {
      #grab ids that have more than 3 partners in last year
      ids <- as.numeric(names(which(tabP1>3))) 
      
      # drop all rels associated with that id from main dat and make separate dataframe
      relsToLimit <- allRels[which(allRels$p1 %in% ids),]
      Rels <- allRels[-which(allRels$p1 %in% ids),]
      
      # for each id grab 3 most recent rels
      for (i in length(ids)) {
        z <- relsToLimit[relsToLimit$p1 %in% ids[i],]
        mostRecent <- nrow(z)
        third <- mostRecent-2
        z <- z[third:mostRecent,]
        relsp1 <- rbind(relsp1, z)
      }
    }
    
    if (length(which(tabP2>3)) > 0) {
      #grab ids that have more than 3 partners in last year
      ids <- as.numeric(names(which(tabP2>3))) 
      
      # drop all rels associated with that id from main dat and make separate dataframe
      relsToLimit <- allRels[which(allRels$p2 %in% ids),]
      Rels <- allRels[-which(allRels$p2 %in% ids),]
      
      # for each id grab 3 most recent rels
      for (i in length(ids)) {
        z <- relsToLimit[relsToLimit$p2 %in% ids[i],]
        mostRecent <- nrow(z)
        third <- mostRecent-2
        z <- z[third:mostRecent,]
        relsp2 <- rbind(relsp2, z)
      }
    }
  
  relsLimited <- rbind(relsp1, relsp2)
  
  allRelsThisSim <- rbind(allRels, relsLimited)
  
  allRelsComplete <- rbind(allRelsComplete, allRelsThisSim)
  
  }
  
  return(allRelsComplete)
  
}




