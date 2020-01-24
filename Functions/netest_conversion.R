
#' @title Convert ergm.ego to netest
#'
#' @description Converts an ergm.ego output object to a netest object for use in EpiModel.
#' 
#' code developed by JKB for SHAMP
#'
#' @param x is a fit ergm.ego object.
#' @param y An object of class \code{disscoef} output from the dissolution_coefs function
#' @param set.formula.nw Logical indicating whether to add the (empty) starting network to the environment of the formation model (fit$formula). Default is TRUE, to mimic netest. Having a network attached to the formula facilitates using defaults for other ergm functions, e.g. see the help file for simulate.ergm, "basis" parameter.
#' @param start.net An object of class network that represents the starting network for the fit. May have edges (they will be deleted). Default is NULL, in which case it will be inferred from the fitted network. Use this option to standardize the starting network for ee.netest calls across different ergm.ego fits where the input data was the same but the pseudopopulation ended up with slight differences. See example.
#'
#'
#' @return
#' This function returns a \code{netest} for use with EpiModel functions.
#' @examples
#'
#' \dontrun{
#' # ergm.ego model
#' library(ergm.ego)
#' data(faux.mesa.high)
#' fmh.ego <- as.egodata(faux.mesa.high)
#' egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'                    +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"),
#'                    popsize=network.size(faux.mesa.high))
#' # This also works - use to test having an offset
#' egofit <- ergm.ego(fmh.ego~edges+offset(nodematch("Sex", diff=FALSE)),
#'                    offset.coef=c(-Inf),
#'                    popsize=network.size(faux.mesa.high))
#'
#' #----------------------------
#' # Single dissolution
#' library(EpiModel)
#' diss = ~offset(edges)
#' coef.diss <- dissolution_coefs(dissolution = diss,
#'                                duration = 400,
#'                                d.rate = 3.5e-5)
#'
#' est <- ee.netest(egofit, coef.diss)
#'
#' # Netdx, static
#' dx1 <- netdx(est, nsims = 1e4, dynamic = FALSE,
#'              nwstats.formula = ~edges + meandeg + concurrent)
#' dx1
#' plot(dx1, method = "b", stats = c("edges", "concurrent"))
#'
#' # netdx, dynamic
#' dx2 <- netdx(est, nsims = 5, nsteps = 500,
#'              nwstats.formula = ~edges + meandeg + concurrent,
#'              set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
#' dx2
#' plot(dx2, stats = c("edges", "meandeg"), plots.joined = FALSE)
#' plot(dx2, type = "duration")
#' plot(dx2, type = "dissolution", qnts.col = "orange2")
#' plot(dx2, type = "dissolution", method = "b", col = "bisque")
#'
#' 
#' #----------------------------
#' # Dissolution varies by Race
#' dissR <- ~offset(edges) + offset(nodematch("Race", diff=TRUE))
#' coef.dissR <- dissolution_coefs(dissolution = dissR,
#'                                 duration = c(300, 100, 200, 300, 100, 400))
#' estR <- ee.netest(egofit, coef.dissR)
#' dxR1 <- netdx(estR, nsims=1e4, dynamic=FALSE)
#' dxR2 <- netdx(estR, nsims=5, nsteps=500, sequential=FALSE,
#'              set.control.ergm=control.simulate.ergm(MCMC.burnin=1e6))
#' dxR2
#' # This gives the error "Dissolution plots for heterogeneous dissolution models not currently available"
#' # plot(dxR2, type='dissolution')
#'
#' #----------------------------
#' # Standardize starting network across fits (trivial example)
#' # To come: meaningful example showing utility for netsim() with 
#' # a list of networks
#' # Add weights to the egodata
#' fmh.ego$egoWt <- sample(1:100, nrow(fmh.ego$egos), replace=TRUE)
#' 
#' # Fit twice using ppop.wt='sample' to get sampling variability
#' fit1 <- ergm.ego(fmh.ego~edges+offset(nodematch("Sex", diff=FALSE)),
#'                    offset.coef=c(-Inf),
#'                    popsize=1000,
#'                    control=control.ergm.ego(ppop.wt='sample'))
#' fit2 <- ergm.ego(fmh.ego~edges+offset(nodematch("Sex", diff=FALSE)),
#'                    offset.coef=c(-Inf),
#'                    popsize=1000,
#'                    control=control.ergm.ego(ppop.wt='sample'))
#' 
#' # Show that the networks are not identical
#' table(fit1$network %v% 'Sex')
#' table(fit2$network %v% 'Sex')
#' 
#' # Standardize both ee.netest() objects to fit1's network
#' est1 <- ee.netest(fit1, coef.diss)
#' est2 <- ee.netest(fit2, coef.diss, start.net=fit1$network)
#' 
#' # Show that retrieved networks are idential
#' table(ergm.getnetwork(est1$fit$formula) %v% 'Sex')
#' table(ergm.getnetwork(est2$fit$formula) %v% 'Sex')
#' }
#'
#' @keywords module convert object class
#' @export


ee.netest <- function (x,y, set.formula.nw=TRUE, start.net=NULL) {
  
  ## the netest object contains
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  
  ##Initialize output list
  out <- list()
  
  ## fit is an ergm object that contains.
  ##[1] "coef"          "sample"        "iterations"    "MCMCtheta"     "loglikelihood" "gradient"      "hessian"       "covar"
  ##[9] "failure"       "network"       "newnetworks"   "newnetwork"    "coef.init"     "est.cov"       "coef.hist"     "stats.hist"
  ##[17] "steplen.hist"  "control"       "etamap"        "formula"       "target.stats"  "target.esteq"  "constraints"   "reference"
  ##[25] "estimate"      "offset"        "drop"          "estimable"     "null.lik"
  
  ##Delete egodata
  x$egodata<-NULL
  
  ##Assign the ergm.ego object to fit in netest
  out$fit <- x
  class(out$fit) <- "ergm"
  
  ##Adjust the ergm coeficient to get rid of the network size adjustment
  out$fit$coef<-x$coef[2:length(x$coef)]
  out$fit$coef[1]<-out$fit$coef[1] + x$coef[1]
  
  #Remove the netsize.adj colum from "sample"
  ncol<-length(x$sample[1,])
  out$fit$sample<-x$sample[,-1]
  
  #Iterations is copied directly from the ergm.ego fit to the netest fit object.
  out$fit$iterations <- x$iterations
  
  ##Adjust the MCMCtheta coeficient to get rid of the network size adjustment
  out$fit$MCMCtheta<-x$MCMCtheta[2:length(x$MCMCtheta)]
  out$fit$MCMCtheta[1]<-out$fit$MCMCtheta[1] + x$MCMCtheta[1]
  
  
  #loglikelihood is copied directly from the ergm.ego fit to the netest fit object.
  out$fit$loglikelihood <- x$loglikelihood
  
  ##Remove the first element (NA) of gradient associated with network size adjust.
  out$fit$gradient<-x$gradient[-1]
  
  ##Remove row 1 and col 1 for network size adjust from hessian matrix.
  ncol <- length(x$hessian[1,])
  out$fit$hessian <- x$hessian[2:ncol,2:ncol]
  
  ##Remove row 1 and col 1 for network size adjust from covar matrix.
  ncol <- length(x$covar[1,])
  out$fit$covar <- x$covar[2:ncol,2:ncol]
  
  #failure is copied directly from the ergm.ego fit to the netest fit object.
  out$fit$failure <- x$failure
  
  #network, newnetworks and newnetwork are copied directly from the ergm.ego fit to the netest fit object.
  out$fit$network <- x$network
  out$fit$newnetworks <- x$newnetworks
  out$fit$newnetwork <- x$newnetwork
  
  #coef.init is adjusted to remove the netsize.adj and adjust the edges term
  out$fit$coef.init<-x$coef.init[2:length(x$coef.init)]
  out$fit$coef.init[1]<-out$fit$coef.init[1] + x$coef.init[1]
  
  #est.cov: Remove row 1 and col 1 for network size adjust from est.cov
  ncol <- length(x$est.cov[1,])
  out$fit$est.cov <- x$est.cov[2:ncol,2:ncol]
  
  #coef.hist: combine edges and nesize adjustment, drop netsize adjustment.
  ncol <- length(x$coef.hist[1,])
  x$coef.hist[,3] <- x$coef.hist[,2] + x$coef.hist[,3]
  out$fit$coef.hist <- x$coef.hist[,-1]
  
  #stats.hist: drop netsize adjustment.
  out$fit$stats.hist <- x$stats.hist[,-1]
  
  #steplen.hist: Direct carryover
  out$fit$steplen.hist <- x$steplen.hist
  
  #control:  Merge edges and netsize.adj in $init
  out$fit$control <- x$control
  out$fit$control$init <- x$control$init[-1]
  out$fit$control$init[1] <- x$control$init[1] + x$control$init[2]
  
  #etamap: drop first element of each attribute
  out$fit$etamap <- x$etamap
  out$fit$etamap$cononical <- x$etamap$cononical[-(length(x$etamap$cononical))]
  out$fit$etamap$offsetmap <- x$etamap$offsetmap[-1]
  out$fit$etamap$offset <- x$etamap$offset[-1]
  out$fit$etamap$offsettheta <- x$etamap$offsettheta-1
  out$fit$etamap$curved <- x$etamap$curved
  out$fit$etamap$etalength <- x$etamap$etalength - 1
  
  #formula
  out$fit$formula <- x$formula
  z<-(as.formula("nw ~ o"))
  out$fit$formula[2] <- z[2]
  
  #JKB 4/24/18 - nw object needs to be added to the environment of the formation formula
  #patterned after netest() lines 182-183
  #DTH 4/46/2018 - testing cause of Males-Male ties
  # deleting all ties to creat an empty network
  
  #JKB 6/5/18: allow for an alternate starting network
  if (is.null(start.net)) {
    nw <- x$network
  } else nw <- start.net
  count <- network.edgecount(nw)
  delete.edges(nw,1:count)
  if (set.formula.nw) environment(out$fit$formula) <- environment()
  
  
  #target stats: drop the first element for netsizw.adj
  out$fit$target.stats <- x$target.stats[-1]
  
  #target.esteq: Direct copy
  out$fit$target.esteq <- x$target.esteq
  
  #constraints: Direct copy
  out$fit$constraints <- x$constraints
  
  #reference: Direct copy
  out$fit$reference <- x$reference
  
  #estimate: Direct copy
  out$fit$estimate <- x$estimate
  
  #offset: Drop the first element
  out$fit$offset <- x$offset[-1]
  
  #drop: Drop the first element
  out$fit$drop <- x$drop[-1]
  
  #estimable: Drop the first element
  out$fit$estimable <- x$estimable[-1]
  
  #null.lik: Direct copy
  out$fit$null.lik <- x$null.lik[-1]
  
  
  ## the netest object contains
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  #The objects other than fit.
  
  formation <- x$formula[-2]
  target.stats <- as.numeric(x$target.stats[2:length(x$target.stats)])
  # JKB 4/25/18 - changed next 2 lines, which turned target.stats into integer(0)
  # I'm presuming the goal was to remove the offset term(s)
  # DTH This works for me
  
  target.stats <- target.stats[!out$fit$offset]
  target.stats.names <- names(out$fit$target.stats)[!out$fit$offset]
  
  coef.form.crude <- x$coef[2:length(x$coef)]
  coef.form <- coef.form.crude
  coef.form[1]<-coef.form[1]+x$coef[1]
  coef.form.crude <- coef.form
  
  if (y$coef.adj[1] != -Inf){coef.form[1]<- coef.form[1] - y$coef.adj[1]}
  
  
  constraints <- x$constraints
  edapprox <- TRUE
  
  
  ##assighn the elements of the netest object
  out$formation <- formation
  out$target.stats <- target.stats
  
  ##This is redundant with the above section on target stats.
  
  #if (length(names(out$fit$coef)) == length(target.stats)) {
  #    out$target.stats.names <- names(out$fit$coef)
  #}
  #else {
  #    out$target.stats.names <- names(out$fit$coef)[!out$fit$offset]}
  
  out$target.stats.names <- target.stats.names
  
  
  out$coef.form <- coef.form
  out$dissolution <- y$dissolution
  out$coef.diss <- y
  out$constraints <- constraints
  out$edapprox <- edapprox
  
  
  
  class(out) <- "netest"
  return(out)
}


ee.netest.emily <- function (x,y, set.formula.nw=TRUE, start.net=NULL) {
  
  ## the netest object contains
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  
  ##Initialize output list
  out <- list()
  
  ## fit is an ergm object that contains.
  ##[1] "coef"          "sample"        "iterations"    "MCMCtheta"     "loglikelihood" "gradient"      "hessian"       "covar"
  ##[9] "failure"       "network"       "newnetworks"   "newnetwork"    "coef.init"     "est.cov"       "coef.hist"     "stats.hist"
  ##[17] "steplen.hist"  "control"       "etamap"        "formula"       "target.stats"  "target.esteq"  "constraints"   "reference"
  ##[25] "estimate"      "offset"        "drop"          "estimable"     "null.lik"
  
  ##Delete egodata
  x$egodata<-NULL
  
  ##Assign the ergm.ego object to fit in netest
  out$fit <- x
  class(out$fit) <- "ergm"
  
  ##Adjust the ergm coeficient to get rid of the network size adjustment
  out$fit$coef<-x$coef[2:length(x$coef)]
  out$fit$coef[1]<-out$fit$coef[1] + x$coef[1]
  
  #Remove the netsize.adj colum from "sample"
  ncol<-length(x$sample[1,])
  out$fit$sample<-x$sample[,-1]
  
  #Iterations is copied directly from the ergm.ego fit to the netest fit object.
  out$fit$iterations <- x$iterations
  
  ##Adjust the MCMCtheta coeficient to get rid of the network size adjustment
  out$fit$MCMCtheta<-x$MCMCtheta[2:length(x$MCMCtheta)]
  out$fit$MCMCtheta[1]<-out$fit$MCMCtheta[1] + x$MCMCtheta[1]
  
  
  #loglikelihood is copied directly from the ergm.ego fit to the netest fit object.
  out$fit$loglikelihood <- x$loglikelihood
  
  ##Remove the first element (NA) of gradient associated with network size adjust.
  out$fit$gradient<-x$gradient[-1]
  
  ##Remove row 1 and col 1 for network size adjust from hessian matrix.
  ncol <- length(x$hessian[1,])
  out$fit$hessian <- x$hessian[2:ncol,2:ncol]
  
  ##Remove row 1 and col 1 for network size adjust from covar matrix.
  ncol <- length(x$covar[1,])
  out$fit$covar <- x$covar[2:ncol,2:ncol]
  
  #failure is copied directly from the ergm.ego fit to the netest fit object.
  out$fit$failure <- x$failure
  
  #network, newnetworks and newnetwork are copied directly from the ergm.ego fit to the netest fit object.
  out$fit$network <- x$network
  out$fit$newnetworks <- x$newnetworks
  out$fit$newnetwork <- x$newnetwork
  
  #coef.init is adjusted to remove the netsize.adj and adjust the edges term
  out$fit$coef.init<-x$coef.init[2:length(x$coef.init)]
  out$fit$coef.init[1]<-out$fit$coef.init[1] + x$coef.init[1]
  
  #est.cov: Remove row 1 and col 1 for network size adjust from est.cov
  ncol <- length(x$est.cov[1,])
  out$fit$est.cov <- x$est.cov[2:ncol,2:ncol]
  
  #coef.hist: combine edges and nesize adjustment, drop netsize adjustment. 
  # THIS CHANGES THE COEF FOR THE THIRD ITEM, IT DOESN'T CHANGE THE EDGES. CHANGING to edit the edges term. 
  ncol <- length(x$coef.hist[1,])
  #x$coef.hist[,3] <- x$coef.hist[,2] + x$coef.hist[,3]
  x$coef.hist[,2] <- x$coef.hist[,1] + x$coef.hist[,2]
  out$fit$coef.hist <- x$coef.hist[,-1]
  
  #stats.hist: drop netsize adjustment.
  out$fit$stats.hist <- x$stats.hist[,-1]
  
  #steplen.hist: Direct carryover
  out$fit$steplen.hist <- x$steplen.hist
  
  #control:  Merge edges and netsize.adj in $init
  out$fit$control <- x$control
  out$fit$control$init <- x$control$init[-1]
  out$fit$control$init[1] <- x$control$init[1] + x$control$init[2]
  
  #etamap: drop first element of each attribute
  out$fit$etamap <- x$etamap
  out$fit$etamap$cononical <- x$etamap$cononical[-(length(x$etamap$cononical))]
  out$fit$etamap$offsetmap <- x$etamap$offsetmap[-1]
  out$fit$etamap$offset <- x$etamap$offset[-1]
  out$fit$etamap$offsettheta <- x$etamap$offsettheta-1
  out$fit$etamap$curved <- x$etamap$curved
  out$fit$etamap$etalength <- x$etamap$etalength - 1
  
  #formula
  out$fit$formula <- x$formula
  z<-(as.formula("nw ~ o"))
  out$fit$formula[2] <- z[2]
  
  #JKB 4/24/18 - nw object needs to be added to the environment of the formation formula
  #patterned after netest() lines 182-183
  #DTH 4/46/2018 - testing cause of Males-Male ties
  # deleting all ties to creat an empty network
  
  #JKB 6/5/18: allow for an alternate starting network
  if (is.null(start.net)) {
    nw <- x$network
  } else nw <- start.net
  count <- network.edgecount(nw)
  delete.edges(nw,1:count)
  if (set.formula.nw) environment(out$fit$formula) <- environment()
  
  
  #target stats: drop the first element for netsizw.adj
  out$fit$target.stats <- x$target.stats[-1]
  
  #target.esteq: Direct copy
  out$fit$target.esteq <- x$target.esteq
  
  #constraints: Direct copy
  out$fit$constraints <- x$constraints
  
  #reference: Direct copy
  out$fit$reference <- x$reference
  
  #estimate: Direct copy
  out$fit$estimate <- x$estimate
  
  #offset: Drop the first element
  out$fit$offset <- x$offset[-1]
  
  #drop: Drop the first element
  out$fit$drop <- x$drop[-1]
  
  #estimable: Drop the first element
  out$fit$estimable <- x$estimable[-1]
  
  #null.lik: Direct copy
  # out$fit$null.lik <- x$null.lik[-1]
  # the above doesn't copy the null.lik - do we want this to be `numeric(0)`?
  out$fit$null.lik <- x$null.lik[1]
  
  ## the netest object contains
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  #The objects other than fit.
  
  formation <- x$formula[-2]
  target.stats <- as.numeric(x$target.stats[2:length(x$target.stats)])
  # JKB 4/25/18 - changed next 2 lines, which turned target.stats into integer(0)
  # I'm presuming the goal was to remove the offset term(s)
  # DTH This works for me
  
  target.stats <- target.stats[!out$fit$offset]
  target.stats.names <- names(out$fit$target.stats)[!out$fit$offset]
  
  coef.form.crude <- x$coef[2:length(x$coef)]
  coef.form <- coef.form.crude
  coef.form[1]<-coef.form[1]+x$coef[1]
  coef.form.crude <- coef.form
  
  if (y$coef.adj[1] != -Inf){coef.form[1]<- coef.form[1] - y$coef.adj[1]}
  
  
  constraints <- x$constraints
  edapprox <- TRUE
  
  
  ##assighn the elements of the netest object
  out$formation <- formation
  out$target.stats <- target.stats
  
  ##This is redundant with the above section on target stats.
  
  #if (length(names(out$fit$coef)) == length(target.stats)) {
  #    out$target.stats.names <- names(out$fit$coef)
  #}
  #else {
  #    out$target.stats.names <- names(out$fit$coef)[!out$fit$offset]}
  
  out$target.stats.names <- target.stats.names
  
  
  out$coef.form <- coef.form
  out$dissolution <- y$dissolution
  out$coef.diss <- y
  out$constraints <- constraints
  out$edapprox <- edapprox
  
  
  
  class(out) <- "netest"
  return(out)
}
