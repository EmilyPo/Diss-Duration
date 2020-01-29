new.netest <- function (x, y, set.formula.nw = TRUE, start.net = NULL, O.T.adjust = FALSE) {
  out <- list()
  x$egodata <- NULL
  out$fit <- x
  class(out$fit) <- "ergm"
  out$fit$coef <- x$coef[2:length(x$coef)]
  out$fit$coef[1] <- out$fit$coef[1] + x$coef[1]
  ncol <- length(x$sample[1, ])
  out$fit$sample <- x$sample[, -1]
  out$fit$iterations <- x$iterations
  out$fit$MCMCtheta <- x$MCMCtheta[2:length(x$MCMCtheta)]
  out$fit$MCMCtheta[1] <- out$fit$MCMCtheta[1] + x$MCMCtheta[1]
  out$fit$loglikelihood <- x$loglikelihood
  out$fit$gradient <- x$gradient[-1]
  ncol <- length(x$hessian[1, ])
  out$fit$hessian <- x$hessian[2:ncol, 2:ncol]
  ncol <- length(x$covar[1, ])
  out$fit$covar <- x$covar[2:ncol, 2:ncol]
  out$fit$failure <- x$failure
  out$fit$network <- x$network
  out$fit$newnetworks <- x$newnetworks
  out$fit$newnetwork <- x$newnetwork
  out$fit$coef.init <- x$coef.init[2:length(x$coef.init)]
  out$fit$coef.init[1] <- out$fit$coef.init[1] + x$coef.init[1]
  ncol <- length(x$est.cov[1, ])
  out$fit$est.cov <- x$est.cov[2:ncol, 2:ncol]
  ncol <- length(x$coef.hist[1, ])
  x$coef.hist[, 2] <- x$coef.hist[, 1] + x$coef.hist[, 2]
  out$fit$coef.hist <- x$coef.hist[, -1]
  out$fit$stats.hist <- x$stats.hist[, -1]
  out$fit$steplen.hist <- x$steplen.hist
  out$fit$control <- x$control
  out$fit$control$init <- x$control$init[-1]
  out$fit$control$init[1] <- x$control$init[1] + x$control$init[2]
  out$fit$etamap <- x$etamap
  out$fit$etamap$cononical <- x$etamap$cononical[-(length(x$etamap$cononical))]
  out$fit$etamap$offsetmap <- x$etamap$offsetmap[-1]
  out$fit$etamap$offset <- x$etamap$offset[-1]
  out$fit$etamap$offsettheta <- x$etamap$offsettheta - 1
  out$fit$etamap$curved <- x$etamap$curved
  out$fit$etamap$etalength <- x$etamap$etalength - 1
  out$fit$formula <- x$formula
  z <- (as.formula("nw ~ o"))
  out$fit$formula[2] <- z[2]
  if (is.null(start.net)) {
    nw <- x$network
  }
  else nw <- start.net
  count <- network.edgecount(nw)
  delete.edges(nw, 1:count)
  if (set.formula.nw) 
    environment(out$fit$formula) <- environment()
  if (is.na(x$target.stats[1]) == TRUE) {
    out$fit$target.stats <- x$target.stats[-1]
  }
  else out$fit$target.stats <- x$target.stats
  out$fit$target.esteq <- x$target.esteq
  out$fit$constraints <- x$constraints
  out$fit$reference <- x$reference
  out$fit$estimate <- x$estimate
  out$fit$offset <- x$offset[-1]
  out$fit$drop <- x$drop[-1]
  out$fit$estimable <- x$estimable[-1]
  out$fit$null.lik <- x$null.lik[-1]
  formation <- x$formula[-2]
  target.stats <- as.numeric(out$fit$target.stats)
  target.stats <- target.stats[!out$fit$offset]
  target.stats.names <- names(out$fit$target.stats)[!out$fit$offset]
  coef.form.crude <- x$coef[2:length(x$coef)]
  coef.form <- coef.form.crude
  coef.form[1] <- coef.form[1] + x$coef[1]
  coef.form.crude <- coef.form
  if (y$coef.adj[1] != -Inf) {
    coef.form[1] <- coef.form[1] - y$coef.adj[1]
  }
  if (O.T.adjust == TRUE) {
    coef.form[1] <- coef.form[1] + -log(52)
    out$fit$coef[1] <- out$fit$coef[1] + -log(52)
  }
  constraints <- x$constraints
  edapprox <- TRUE
  out$formation <- formation
  out$target.stats <- target.stats
  out$target.stats.names <- target.stats.names
  out$coef.form <- coef.form
  out$dissolution <- y$dissolution
  out$coef.diss <- y
  out$constraints <- constraints
  out$edapprox <- edapprox
  class(out) <- "netest"
  return(out)
}
