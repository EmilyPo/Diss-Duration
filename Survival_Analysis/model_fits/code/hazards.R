# old hazard functions
# need to find way to plot hazard funcations of flexsurv objects without default plotting KM


# h(t)=p(t)/S(t)

# full plot
plot(NULL,
     xlab = "t, months", ylab = "H(t)",
     xlim=c(0, max(x$t_c)), ylim=c(0,0.15),
     main="Hazard of Dissolution"
)
#exponential
lines(x = durVec,
      y = dexp(durVec, coef(expFit))/pexp(durVec, coef(expFit), lower.tail = FALSE),
      col="blue")
#weibull
lines(x= durVec,
      y=dweibull(durVec, shape=coef(weibFit)[1], scale=coef(weibFit)[2])/pweibull(durVec, shape=coef(weibFit)[1], scale=coef(weibFit)[2], lower.tail = FALSE),
      col="red")
#gamma
lines(x= durVec,
      y=dgamma(durVec, shape=coef(gammaFit)[1], scale=coef(gammaFit)[2])/pgamma(durVec, shape=coef(gammaFit)[1], scale=coef(gammaFit)[2], lower.tail = FALSE),
      col="darkgreen")
#mixed exponential
lines(x = durVec,
      y = (invLogit(coef(mixExpFit)[3]) * dexp(durVec, exp(coef(mixExpFit)[1])) +
             (1-invLogit(coef(mixExpFit)[3])) * dexp(durVec, exp(coef(mixExpFit)[2]))) / 
        (invLogit(coef(mixExpFit)[3]) * pexp(durVec, exp(coef(mixExpFit)[1]), lower.tail = FALSE) +
           (1-invLogit(coef(mixExpFit)[3])) * pexp(durVec, exp(coef(mixExpFit)[2]), lower.tail = FALSE)),
      col ="purple"
)
legend(300,.15, legend = c("Kaplan-Meier", "Exponential", "Mixed Exponential", "Gamma", "Weibull"),
       col = c("black", "blue", "purple", "darkgreen", "red"), cex=0.7, lwd=1)

# zoomed 
plot(NULL,
     xlab = "t, months", ylab = "H(t)",
     xlim=c(0,150), ylim=c(0,0.15),
     main="(Zoomed) Hazard of Dissolution"
)
#exponential
lines(x = durVec,
      y = dexp(durVec, coef(expFit))/pexp(durVec, coef(expFit), lower.tail = FALSE),
      col="blue")
#weibull
lines(x= durVec,
      y=dweibull(durVec, shape=coef(weibFit)[1], scale=coef(weibFit)[2])/pweibull(durVec, shape=coef(weibFit)[1], scale=coef(weibFit)[2], lower.tail = FALSE),
      col="red")
#gamma
lines(x= durVec,
      y=dgamma(durVec, shape=coef(gammaFit)[1], scale=coef(gammaFit)[2])/pgamma(durVec, shape=coef(gammaFit)[1], scale=coef(gammaFit)[2], lower.tail = FALSE),
      col="darkgreen")
#mixed exponential
lines(x = durVec,
      y = (invLogit(coef(mixExpFit)[3]) * dexp(durVec, exp(coef(mixExpFit)[1])) +
             (1-invLogit(coef(mixExpFit)[3])) * dexp(durVec, exp(coef(mixExpFit)[2]))) / 
        (invLogit(coef(mixExpFit)[3]) * pexp(durVec, exp(coef(mixExpFit)[1]), lower.tail = FALSE) +
           (1-invLogit(coef(mixExpFit)[3])) * pexp(durVec, exp(coef(mixExpFit)[2]), lower.tail = FALSE)),
      col ="purple"
)
legend(75,.15, legend = c("Kaplan-Meier", "Exponential", "Mixed Exponential", "Gamma", "Weibull"),
       col = c("black", "blue", "purple", "darkgreen", "red"), cex=0.7, lwd=1)

# using the flexsurv objects
plot(km_weighted, 
     conf.int = FALSE,
     xlab = "t, months", ylab = "S(t)",
     main="Survival of Relationships"
)
rug(dat$t_c, lwd=0.07)

plot(expFit, col="blue", type = "hazard")
lines(weibFit, col="red", type = "hazard")
lines(gammaFit, col="darkgreen", type = "hazard")

lines(x = durVec,
      y = invLogit(coef(mixExpFit)[3]) *
        pexp(q = durVec, exp(coef(mixExpFit)[1]), lower.tail = FALSE) +
        (1-invLogit(coef(mixExpFit)[3])) *
        pexp(q = durVec, exp(coef(mixExpFit)[2]), lower.tail = FALSE),
      col ="purple", lwd=2
)
legend(300,1, legend = c("Kaplan-Meier", "Exponential", "Mixed Exponential", "Weibull", "Gamme"),
       col = c("black", "blue", "purple", "red", "darkgreen"), cex=0.7, lwd=2)
