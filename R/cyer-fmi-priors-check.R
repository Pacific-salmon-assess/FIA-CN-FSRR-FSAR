logit_marine <- seq(-3.6, 1, 0.1)  # Reasonable range for logit(ER)

set.seed(123)
for(i in 1:100) {
  b0 <- rnorm(1, 0, 0.25)
  b1 <- rnorm(1, 1, 0.25)
  logit_can <- b0 + b1 * logit_marine
  
  if(i == 1) {
    plot(plogis(logit_marine), plogis(logit_can), type='l', 
         xlim=c(0,1), ylim=c(0,1), col=rgb(0,0,0,0.1),
         xlab="Marine FMI ER", ylab="Canadian ER")
    abline(0, 1, col='red', lwd=2)
  } else {
    lines(plogis(logit_marine), plogis(logit_can), col=rgb(0,0,0,0.1))
  }
}



marine_fmi <- seq(0.001, 0.999, 0.01)  # Proportions (avoiding 0 and 1)

set.seed(123)
plot(NULL, xlim=c(0,1), ylim=c(0,1), 
     xlab="Marine FMI ER", ylab="Canadian ER",
     main="Prior predictive: untransformed predictor")
abline(0, 1, col='red', lwd=2)

for(i in 1:100) {
  b0 <- rnorm(1, 0, 3)      # Your intercept prior
  b1 <- rnorm(1, 0, 5)      # Your slope prior
  logit_can <- b0 + b1 * marine_fmi
  can_er <- plogis(logit_can)
  
  lines(marine_fmi, can_er, col=rgb(0,0,0,0.1))
}
