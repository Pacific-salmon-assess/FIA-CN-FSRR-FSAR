logit_marine <- qlogis(seq(0.001, 0.4, 0.01))  # Reasonable range for logit(ER)
logit_marine_scaled <- logit_marine - mean(logit_marine)

set.seed(123)
for(i in 1:100) {
  
  b0 <- rnorm(1, mean(logit_marine), 1)
  b1 <- rnorm(1, 1, 1)
  logit_can <- b0 + b1 * logit_marine_scaled
  
  if(i == 1) {
    plot(plogis(logit_marine), plogis(logit_can), type='l', 
         xlim=c(0,0.4), ylim=c(0,0.4), col=rgb(0,0,0,0.2),
         xlab="Marine FMI ER", ylab="Canadian ER")
    abline(0, 1, col='red', lwd=2)
  } else {
    lines(plogis(logit_marine), plogis(logit_can), col=rgb(0,0,0,0.2))
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


library(ggplot2)

# Simulate what different priors imply for stock-level variation
set.seed(123)
n_stocks <- 6
n_draws <- 1000

# Different exponential rates
lambda_values <- c(0.5, 1, 2, 5)

results <- data.frame()

for(lambda in lambda_values) {
  for(i in 1:n_draws) {
    sd_stock <- rexp(1, rate = lambda)
    stock_effects <- rnorm(n_stocks, 0, sd_stock)
    
    # At mean logit_fmi (let's say -0.9)
    mean_intercept <- -0.9
    stock_intercepts <- mean_intercept + stock_effects
    
    # Convert to probability scale
    stock_probs <- plogis(stock_intercepts)
    
    results <- rbind(results, data.frame(
      lambda = paste0("exp(", lambda, ")"),
      draw = i,
      stock = 1:n_stocks,
      logit_intercept = stock_intercepts,
      prob = stock_probs
    ))
  }
}

# Visualize on probability scale
ggplot(results, aes(x = prob, fill = lambda)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~lambda, ncol = 1) +
  labs(title = "Prior predictive: Stock-level variation in exploitation rate",
       x = "Exploitation rate (at mean marine_fmi)",
       y = "Density") +
  theme_minimal()

