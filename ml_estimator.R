# Define the log-likelihood function
log_likelihood <- function(beta, x, y) {
  p <- 1 / (1 + exp(-beta[1] - beta[2]*x))
  sum(y * log(p) + (1 - y) * log(1 - p))
}

# Define the gradient of the log-likelihood function
gradient <- function(beta, x, y) {
  p <- 1 / (1 + exp(-beta[1] - beta[2]*x))
  db0 <- sum(y - p)
  db1 <- sum((y - p) * x)
  c(db0, db1)
}

# Steepest ascent function with alpha0 when new iteration
steepest_ascent_constant <- function(beta_start, x, y, tol=1e-5, alpha0=1) {
  beta <- beta_start
  log_likelihood_count <- 0
  gradient_count <- 0
  conv <- 999
  while (conv > tol) {
    alpha <- alpha0
    beta_old <- beta
    beta <- beta_old + alpha * gradient(beta_old, x, y)
    gradient_count <- gradient_count + 1
    log_likelihood_count <- log_likelihood_count + 2
    while (log_likelihood(beta, x , y) < log_likelihood(beta_old, x , y)) {
      alpha <- alpha/2
      beta <- beta_old + alpha * gradient(beta_old, x , y)
      gradient_count <- gradient_count + 1
      log_likelihood_count <- log_likelihood_count + 2
    }
    conv <- sum((beta - beta_old) * (beta - beta_old))
  }
  result <- list(coefficients=c(beta0=beta[1],beta1=beta[2]), counts=c(func_count=log_likelihood_count, gradient_count=gradient_count))
  return(result)
}

# Steepest ascent function with decreasing alpha when new iteration
steepest_ascent_decrease <- function(beta_start, x, y, tol=1e-5, alpha0=1) {
  beta <- beta_start
  log_likelihood_count <- 0
  gradient_count <- 0
  conv <- 999
  alpha <- alpha0
  while (conv > tol) {
    beta_old <- beta
    beta <- beta_old + alpha * gradient(beta_old, x, y)
    gradient_count <- gradient_count + 1
    log_likelihood_count <- log_likelihood_count + 2
    while (log_likelihood(beta, x , y) < log_likelihood(beta_old, x , y)) {
      alpha <- alpha/2
      beta <- beta_old + alpha * gradient(beta_old, x , y)
      gradient_count <- gradient_count + 1
      log_likelihood_count <- log_likelihood_count + 2
    }
    conv <- sum((beta - beta_old) * (beta - beta_old))
  }
  result <- list(coefficients=c(beta0=beta[1],beta1=beta[2]), counts=c(func_count=log_likelihood_count, gradient_count=gradient_count))
  return(result)
}

# drug and placebo data
x <- c(0, 0, 0, 0.1, 0.1, 0.3, 0.3, 0.9, 0.9, 0.9)
y <- c(0, 0, 1, 0, 1, 1, 1, 0, 1, 1)

# Initial values for beta
beta_start <- c(-0.2, 1)

# Compute the ML estimator
steepest_ascent_constant_result <- steepest_ascent_constant(beta_start, x, y)

steepest_ascent_decrease_result <- steepest_ascent_decrease(beta_start, x, y)

# Print the results
cat("=====Using constant alpha when new iteration for steepest ascent algorithm======\n")
print(steepest_ascent_constant_result)
cat("=====Using decrease alpha when new iteration for steepest ascent algorithm======\n")
print(steepest_ascent_decrease_result)


# BFGS
result_bfgs <- optim(beta_start, log_likelihood, gradient, x=x, y=y, method="BFGS", control = list(fnscale = -1))
print(result_bfgs)

# Nelder-Mead
result_nelder <- optim(beta_start, log_likelihood, gradient, x=x, y=y, method="Nelder-Mead", control = list(fnscale = -1))
print(result_nelder)


# Fit the model
fit <- glm(y ~ x, family=binomial(link="logit"))

# Print the results
print(summary(fit))

