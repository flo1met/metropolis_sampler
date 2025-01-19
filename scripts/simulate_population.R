### create population data ###
set.seed(1337)
N <- 1e6 # sample size
X1 <- rnorm(N, 1, 0.5) # create continuous covariate
X2 <- rbinom(N, 1, 0.5) # create binary covariate

Z <- 0.5 + 0.8*X1 + 0.3*X2 # linear combination
Pr <- 1/(1+exp(-Z)) # inv-logit function

Y <- rbinom(N, 1, Pr) # create outcome variabele

# see if it worked:
glm(Y~X1+X2, family = "binomial") # -> looks good

pop <- data.frame(X1 = X1, X2 = X2, Y = Y)

save(pop, file = "data/population.RData")
