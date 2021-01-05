#example 3.3
library(MCMCglmm)

f <- function(x) {
 (cos(50*x) + sin(20*x))^2 
}

x <- seq(0,1,0.01)
fx <- f(x)
plot(x,fx,'l')

#monte carlo using the uniform distribution
n <- 10000
x <- runif(n)
hist(x)
mean(f(x)/dunif(x))


#monte carlo using the normal distribution
n <- 20000000
x <- rtnorm(n,0.5,0.15,0,1)
hist(x)
mean(f(x)/dnorm(x,0.5,0.15))

