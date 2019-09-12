##
## Simulação de distribuições a priori
##
nsim <- 10000

## Exemplo 1:
## Y \sim B(\theta)
## \theta \sim B(1, 1)

th <- rbeta(nsim, 1, 1)
hist(th, prob=TRUE); lines(density(th))
curve(dbeta(x, 1, 1), from=0, to=1, col=4, add=TRUE)

th2 <- th^2
hist(th2, prob=TRUE); lines(density(th2))
dth2 <- function(th) 1/(2*sqrt(th))
curve(dth2, from=0, to=1, col=4, add=TRUE)

th.logit <- log(th/(1-th))
hist(th.logit, prob=TRUE); lines(density(th.logit))

##
th <- rbeta(nsim, 5, 5)
hist(th, prob=TRUE); lines(density(th))
curve(dbeta(x, 5, 5), from=0, to=1, col=4, add=TRUE)

th.logit <- log(th/(1-th))
hist(th.logit, prob=TRUE); lines(density(th.logit))

##
th <- rbeta(nsim, 6, 2)
hist(th, prob=TRUE); lines(density(th))
curve(dbeta(x, 6, 2), from=0, to=1, col=4, add=TRUE)

th.logit <- log(th/(1-th))
hist(th.logit, prob=TRUE); lines(density(th.logit))

##
th <- rbeta(nsim, 1/2, 1/2)
hist(th, prob=TRUE); lines(density(th))
curve(dbeta(x, 1/2, 1/2), from=0, to=1, col=4, add=TRUE)

th2 <- th^2
hist(th2, prob=TRUE); lines(density(th2))
dth2 <- function(th) 1/(2*sqrt(th))
curve(dth2, from=0, to=1, col=4, add=TRUE)

## experimentar outras possibilidades!!!


## Exemplo 2:
## Y \sim B(\theta)
## \phi = logit(\theta) = log(\theta/(1-\theta)) \sim N(mu, \sigma^2)
th.logit <- rnorm(nsim, m=0, sd=3)
hist(th.logit, prob=TRUE); lines(density(th.logit))
curve(dnorm(x, 0, 3), from=-10, to=10, col=4, add=TRUE)

th <- exp(th.logit)/(1+exp(th.logit))
hist(th, prob=TRUE); lines(density(th))
MASS:::fitdistr(th, densfun="beta", start=list(shape1=1, shape2=1))

## experimentar outras possibilidades!!!


## Exemplo 3:
## Y \sim P(\theta)
## \theta \sim G(shape=\alpha, scale=\beta)
## f(\theta) = \frac{\beta^\alpha}{\Gamma(\alpha)} 

th <- rgamma(nsim, 2, 2)
hist(th, prob=TRUE); lines(density(th))
curve(dgamma(x, 2, 2), from=0, to=8, col=4, add=TRUE)

## \phi = P(Y > 0) = 1 - P[Y=0] = 1 - exp{-\theta}
phi <- 1 - exp(-th)
hist(phi, prob=TRUE); lines(density(phi))

## \psi = log(theta)
psi <- log(th)
hist(psi, prob=TRUE); lines(density(psi))

## Como simular da priori de Jeffreys neste caso?


## experimentar outras possibilidades!!!


## Exemplo 4:
## Y \sim P(\theta)
## psi = log(\theta) \sim N(\mu, \sigma^2)

psi <- rnorm(nsim, m=0, sd=0.8)
hist(psi, prob=TRUE); lines(density(psi))

th <- exp(psi)
hist(th, prob=TRUE); lines(density(th))
curve(dlnorm(x, 0, 0.8), from=0, to=8, col=4, add=TRUE)

## \phi = P(Y > 0) = 1 - P[Y=0] = 1 - exp{-\theta}
phi <- 1 - exp(-th)
hist(phi, prob=TRUE); lines(density(phi))




