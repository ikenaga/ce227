#-----------------------------------------------------------------------
# Introdução à inferência Bayesiana · CE 227
# web.leg.ufpr.br/ensino/EMR
#
#                                            Prof. Dr. Fernando Mayer
#                                       Prof. Dr. Paulo Justiniano R. Jr
#
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2019-Ago-29 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#



##
## Simulação de distribuições a priori
##
par.ori <- par(no.readonly=TRUE)
nsim <- 10000

## Exemplo 1:
## Y \sim B(\theta)
## \theta \sim B(1, 1)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,1,0))
th <- rbeta(nsim, 1, 1)
hist(th, prob=TRUE); lines(density(th))
curve(dbeta(x, 1, 1), from=0, to=1, col=4, add=TRUE)

th2 <- th^2
hist(th2, prob=TRUE); lines(density(th2))
dth2 <- function(th) 1/(2*sqrt(th))
curve(dth2, from=0, to=1, col=4, add=TRUE)

th.odds <- th/(1-th)
hist(th.odds, prob=TRUE); lines(density(th.odds))

th.logit <- log(th/(1-th))
hist(th.logit, prob=TRUE); lines(density(th.logit))
sd(th.logit)

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
##

th.logit <- rgamma(nsim, sh=1.5, sc=0.6)
hist(th.logit, prob=TRUE); lines(density(th.logit))
curve(dgamma(x, sh=1.5, sc=0.6), from=0, to=5, col=4, add=TRUE)

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
## Uma opção é adotar uma priori discreta (aproximar a contínua por discreta)
th <- seq(0.1, 5, by=0.1)
pth <- 1/sqrt(th)
pth <- pth/sum(pth)
sum(pth)
plot(th, pth, type="h")

sample(th, 10, prob=pth)


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

