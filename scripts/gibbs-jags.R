##======================================================================
## Exemplo 3.4
## Y ~ N(mu, sigma^2)   com tau = sigma^{-2}
## mu ~ N(0, s^2)       com tau_0 = s^{-2}
## tau ~ Gama(a,b)

## Condicionais são
## mu | tau ~ N(m, C)
## tau | mu ~ Gama(a + n/2, b + 1/2 sum((y_i - mu)^2))
## Onde:
## m = C\bar{y}
## C = 1/(n\tau + tau_0)

## Simula valores
set.seed(2)
y <- rnorm(30, 0, 1)
mean(y)
var(y)

library(runjags)

## Dados
datalist <- dump.format(list(y = y))
params <- c("mu", "tau")
inicial <- dump.format(list(mu = 0, tau = 1))

## Modelo
mod <- "model{
mu ~ dnorm(0, .001)
tau ~ dgamma(2, 2)
for(i in 1:length(y)){
y[i] ~ dnorm(mu, tau)
}
}"

## Ajuste
m.jags <- run.jags(model = mod, monitor = params,
                   data = datalist, inits = inicial,
                   n.chains = 1, burnin = 50000, thin = 5,
                   sample = 10000)

## Resultados
m.jags
names(m.jags)
head(m.jags$mcmc)
str(m.jags$mcmc)
mu.jags <- m.jags$mcmc[[1]][, 1]
tau.jags <- m.jags$mcmc[[1]][, 2]
str(mu.jags)
str(tau.jags)
summary(mu.jags)
summary(tau.jags)
par(mfrow = c(1,2))
plot(density(mu.jags))
plot(density(tau.jags))

## Compara amostra com distribuicoes teoricas
par(mfrow = c(2, 2))
hist(x[,1], freq = FALSE, main = "", xlab = expression(mu))
curve(dnorm(x, ybar/(n * 1 + tau0), sqrt(1/(n * 1 + tau0))),
      add = TRUE, col = 2)
hist(x[,2], freq = FALSE, main = "", xlab = expression(sigma^2))
curve(dgamma(x, a + n/2, b + sum((y - 0)^2)/2),
      add = TRUE, col = 2)

hist(mu.jags, freq = FALSE, main = "", xlab = expression(mu))
curve(dnorm(x, ybar/(n * 1 + tau0), sqrt(1/(n * 1 + tau0))),
      add = TRUE, col = 2)
hist(tau.jags, freq = FALSE, main = "", xlab = expression(sigma^2))
curve(dgamma(x, a + n/2, b + sum((y - 0)^2)/2),
      add = TRUE, col = 2)

par(mfrow = c(1, 1))


##======================================================================
## Exemplo 3.1 - Stuart Coles
## Y ~ N(mu, sigma^2)   com tau = sigma^{-2}
## mu ~ N(mu_0, s^2)    com tau_0 = s^{-2}
## tau ~ Gama(a,b)

## Condicionais são
## mu | sigma^2 ~ N(m, C)
## tau | mu ~ Gama(a + n/2, b + 1/2 sum((y_i - mu)^2))
## Onde:
## m = C \tau sum(x_i) + \mu_0 \tau_0
## C = 1/(n\tau + \tau_0)

## Simula valores
set.seed(2)
y <- rnorm(1000, 0, 5)

## Dados
datalist <- dump.format(list(y = y))
params <- c("mu", "tau")
inicial <- dump.format(list(mu = 0, tau = 1))

## Modelo
mod <- "model{
for(i in 1:length(y)){
y[i] ~ dnorm(mu, tau)
}
mu ~ dnorm(0, .0001)
tau ~ dgamma(2, 2)
}"

## Ajuste
m.jags <- run.jags(model = mod, monitor = params,
                   data = datalist, inits = inicial,
                   n.chains = 1, burnin = 50000, thin = 5,
                   sample = 10000)

## Resultados
m.jags
names(m.jags)
head(m.jags$mcmc)
str(m.jags$mcmc)
mu.jags <- m.jags$mcmc[[1]][, 1]
tau.jags <- m.jags$mcmc[[1]][, 2]
str(mu.jags)
str(tau.jags)
summary(mu.jags)
summary(tau.jags)
par(mfrow = c(1,2))
plot(density(mu.jags))
plot(density(tau.jags))

## Compara amostra com distribuicoes teoricas
par(mfrow = c(2, 2))
hist(x[,1], freq = FALSE, main = "", xlab = expression(mu))
C <- 1/(n * tau + tau0)
m <- C * (tau * ysum + mu0 * tau0)
curve(dnorm(x, m, sqrt(C)), add = TRUE, col = 2)
hist(x[,2], freq = FALSE, main = "", xlab = expression(sigma^2))
curve(dgamma(x, a + n/2, b + sum((y - 0)^2)/2),
      add = TRUE, col = 2)

hist(mu.jags, freq = FALSE, main = "", xlab = expression(mu))
C <- 1/(n * tau + tau0)
m <- C * (tau * ysum + mu0 * tau0)
curve(dnorm(x, m, sqrt(C)), add = TRUE, col = 2)
hist(x[,2], freq = FALSE, main = "", xlab = expression(sigma^2))
curve(dgamma(x, a + n/2, b + sum((y - 0)^2)/2),
      add = TRUE, col = 2)
par(mfrow = c(1, 1))
