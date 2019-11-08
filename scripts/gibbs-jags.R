##======================================================================
## Exemplo da baleia (ver slides)
x <- 10
n <- 15
teta <- seq(0, 2, length = 200)
alfa <- 1
beta <- 0.1
## Calcula a densidade da priori
priori.ni <- dgamma(teta, alfa, beta)
(alfa.star <- alfa + x)
(beta.star <- beta + n)

post.ni <- dgamma(teta, alfa.star, beta.star)
plot(teta, post.ni, type = "l", xlab = expression(theta),
     ylab = "Densidade de probabilidade")
lines(teta, priori.ni, lty = 2)
lines(teta, dgamma(teta, 1 + x, n), col = 2, lty = 2)
legend("topright",  lty = c(1, 2, 2), col = c(1, 1, 2),
       legend = c("Posterior", "Priori", "Verossimilhança"))

## Dados
datalist <- dump.format(list(x = x, n = n))
params <- c("theta")
inicial <- dump.format(list(theta = 0.5))

## Modelo
mod <- "model{
x ~ dpois(n * theta)
theta ~ dgamma(1, 0.1)
}"

## Ajuste
library(runjags)
m.jags <- run.jags(
    model = mod, monitor = params, data = datalist, inits = inicial,
    n.chains = 1, burnin = 5000, thin = 5, sample = 10000
)

## Resultados
m.jags
qgamma(c(0.025, 0.5, 0.975), alfa.star, beta.star)
plot(m.jags)

## Usando prioris informativas
alfa.i <- 4.5
beta.i <- 10
## Calcula a densidade da priori
priori.i <- dgamma(teta, alfa.i, beta.i)
alfa.star.i <- alfa.i + x
beta.star.i <- beta.i + n
## Cálculo da densidade da posterior com a priori informativa
post.i <- dgamma(teta, alfa.star.i, beta.star.i)
## Visualização
plot(teta, post.i, type = "l", xlab = expression(theta),
     ylab = "Densidade de probabilidade")
lines(teta, priori.i, lty = 2)
lines(teta, dgamma(teta, 1 + x, n), col = 2, lty = 2)
legend("topright",  lty = c(1, 2, 2), col = c(1, 1, 2),
       legend = c("Posterior", "Priori", "Verossimilhança"))

## Modelo
mod <- "model{
x ~ dpois(n * theta)
theta ~ dgamma(4.5, 10)
}"

## Ajuste
m.jags <- run.jags(
    model = mod, monitor = params, data = datalist, inits = inicial,
    n.chains = 1, burnin = 5000, thin = 5, sample = 10000
)

## Resultados
m.jags
qgamma(c(0.025, 0.5, 0.975), alfa.star.i, beta.star.i)
plot(m.jags)

##======================================================================
## Exemplo Poisson-gama (PJ)

## Dados
set.seed(2018)
ctes <- list(a=3, c=2.5, d=0.8, n=50)
betas <- with(ctes, 1/rgamma(n, shape=c, scale=d))
c(mean(betas),var(betas))
lambdas <- with(ctes, rgamma(n, shape=a, rate=betas))
(ctes$y <- rpois(ctes$n, lambda=lambdas))
with(ctes, c(media=mean(y), var=var(y)))
with(ctes, plot(prop.table(table(y)), type="h", ylim=c(0,0.3)))
with(ctes,lines((0:max(y))+0.1, dpois(0:max(y), lambda=mean(y)),
                type="h", col=2))

EVIG <- function(a,b){
    E <- ifelse(a > 1, 1/(b * (a-1)), "não pode ser calculada para este valor de a")
    V  <- ifelse(a > 2, 1/((b^2) * ((a-1)^2) * (a-2)), "não pode ser calculada para este valor de a")
    return(c(E,V))
}
set.seed(2018)
ctes <- list(a=3, c=2.5, d=0.8, n=50)
with(ctes, EVIG(c, d))
betas <- with(ctes, 1/rgamma(n, shape=c, scale=d))
c(mean(betas),var(betas))
lambdas <- with(ctes, rgamma(n, shape=a, rate=betas))
(ctes$y <- rpois(ctes$n, lambda=lambdas))
with(ctes, c(media=mean(y), var=var(y)))
with(ctes, plot(prop.table(table(y)), type="h", ylim=c(0,0.3)))
with(ctes,lines((0:max(y))+0.1, dpois(0:max(y), lambda=mean(y)), type="h", col=2))
##
ctes$sumY <- sum(ctes$y)
##
N <- 11000
B <- 1000
beta.sam <- lambda.sam <- numeric(N)
beta.sam[1] <- lambda.sam[1] <- 10
{
    for(i in 2:N){
        beta.sam[i] <- with(ctes, 1/rgamma(1, shape=a+c, scale=d/(d*lambda.sam[i-1]+1)))
        lambda.sam[i] <- with(ctes, rgamma(1, shape=ctes$a+sumY, scale=beta.sam[i]/(n*beta.sam[i]+1)))
    }
}

par(mfrow=c(2,1))
plot(beta.sam, type="l")
plot(lambda.sam, type="l")
## retirando amostras consideradas aquecimento
beta.sam <- beta.sam[-(1:B)]
lambda.sam <- lambda.sam[-(1:B)]
plot(beta.sam, type="l")
plot(lambda.sam, type="l")
plot(log(beta.sam), type="l")
plot(lambda.sam, type="l")


datalist <- dump.format(list(y = ctes$y))
params <- c("lambda", "beta")
inicial <- dump.format(list(lambda = 10, xi = 10))

## Modelo
mod <- "model{
for(i in 1:length(y)){
y[i] ~ dpois(lambda)
}
lambda ~ dgamma(3, xi)
xi ~ dgamma(2.5, 0.8)
beta <- 1/xi
}"

## Ajuste
m.jags <- run.jags(
    model = mod, monitor = params, data = datalist, inits = inicial,
    n.chains = 1, burnin = 5000, thin = 5, sample = 10000
)

## Resultados
m.jags
quantile(lambda.sam, probs = c(.025, .5, .975))
quantile(beta.sam, probs = c(.025, .5, .975))

par(mfrow=c(1,2))
plot(density(lambda.sam)); abline(v=mean(lambdas)); rug(lambdas)
plot(density(beta.sam)); abline(v=mean(betas)); rug(betas)
par(mfrow=c(1,1))

plot(m.jags)

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

##----------------------------------------------------------------------
## Implementação de Gibbs

## Define constantes
N <- 10000
burn <- 5000
thin <- 5
Nsim <- burn + (N * thin)
## Matriz para armazenar as amostras
X <- matrix(0, Nsim, 2)

## Simula valores
set.seed(2)
n <- 30
mu <- 0
sigma2 <- 25
tau <- 1/sigma2
y <- rnorm(n, mu, sqrt(sigma2))
mean(y)
var(y)
ysum <- sum(y)

## Prioris
## Valores da priori para mu ~ N(mu0, tau0)
mu0 <- 0
sigma20 <- 1
tau0 <- 1/sigma20
## Valores da priori para tau ~ Gama(a, b)
## Aqui da para testar varios valores para ver o impacto da priori
a <- .001
b <- .001
curve(dgamma(x, a, b), from = 0, to = 5)

## Gera a cadeia
X[1, ] <- c(mean(y), 1/var(y)) # Valores iniciais
for (i in 2:Nsim) {
    tau.iter <- X[i - 1, 2]
    C <- 1/(n * tau.iter + tau0)
    m <- C * (tau.iter * ysum + mu0 * tau0)
    X[i, 1] <- rnorm(1, m, sqrt(C))
    mu.iter <- X[i, 1]
    X[i, 2] <- rgamma(1, a + n/2, b + sum((y - mu.iter)^2)/2)
}

## Burnin
Xsim <- X[-(1:burn), ]
dim(Xsim)
## Thinning
Xsim <- Xsim[seq(1, nrow(Xsim), thin), ]
dim(Xsim)

## Cadeias (mude os valores iniciais para ver convergencia)
par(mfrow = c(2, 1))
plot.ts(Xsim[, 1], ylab = expression(mu))
abline(h = mu, col = 2)
plot.ts(Xsim[, 2], ylab = expression(tau))
abline(h = tau, col = 2)
par(mfrow = c(1, 1))

## Conjunta
plot(Xsim, main = "", xlab = expression(mu),
     ylab = expression(tau), ylim = range(X[, 2]))
abline(v = mu, h = tau, lty = 2, col = 2)

## Compara estatisticas
summary(Xsim)

## Compara amostra com distribuicoes teoricas (condicionais)
par(mfrow = c(1, 2))
hist(Xsim[,1], freq = FALSE, main = "", xlab = expression(mu))
C <- 1/(n * tau + tau0)
m <- C * (tau * ysum + mu0 * tau0)
curve(dnorm(x, m, sqrt(C)), add = TRUE, col = 2)
hist(Xsim[,2], freq = FALSE, main = "", xlab = expression(tau))
curve(dgamma(x, a + n/2, b + sum((y - mu)^2)/2),
      add = TRUE, col = 2)
par(mfrow = c(1, 1))

##----------------------------------------------------------------------
## Usando JAGS
library(runjags)

## Dados
datalist <- dump.format(list(y = y))
params <- c("mu", "tau")
inicial <- dump.format(list(mu = mean(y), tau = 1/var(y)))

## Modelo
mod <- "model{
for(i in 1:length(y)){
y[i] ~ dnorm(mu, tau)
}
mu ~ dnorm(0, 1)
tau ~ dgamma(0.001, 0.001)
}"

## Ajuste
m.jags <- run.jags(
    model = mod, monitor = params, data = datalist, inits = inicial,
    n.chains = 1, burnin = 5000, thin = 5, sample = 10000,
    mutate = list("prec2sd", vars = "tau")
    ## method = "parallel"
)

## Resultados
m.jags
names(m.jags)
plot(m.jags)

## Extrai cadeia
head(m.jags$mcmc)
str(m.jags$mcmc)
mu.jags <- m.jags$mcmc[[1]][, 1]
tau.jags <- m.jags$mcmc[[1]][, 2]
str(mu.jags)
str(tau.jags)
summary(mu.jags)
summary(tau.jags)
par(mfrow = c(1, 2))
plot(density(mu.jags))
plot(density(tau.jags))
par(mfrow = c(1, 1))

## Compara amostra com distribuicoes teoricas (condicionais)
par(mfrow = c(2, 2))
## Implementação manual
hist(Xsim[,1], freq = FALSE, main = "", xlab = expression(mu))
C <- 1/(n * tau + tau0)
m <- C * (tau * ysum + mu0 * tau0)
curve(dnorm(x, m, sqrt(C)), add = TRUE, col = 2)
hist(Xsim[,2], freq = FALSE, main = "", xlab = expression(tau))
curve(dgamma(x, a + n/2, b + sum((y - mu)^2)/2),
      add = TRUE, col = 2)
## JAGS
hist(mu.jags, freq = FALSE, main = "", xlab = expression(mu))
C <- 1/(n * tau + tau0)
m <- C * (tau * ysum + mu0 * tau0)
curve(dnorm(x, m, sqrt(C)), add = TRUE, col = 2)
hist(tau.jags, freq = FALSE, main = "", xlab = expression(tau))
curve(dgamma(x, a + n/2, b + sum((y - mu)^2)/2),
      add = TRUE, col = 2)
par(mfrow = c(1, 1))

##======================================================================
## Modelo de regressao

## Dados
n <- 20
x <- sort(runif(n, 0, 20))
epsilon <- rnorm(n, 0, 2.5)
y <- 2 + 0.5*x + epsilon
plot(x, y)

##----------------------------------------------------------------------
## Usando rjags
cat( "model {
      	for (i in 1:n){
		y[i] ~ dnorm(mu[i], tau)
		mu[i] <- b0 + b1 * x[i]
	}
	b0 ~ dnorm(0, .0001)
	b1 ~ dnorm(0, .0001)
	tau <- pow(sigma, -2)
	sigma ~ dunif(0, 100)
}", file="reglin.model")

## alternativa:
## tau ~dgamma(0.001, 0.001)
## sigma2 <- 1/tau

library(rjags)

## Iniciais
inis <- list(list(b0=0, b1=1, sigma=1),
             list(b0=2, b1=0.1, sigma=5))

jags <- jags.model(
    'reglin.model',
    data = list('x' = x,'y' = y,'n' = n),
    n.chains = 2,
    inits = inis,
    n.adapt = 100
)
class(jags)
jags

## Amostra da posterior
sam <- jags.samples(jags,
             c('b0', 'b1', 'sigma'),
             1000)
class(sam)
str(sam)

## Gera classe mais apropriada
sam <- coda.samples(jags,
             c('b0', 'b1', 'sigma', 'y'),
             1000)
class(sam)
str(sam)
plot(sam)

## Intervalo HPD
int <- HPDinterval(sam)
str(int)
int

##----------------------------------------------------------------------
## Usando runjags

datalist <- dump.format(list(y = y, x = x, n = n))
params <- c("b0", "b1", "tau")
inicial <- dump.format(list(b0 = 0, b1 = 1, sigma = 1))

## Modelo
mod <- "model {
      	for (i in 1:n){
		y[i] ~ dnorm(mu[i], tau)
		mu[i] <- b0 + b1 * x[i]
	}
	b0 ~ dnorm(0, .0001)
	b1 ~ dnorm(0, .0001)
	tau <- pow(sigma, -2)
	sigma ~ dunif(0, 100)
}"

## Ajuste
m.jags <- run.jags(
    model = mod, monitor = params, data = datalist, inits = inicial,
    n.chains = 1, burnin = 5000, thin = 5, sample = 10000,
    mutate = list("prec2sd", vars = "tau")
)

## Resultados
m.jags
quantile(lambda.sam, probs = c(.025, .5, .975))
quantile(beta.sam, probs = c(.025, .5, .975))
plot(m.jags)
