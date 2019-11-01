##======================================================================
## Exemplo 7.2 do Casella e Robert - Beta-binomial
## X | theta ~ Bin(n, theta)
## theta ~ Beta(a, b)

## Condicionais são
## X | theta ~ Bin(mu, theta)
## theta | X ~ Beta(x + a, n - x + b)

## Define constantes
N <- 5000
## Burnin
burn <- 1000
## Vetores para armazenar as amostras
T <- numeric(N)
X <- numeric(N)

## Aqui da para testar varios valores para ver o impacto da priori
a <- 3
b <- 7
curve(dbeta(x, a, b), from = 0, to = 1)

## Define valores
n <- 15
## Valores iniciais
T[1] <- rbeta(1, a, b)
X[1] <- rbinom(1, n, T[1])

## Amostrador de Gibbs
for (i in 2:N) {
    X[i] <- rbinom(1, n, T[i - 1])
    T[i] <- rbeta(1, X[i] + a, n - X[i] + b)
}

## Cadeias (mude os valores iniciais para ver convergencia)
par(mfrow = c(2, 1))
plot(X, type = "l")
plot(T, type = "l")
par(mfrow = c(1, 2))

## Correlacao entre os valores
acf(X)
acf(T)
par(mfrow = c(1, 1))
## Conjunta
plot(X, T)

## Integrando a conjunta em relação a theta, chega na marginal de X, que
## é uma Beta-Binomial
betabinom <- function(x, a, b, n) {
    choose(n, x) * (gamma(a + b)/(gamma(a) * gamma(b))) *
        ((gamma(x + a) * gamma(n - x + b))/
        gamma(a + b + n))
}

## Compara amostra com distribuicoes teoricas
par(mfrow = c(1, 2))
## A Beta-binomial é uma distribuição discreta
plot(prop.table(table(X)), type = "h")
points((0:14) + 0.2,
       betabinom(x = 0:14, a = a, b = b, n = n), type = "h", col = 2)
hist(T, freq = FALSE, main = "", xlab = expression(theta))
curve(dbeta(x, a, b), add = TRUE, col = 2)
par(mfrow = c(1, 1))

## Fazendo o thinning
x <- X[seq(1, length(X), 5)]
t <- T[seq(1, length(T), 5)]

## Correlacao entre os valores
par(mfrow = c(1, 2))
acf(x)
acf(t)
par(mfrow = c(2, 1))

## Cadeias (mude os valores iniciais para ver convergencia)
par(mfrow = c(2, 1))
plot(x, type = "l")
plot(t, type = "l")
par(mfrow = c(1, 1))

## Conjunta
plot(x, t)

## Compara amostra com distribuicoes teoricas
par(mfrow = c(1, 2))
## A Beta-binomial é uma distribuição discreta
plot(prop.table(table(x)), type = "h")
points((0:14) + 0.2,
       betabinom(x = 0:14, a = a, b = b, n = n), type = "h", col = 2)
hist(t, freq = FALSE, main = "", xlab = expression(theta))
curve(dbeta(x, a, b), add = TRUE, col = 2)
par(mfrow = c(1, 1))
##======================================================================

##======================================================================
## Exemplo de Rizzo, pg. 320
## Example 11.10 (Gibbs sampler: Bivariate distribution)

## Define constantes
N <- 5000
## Burnin
burn <- 1000
## Matriz para armazenar as amostras
X <- matrix(0, N, 2)

## Define parametros da Normal bivariada
rho <- -.75
mu1 <- 0
mu2 <- 2
sigma1 <- 1
sigma2 <- .5
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

## Gera a cadeia

## Valores iniciais
X[1, ] <- c(mu1, mu2)

for (i in 2:N) {
    x2 <- X[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] <- rnorm(1, m1, s1)
    x1 <- X[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] <- rnorm(1, m2, s2)
}

## Cadeias (mude os valores iniciais para ver convergencia)
matplot(X, type = "l")
## Correlacao entre os valores
acf(X[,1])
acf(X[,2])
## Conjunta
plot(X, main = "", xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(X[, 2]))

## Descarta os primeiros 1000 valores
b <- burn + 1
x <- X[b:N, ]

matplot(x, type = "l")
## Nao elimina o problema de autocorrelacao
acf(x[,1])
acf(x[,2])
## Conjunta
plot(x, main = "", xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(x[, 2]))

## Compara estatisticas
colMeans(x)
cov(x)
cor(x)
##======================================================================

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

## Define constantes
N <- 5000
## Burnin
burn <- 1000
## Matriz para armazenar as amostras
X <- matrix(0, N, 2)

## Simula valores
set.seed(2)
y <- rnorm(30, 0, 1)
mean(y)
var(y)

## Define parametros
ybar <- mean(y)
n <- length(y)
tau0 <- 1
## Aqui da para testar varios valores para ver o impacto da priori
a <- 2
b <- 2
curve(dgamma(x, a, b), from = 0, to = 5)

## Gera a cadeia

## Valores iniciais
X[1, ] <- c(0, 1)
for (i in 2:N) {
    tau.iter <- X[i - 1, 2]
    C <- 1/(n * tau.iter + tau0)
    m <- C * ybar
    X[i, 1] <- rnorm(1, m, sqrt(C))
    mu.iter <- X[i, 1]
    X[i, 2] <- rgamma(1, a + n/2, b + sum((y - mu.iter)^2)/2)
}

## Cadeias (mude os valores iniciais para ver convergencia)
matplot(X, type = "l")

## Correlacao entre os valores
acf(X[,1])
acf(X[,2])
## Conjunta
plot(X, main = "", xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(X[, 2]))
abline(v = 0, h = 1, lty = 2, col = 2)

## Descarta os primeiros 1000 valores
burn <- burn + 1
x <- X[burn:N, ]

matplot(x, type = "l")
acf(x[,1])
acf(x[,2])
## Conjunta
plot(x, main = "", xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(x[, 2]))
abline(v = 0, h = 1, lty = 2, col = 2)

## Compara estatisticas
colMeans(x)

## Compara amostra com distribuicoes teoricas
par(mfrow = c(1, 2))
hist(x[,1], freq = FALSE, main = "", xlab = expression(mu))
curve(dnorm(x, ybar/(n * 1 + tau0), sqrt(1/(n * 1 + tau0))),
      add = TRUE, col = 2)
hist(x[,2], freq = FALSE, main = "", xlab = expression(sigma^2))
curve(dgamma(x, a + n/2, b + sum((y - 0)^2)/2),
      add = TRUE, col = 2)
par(mfrow = c(1, 1))

##======================================================================

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

## Define constantes
N <- 5000
## Burnin
burn <- 1000
## Matriz para armazenar as amostras
X <- matrix(0, N, 2)

## Simula valores
set.seed(2)
y <- rnorm(30, 0, 1)

## Define parametros
ysum <- sum(y)
n <- length(y)

## Prioris
mu0 <- 0
tau0 <- 1
## Aqui da para testar varios valores para ver o impacto da priori
a <- 2
b <- 2
curve(dgamma(x, a, b), from = 0, to = 5)

## Gera a cadeia

## Valores iniciais
X[1, ] <- c(0, 1)
for (i in 2:N) {
    tau.iter <- X[i - 1, 2]
    C <- 1/(n * tau.iter + tau0)
    m <- C * (tau.iter * ysum + mu0 * tau0)
    X[i, 1] <- rnorm(1, m, sqrt(C))
    mu.iter <- X[i, 1]
    X[i, 2] <- rgamma(1, a + n/2, b + sum((y - mu.iter)^2)/2)
}

## Cadeias (mude os valores iniciais para ver convergencia)
matplot(X, type = "l")

## Correlacao entre os valores
acf(X[,1])
acf(X[,2])
## Conjunta
plot(X, main = "", xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(X[, 2]))
abline(v = 0, h = 1, lty = 2, col = 2)

## Descarta os primeiros 1000 valores
burn <- burn + 1
x <- X[burn:N, ]

matplot(x, type = "l")
acf(x[,1])
acf(x[,2])
## Conjunta
plot(x, main = "", xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(x[, 2]))
abline(v = 0, h = 1, lty = 2, col = 2)

## Compara estatisticas
colMeans(x)

## Compara amostra com distribuicoes teoricas
par(mfrow = c(1, 2))
hist(x[,1], freq = FALSE, main = "", xlab = expression(mu))
C <- 1/(n * tau + tau0)
m <- C * (tau * ysum + mu0 * tau0)
curve(dnorm(x, m, sqrt(C)), add = TRUE, col = 2)
hist(x[,2], freq = FALSE, main = "", xlab = expression(sigma^2))
curve(dgamma(x, a + n/2, b + sum((y - 0)^2)/2),
      add = TRUE, col = 2)
par(mfrow = c(1, 1))
##======================================================================
