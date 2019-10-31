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
## mu ~ N(0, s^2)
## sigma^2 ~ Gama(a,b)

## Condicionais são
## mu | sigma^2 ~ N(m, C)
## sigma^2 | mu ~ Gama(a + n/2, b + 1/2 sum((y_i - mu)^2))
## Onde:
## m = C\bar{y}
## C = 1/(n\tau + s^{-2})

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
ybar <- mean(y)
n <- length(y)
s2 <- 1
## Aqui da para testar varios valores para ver o impacto da priori
a <- 2
b <- 2
curve(dgamma(x, a, b), from = 0, to = 5)

## Gera a cadeia

## Valores iniciais
X[1, ] <- c(0, 1)
for (i in 2:N) {
    X[i, 1] <- rnorm(1, ybar/(n * X[i - 1, 2] + 1/s2),
                     1/(n * X[i - 1, 2] + 1/s2))
    X[i, 2] <- rgamma(1, a + n/2, b + sum((y - X[i,1])^2)/2)
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
burn <- burn + 1
x <- X[burn:N, ]

matplot(x, type = "l")
acf(x[,1])
acf(x[,2])
## Conjunta
plot(x, main = "", xlab = bquote(X[1]),
     ylab = bquote(X[2]), ylim = range(x[, 2]))

## Compara estatisticas
colMeans(x)

## Compara amostra com distribuicoes teoricas
hist(x[,1], freq = FALSE)
curve(dnorm(x, ybar/(n * 1 + 1/s2), 1/(n * 1 + 1/s2)),
      add = TRUE, col = 2)
hist(x[,2], freq = FALSE)
curve(dgamma(x, a + n/2, b + sum((y - 0)^2)/2),
      add = TRUE, col = 2)
