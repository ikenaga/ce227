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
