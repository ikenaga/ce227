## Dados
n <- 250
y <- 32

## Proporção amostral
(theta.hat <- y/n)

## Intervalo de confiança para a proporção

## Com alpha = 0.05, define valor crítico
(zc <- qnorm(.025))
## Calcula erro padrão
(ep <- sqrt((.5 * .5)/n))
## Determina intervalo
theta.hat - zc * ep * c(-1, 1)

## Teste de hipótese
## H0: theta = 0.15
## Ha: theta < 0.15
theta0 <- 0.15

## Estatistica de teste
(zcalc <- (theta.hat - theta0)/sqrt((theta0 * (1-theta0))/n))
## Com alpha = 0.05, o valor cítico é
(zcrit <- qnorm(.025))
ifelse(abs(zcalc) > abs(zcrit), "Rejeita H0", "Não rejeita H0")

## p-valor
pnorm(zcalc)

## Por simulação, sob theta0

## Gera k amostras com n = 250 e assumindo que theta0 é verdadeiro
set.seed(12)
k <- 10000
am <- rbinom(k, size = 250, prob = theta0)

## Proporção amostral
theta.hat.am <- am/n
hist(theta.hat.am)

## Padroniza a distribuição das propoções amostrais
z <- (theta.hat.am - theta0)/sqrt((theta0 * (1-theta0))/n)
hist(z, freq = FALSE)
## Aproximação pela normal
curve(dnorm, -3, 3, add = TRUE)
abline(v = zcalc, col = 2)

## Proporcao da amostra abaixo de theta0 ~ "p-valor"
sum(z < zcalc)/k
