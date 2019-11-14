##======================================================================
## Exemplo de Kinas e Andrade (2010)

dados <- read.table("../dados/martelo.csv", header = TRUE,
                    sep = ",", dec = ".")
str(dados)
plot(dados)

dados$cpue <- dados$capkg * 1000/dados$nanz
dados$lcpue <- log(dados$cpue)

par(mfrow = c(1, 2))
boxplot(cpue ~ area, dados)
boxplot(lcpue ~ area, dados)
par(mfrow = c(1, 1))

fn <- function(x) {
    c(media = mean(x), mediana = median(x), var = var(x))
}
aggregate(cpue ~ area, dados, FUN = fn)
aggregate(lcpue ~ area, dados, FUN = fn)

## Estatísticas necessárias
fnsq <- function(x) sum((x - mean(x))^2)
(ygbar <- with(dados, tapply(lcpue, area, mean)))
(sqg <- with(dados, tapply(lcpue, area, fnsq)))
(ng <- with(dados, tapply(lcpue, area, length)))
(N <- sum(ng))
(G <- length(ng))
(se2G <- sum(sqg)/(N - G))

## Simula da posterior conjunta
Vmu <- diag(1/ng)
MASS::fractions(Vmu)
m <- 3000
tauG <- rgamma(m, (N - G)/2, se2G * (N - G)/2)
varG <- 1/tauG
muGpost <- matrix(NA, nrow = m, ncol = G,
                  dimnames = list(NULL, c(9, 14, 15)))
for(i in 1:m) {
    muGpost[i, ] <- MASS::mvrnorm(1, as.numeric(ygbar), varG[i] * Vmu)
}
matplot(muGpost, type = "l")

## Converte a posterior para a escala original
emuGpost <- exp(muGpost)

## Gráficos com as posteriores
library(gridExtra)
library(latticeExtra)
muGpost.df <- stack(as.data.frame(muGpost))
p1 <- densityplot(~values, muGpost.df, groups = ind, auto.key = TRUE,
                  xlab = expression(mu))
p2 <- densityplot(~varG, xlab = expression(sigma^2))
emuGpost.df <- stack(as.data.frame(emuGpost))
p3 <- densityplot(~values, emuGpost.df, groups = ind, auto.key = TRUE,
            xlab = expression(mu^cpue))
grid.arrange(p1, p2, p3, ncol = 1)

## ECDF das posteriores
ecdfplot(~values, emuGpost.df, groups = ind, auto.key = TRUE,
         xlab = expression(mu^cpue))

## Resumo numérico das posteriores
respost <- function(x) {
    c(quantile(x, probs = c(.025, .5, .975)),
      media = mean(x), dp = sd(x))
}
(tab.res <- t(cbind(apply(emuGpost, 2, respost), s2 = respost(varG))))

## Compara intervalos de credibilidade
segplot(factor(c(9, 14, 15)) ~ tab.res[1:3, 1] + tab.res[1:3, 3],
        centers = tab.res[1:3, 2], draw.bands = FALSE,
        xlab = expression(mu^cpue), ylab = "Áreas")


##----------------------------------------------------------------------
## Usando o JAGS
library(runjags)

## Os dados tem que ser passados como uma matriz com i linhas (grupos) e
## j colunas (observacoes)
## da <- dados[, c("lcpue", "area")]
## da <- da[order(da$area), ]
## ## Modelo
## mod <- "model {
##  	for (i in 1:G){
##     for(j in 1:ng){
## 		  y[i,j] ~ dnorm(mu[i], tau)
##   	}
##   mu[i] ~ dnorm(ygbar[i], tauD)
##   }
##   tau <- pow(sigma, -2)
##   sigma ~ dunif(0, 100)
## }"


datalist <- dump.format(list(y = dados$lcpue,
                             g = as.factor(dados$area),
                             ygbar = as.numeric(ygbar),
                             N = N, G = G,
                             tauD = 1/se2G))
params <- c("mu", "sigma")
inicial <- dump.format(list(sigma = 1))

## Modelo
mod <- "model {
 	for (i in 1:N){
		  y[i] ~ dnorm(mu[g[i]], tau)
 	}
  for(j in 1:G){
    mu[j] ~ dnorm(ygbar[j], tauD)
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
}"

## Ajuste
m.jags <- run.jags(
    model = mod, monitor = params, data = datalist, inits = inicial,
    n.chains = 1, burnin = 5000, thin = 5, sample = 10000
)
## failed.jags(c('model','data','inits'))

## Resultados
m.jags
quantile(lambda.sam, probs = c(.025, .5, .975))
quantile(beta.sam, probs = c(.025, .5, .975))
plot(m.jags)

library(rjags)

cat("model {
 	for (i in 1:N){
		  y[i] ~ dnorm(mu[g[i]], tau)
 	}
  for(j in 1:G){
    mu[j] ~ dnorm(ygbar[j], tauD)
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
}", file = "anova-fix.jags")

datalist <- list(y = dados$lcpue,
                 g = as.factor(dados$area),
                 ygbar = as.numeric(ygbar),
                 N = N, G = G,
                 tauD = 1/se2G)
m.rjags <- jags.model(file = "anova-fix.jags",
                      data = datalist, n.chains = 1)

sam <- coda.samples(m.rjags,
                    c("mu", "sigma"), 20000, thin = 10)
summary(sam)

post.jags <- sam[[1]][, 1:3]
str(post.jags)
epost.jags <- exp(sam[[1]][, 1:3])
str(epost.jags)

## Gráficos com as posteriores
post.jags.df <- stack(as.data.frame(post.jags))
str(post.jags.df)
p1 <- densityplot(~values, post.jags.df, groups = ind, auto.key = TRUE,
                  xlab = expression(mu))
p2 <- densityplot(~sam[[1]][,4], xlab = expression(sigma^2))
epost.jags.df <- stack(as.data.frame(epost.jags))
p3 <- densityplot(~values, epost.jags.df, groups = ind, auto.key = TRUE,
            xlab = expression(mu^cpue))
grid.arrange(p1, p2, p3, ncol = 1)

## ECDF das posteriores
ecdfplot(~values, post.jags.df, groups = ind, auto.key = TRUE,
         xlab = expression(mu^cpue))

## Resumo numérico das posteriores
(tab.res.jags <- t(cbind(apply(epost.jags, 2, respost),
                         s2 = respost(sam[[1]][,4]))))

## Compara intervalos de credibilidade
segplot(factor(c(9, 14, 15)) ~
            tab.res[1:3, 1] + tab.res[1:3, 3],
        centers = tab.res[1:3, 2], draw.bands = FALSE,
        xlab = expression(mu^cpue), ylab = "Áreas")

segplot(factor(c(9, 14, 15)) ~
            tab.res[1:3, 1] + tab.res[1:3, 3],
        centers = tab.res[1:3, 2], draw.bands = FALSE,
        col = 1,
        xlab = expression(mu^cpue), ylab = "Áreas") +
    as.layer(
        segplot(factor(c(9, 14, 15)) ~
                    tab.res.jags[1:3, 1] + tab.res.jags[1:3, 3],
                col = 2, lty = 2,
                centers = tab.res.jags[1:3, 2], draw.bands = FALSE,
                xlab = expression(mu^cpue), ylab = "Áreas")
    )

## Diferenca entre medias do Kinas

##======================================================================
## Exemplo prova PJ
set.seed(22701)
Ng <- 10
Nobs <- 5
N <- Ng*Nobs
delta <- 5
sigma <- 2
mus <- rnorm(Ng, mean=50, sd=delta)
y <- matrix(rnorm(N, mean=mus, sd=sigma), ncol=Ng, byrow=TRUE)
q2.df <- data.frame(y = as.vector(y), grupo = rep(1:Ng, each=Nobs))

require(lme4)
(q2.lmer <- lmer(y ~ 1|grupo, data=q2.df, REML=FALSE))
t(ranef(q2.lmer)$grupo)
t(ranef(q2.lmer)$grupo +  fixef(q2.lmer)[[1]])

cat("model{
  for(i in 1:M){
     for(j in 1:N){
        y[i,j] ~ dnorm(mu[i], tau)
     }
     mu[i] ~ dnorm(theta, tauD)
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
  theta ~ dnorm(0, .001)
  tauD <- pow(delta, -2)
  delta ~ dunif(0, 100)
}", file="av04-q2.jags")

q2.dat <- list(y = t(y), M = Ng, N = Nobs)
q2.model <- jags.model(file="av04-q2.jags", data=q2.dat, n.chains=3)

q2.sam <- coda.samples(q2.model,
                       c("mu", "sigma", "delta", "theta"),
                       20000, thin=10)
summary(q2.sam)
