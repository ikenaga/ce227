##======================================================================
## Aplicações de inferência bayesiana

##----------------------------------------------------------------------
## Pacotes
library(lattice)
library(latticeExtra)
library(gridExtra)
library(lme4)
library(rjags)
library(INLA)

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

##----------------------------------------------------------------------
## Efeito fixo
cat("model{
  for(i in 1:M){
     for(j in 1:N){
        y[i,j] ~ dnorm(mu[i], tau)
     }
     mu[i] ~ dnorm(50, 0.0001)
  }
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
}", file = "av04-q2.jags")

q2.dat <- list(y = t(y), M = Ng, N = Nobs)
q2.model <- jags.model(file="av04-q2.jags", data=q2.dat, n.chains=1)

q2.sam <- coda.samples(q2.model,
                       c("mu", "sigma"),
                       20000, thin=10)
summary(q2.sam)

post.jags <- q2.sam[[1]][, 1:10]
str(post.jags)

## Gráficos com as posteriores
post.jags.df <- stack(as.data.frame(post.jags))
str(post.jags.df)
p1 <- densityplot(~values, post.jags.df, groups = ind, auto.key = TRUE,
                  xlab = expression(mu))
p2 <- densityplot(~q2.sam[[1]][,11], xlab = expression(sigma^2))
grid.arrange(p1, p2, ncol = 1)

## ECDF das posteriores
ecdfplot(~values, post.jags.df, groups = ind, auto.key = TRUE,
         xlab = expression(mu^cpue))

## Resumo numérico das posteriores
respost <- function(x) {
    c(quantile(x, probs = c(.025, .5, .975)),
      media = mean(x), dp = sd(x))
}
(tab.res.jags <- t(cbind(apply(post.jags, 2, respost),
                         s2 = respost(q2.sam[[1]][,11]))))

## Compara intervalos de credibilidade
segplot(factor(1:10) ~
            tab.res.jags[1:10, 1] + tab.res.jags[1:10, 3],
        centers = tab.res.jags[1:10, 2], draw.bands = FALSE,
        xlab = expression(mu), ylab = "Áreas")

## Compara com lm
m <- lm(y ~ -1 + factor(grupo), q2.df)
summary(m)
tab.res.jags[1:10, 4:5]

##----------------------------------------------------------------------
## Efeito aleatório
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
  ## Coeficiente de correlação intra-classe (CCI)
  ## Pode ser calculado diretamente por aqui, ou calculado depois
  ## com as amostras das posteriores
  # CCI <- pow(delta, 2)/(pow(delta, 2) + pow(sigma, 2))
}", file = "av04-q2.jags")

q2.dat <- list(y = t(y), M = Ng, N = Nobs)
q2.model <- jags.model(file="av04-q2.jags", data=q2.dat, n.chains=1)

q2.sam <- coda.samples(q2.model,
                       c("mu", "sigma", "delta", "theta"),
                       20000, thin=10)
summary(q2.sam)

## Coeficiente de correlação intra-classe, calculado usando as amostras
## das posteriores
post.delta <- q2.sam[[1]][, 1]
post.sigma <- q2.sam[[1]][, 12]
CCI <- post.delta^2/(post.delta^2 + post.sigma^2)
hist(CCI)
respost(CCI)

## Posteriores dos mu[i]
post.jags <- q2.sam[[1]][, 2:11]
str(post.jags)

## Gráficos com as posteriores
post.jags.df <- stack(as.data.frame(post.jags))
str(post.jags.df)
p1 <- densityplot(~values, post.jags.df, groups = ind, auto.key = TRUE,
                  xlab = expression(mu))
p1
p2 <- densityplot(~q2.sam[[1]][,12], xlab = expression(sigma))
p3 <- densityplot(~q2.sam[[1]][,13], xlab = expression(theta))
p4 <- densityplot(~q2.sam[[1]][,1], xlab = expression(delta))
grid.arrange(p2, p3, p4, ncol = 1)

## ECDF das posteriores
ecdfplot(~values, post.jags.df, groups = ind, auto.key = TRUE,
         xlab = expression(mu^cpue))

## Resumo numérico das posteriores
(tab.res.jags <- t(cbind(apply(post.jags, 2, respost),
                         s2 = respost(q2.sam[[1]][,12]))))

## Compara intervalos de credibilidade
segplot(factor(1:10) ~
            tab.res.jags[1:10, 1] + tab.res.jags[1:10, 3],
        centers = tab.res.jags[1:10, 2], draw.bands = FALSE,
        xlab = expression(mu), ylab = "Áreas")

## Compara com lme
library(lme4)
(q2.lmer <- lmer(y ~ 1|grupo, data=q2.df, REML=FALSE))
t(ranef(q2.lmer)$grupo)
t(ranef(q2.lmer)$grupo +  fixef(q2.lmer)[[1]])
tab.res.jags[1:10, 4:5]


##----------------------------------------------------------------------
## INLA
library(INLA)
q2.inla <- inla(y ~  f(grupo, model="iid"), family="gaussian",
                data=q2.df)
summary(q2.inla)
## passando de precisão para variância:
sqrt(1/q2.inla$summary.hyperpar[,1])
post.sigma = inla.tmarginal(function(x) sqrt(1/x),
                            q2.inla$marginals.hyperpar[[1]])
post.delta = inla.tmarginal(function(x) sqrt(1/x),
                            q2.inla$marginals.hyperpar[[2]])
rbind(delta=inla.zmarginal(post.delta, T),
      sigma=inla.zmarginal(post.sigma, T))


## Efeitos aleatórios: médias e modas por grupos.
rbind(medias.grupos=q2.inla$summary.fix[1,1] +
          q2.inla$summary.random$grupo$mean,
      modas.grupos=q2.inla$summary.fix[1,6] +
          q2.inla$summary.random$grupo$mode)

## Comparação dos valores das variâncias dos grupos e residual.
par(mfrow=c(1,2))
plot(density(unlist(q2.sam[,"delta",drop=T])), main="",
     xlab=expression(delta),
     ylab=expression(group("[",paste(delta,"|",y),"]")), ylim=c(0, 0.5))
lines(post.delta, lty=2, col=2)
abline(v=c(as.data.frame(VarCorr(q2.lmer))$sdcor[1], delta),
       col=3:4, lty=3:4)
legend("topright", c("JAGS","INLA", "LMER", "Verd."), col=1:4, lty=1:4)
plot(density(unlist(q2.sam[,"sigma",drop=T])), main="",
     xlab=expression(sigma),
     ylab=expression(group("[",paste(sigma,"|",y),"]")), ylim=c(0, 1.5))
lines(post.sigma, lty=2, col=2)
abline(v=c(as.data.frame(VarCorr(q2.lmer))$sdcor[2], sigma),
       col=3:4, lty=3:4)
legend("topright", c("JAGS","INLA", "LMER", "Verd."), col=1:4, lty=1:4)
par(mfrow=c(1,1))

## Médias dos grupos
(q2.grupos <- data.frame(lmer = drop(t(ranef(q2.lmer)$grupo +
                                                     fixef(q2.lmer)[[1]])),
                         JAGS = summary(q2.sam)[[1]][2:11,1],
                         INLA = q2.inla$summary.fix[1,1] +
                             q2.inla$summary.random$grupo$mean,
                         verd = mus))

matplot(1:10, q2.grupos, type="p", xlab="grupo",
        ylab=expression(mu[i]))
legend("topleft", c("lmer","JAGS","INLA", "Verd"),
       pch=as.character(1:4), col=1:4)

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

(tab.res <- t(cbind(apply(emuGpost, 2, respost), s2 = respost(varG))))

## Compara intervalos de credibilidade
segplot(factor(c(9, 14, 15)) ~ tab.res[1:3, 1] + tab.res[1:3, 3],
        centers = tab.res[1:3, 2], draw.bands = FALSE,
        xlab = expression(mu^cpue), ylab = "Áreas")


##----------------------------------------------------------------------
## Usando o JAGS

## A forma de especificar o modelo como está abaixo parece seguir a
## lógica do exemplo anterior. O problema aqui é que o número de
## observações por grupo é diferente (desbalanceado), portanto a matriz
## de dados deve ser construída apropriadamente (com NAs para
## preencher). No entanto isso pode ser complicado para bases grandes.
## Mais abaixo é mostrada uma solução sem precisar mexer na estrutura
## dos dados.
## ## Modelo (aqui apenas para mostrar que dessa forma não funciona)
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

## Usando o rjags

## Definição do modelo mais geral. Aqui, g é o identificador de grupo e
## i o identificador da observação. Portanto
## mu[g[i]]
## é uma forma de identificar a observação i do grupo g. Dessa forma, os
## dados podem ser passados da maneira como estão no data frame.
## Lembrando que isso é necessário pois o número de onservações por
## grupo é diferente nesse caso. De qualquer maneira, o modelo
## especificado dessa forma é mais geral e funciona também no caso
## balanceado.
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

## Posteriores na escala log
post.jags <- sam[[1]][, 1:3]
str(post.jags)
## Posteriores na escala original
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
ecdfplot(~values, epost.jags.df, groups = ind, auto.key = TRUE,
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


## Compara com lm
(tab.res.jags <- t(cbind(apply(post.jags, 2, respost),
                         s2 = respost(sam[[1]][,4]))))

m <- lm(lcpue ~ -1 + factor(area), dados)
summary(m)
anova(m)

## Diferenca entre médias
diff.post <- emuGpost[,3] - emuGpost[,1]
hist(diff.post); abline(v = 0, col = 2)
respost(diff.post)
sum(diff.post > 0)/length(diff.post)
