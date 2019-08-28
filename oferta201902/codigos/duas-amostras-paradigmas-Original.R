##
## Comandos utilizados para gerar slides comparando paradigmas para inferência
## Utilizando a inferẽncia sobre a comparação de duas amostras como exemplo 
## - frequentista (distribuição amostral)
## - verossimilhança (função de verossimilhança) 
## - bayesiano (distribuição à posteriori)
##
## Responsável: Paulo Justiniano Ribeiro Junior


##
## Inferência frequentista - baseada na distribuição amostral
##

## ------------------------------------------------------------------------
source("dados/mandible.R")
#source("http://www.leg.ufpr.br/~paulojus/dados/mandible.R")
mand <- reshape(mandible, varying=list(1:2), direction="long")[,-3]
names(mand) <- c("sex","mandible")
rownames(mand) <- 1:20
mand$sex <- factor(mand$sex, labels=c("female","male"))

## ----echo=F--------------------------------------------------------------
mandible


## ----fig.width=10, fig.height=5, out.width = "0.7\\textwidth"------------
par(mfrow=c(1,2), mar=c(3,2.5, 1, 0), mgp=c(1.8, 0.8, 0))
with(mand, {
    hist(mandible, prob=T, main="", xlab="altura")
    lines(density(mandible))
})
    rug(jitter(mandible$females), col=2, lwd=2)
    rug(jitter(mandible$male), col=4, lwd=2)
    boxplot(mandible)
    stripchart(mandible, vertical=T, add=T, method="jitter", pch=19, cex=0.5, jitter=0.15, col=4) 
#    MASS:::boxcox(mandible~sex, data=mand)


## ----echo=TRUE-----------------------------------------------------------
with(mandible, t.test(females, male, var.equal=TRUE))


## ----echo=TRUE-----------------------------------------------------------
with(mandible, t.test(females, male, var.equal=FALSE))


## ----echo=TRUE-----------------------------------------------------------
with(mandible, t.test(females, male, paired=TRUE))


## ------------------------------------------------------------------------
d.obs <- with(mandible, mean(females) - mean(male))
source("mctest.R")


## ----echo=F,fig.width=5, fig.height=4, out.width = "0.5\\textwidth"------
mctest(mandible$f, mandible$m, paired=T)


## ----echo=TRUE-----------------------------------------------------------
require(stats4)
ll <- function(mu, theta, lsigma, am1, am2){
    sigma <- exp(lsigma)
    l1 <- sum(dnorm(am1, m=mu, sd=sigma, log=T))    
    l2 <- sum(dnorm(am2, m=mu+theta, sd=sigma, log=T))
    return(-(l1+l2))
}    


## ----echo=TRUE-----------------------------------------------------------
fit <- mle(ll, start=list(mu=110, theta=0, lsigma=log(10)), 
           fixed=list(am1 = mandible$fem, am2 = mandible$male))


## ----echo=TRUE-----------------------------------------------------------
summary(fit)
confint(fit, level=0.95)
prof <- profile(fit)


## ----results="hide",fig.width=7, fig.height=3, out.width = "0.99\\textwidth"----
par(mfrow=c(1,3), mar=c(2.5,2.5,0.3, 0.3), mgp=c(1.7, 0.7, 0))
plot(prof, level=c(0.5, 0.8, 0.95, 0.99))


## ----echo=TRUE-----------------------------------------------------------
L1 <- lm(mandible~1, data=mand)
L2 <- lm(mandible~sex, data=mand)
c(logLik(L1), logLik(L2))
anova(L1, L2)


## ----echo=TRUE-----------------------------------------------------------
logLik(L1)
logLik(L2)
##
-2 * logLik(L1) + 2 * 2
-2 * logLik(L2) + 2 * 3
c(AIC(L1), AIC(L2))
##
-2 * logLik(L1) + log(20) * 2
-2 * logLik(L2) + log(20) * 3
c(BIC(L1), BIC(L2))


## ------------------------------------------------------------------------
require(MCMCpack)


## ----echo=TRUE-----------------------------------------------------------
require(MCMCpack)
Btt <- MCMCregress(mandible ~ sex, dat=mand)
summary(Btt)


## ----results="hide",fig.width=7, fig.height=3, out.width = "0.99\\textwidth"----
par(mfrow=c(1,3), mar=c(2.5,2.5,0.3, 0.3), mgp=c(1.7, 0.7, 0))
plot(Btt, trace=F)


## ------------------------------------------------------------------------


## ----alturafit, echo=FALSE, results="hide"-------------------------------
(alt.n <- with(quest, MASS:::fitdistr(altura, densfun="normal")))
(alt.ln <- with(quest, MASS:::fitdistr(altura, densfun="lognormal")))
(alt.g <- with(quest, MASS:::fitdistr(altura, densfun="gamma")))
(alt.bc <- with(quest, MASS:::fitdistr(-((1/altura)-1), densfun="normal")))


## ----pesofit, echo=FALSE, results="hide"---------------------------------
(peso.n <- with(quest, MASS:::fitdistr(peso, densfun="normal")))
(peso.ln <- with(quest, MASS:::fitdistr(peso, densfun="lognormal")))
(peso.g <- with(quest, MASS:::fitdistr(peso, densfun="gamma")))
(peso.bc <- with(quest, MASS:::fitdistr(-((1/peso)-1), densfun="normal")))


## ------------------------------------------------------------------------
attach(quest)


## ----fig.width=8, fig.height=5.5, out.width = "0.9\\textwidth"-----------
par(mfrow=c(2,2), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
hist(altura, prob=T, main=""); lines(density(altura))
curve(dnorm(x, m=mean(altura), sd=sd(altura)), from=140, to=200, col=2, add=T)
plot(ecdf(altura))
curve(pnorm(x, m=mean(altura), sd=sd(altura)), from=140, to=200, col=2, add=T)
#MASS:::boxcox(altura~1, lambda=seq(-3.5, 1, 1/10))
hist(peso, prob=T, main="", ylim=c(0,0.035)); lines(density(peso))
curve(dnorm(x, m=mean(peso), sd=sd(peso)), from=40, to=140, col=2, add=T)
plot(ecdf(peso))
curve(pnorm(x, m=mean(peso), sd=sd(peso)), from=40, to=140, col=2, add=T)
#MASS:::boxcox(peso~1, lambda=seq(-3, 1, 1/10))


## ----fig.width=8, fig.height=5.5, out.width = "0.9\\textwidth"-----------
alturaBC <- -((1/altura)-1)
pesoBC <- -((1/peso)-1)
par(mfrow=c(2,3), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
MASS:::boxcox(altura~1, lambda=seq(-3.5, 1, 1/10))
hist(alturaBC, prob=T, main="", ylim=c(0,1150)); lines(density(-((1/altura)-1)))
curve(dnorm(x, m=mean(alturaBC), sd=sd(alturaBC)), from=0.9930, 
      to=0.9960, col=2, add=T)
plot(ecdf(alturaBC))
curve(pnorm(x, m=mean(alturaBC), sd=sd(alturaBC)), from=0.9930, to=0.9960, col=2, add=T)
MASS:::boxcox(peso~1, lambda=seq(-3, 1, 1/10))
hist(pesoBC, prob=T, main="", ylim=c(0, 130)); lines(density(-((1/peso)-1)))
curve(dnorm(x, m=mean(pesoBC), sd=sd(pesoBC)), from=0.970, to=1, col=2, add=T, n=1000)
plot(ecdf(pesoBC))
curve(pnorm(x, m=mean(pesoBC), sd=sd(pesoBC)), from=0.970, to=1, col=2, add=T)


## ----fig.width=8, fig.height=6, out.width = "0.9\\textwidth"-------------
par(mfrow=c(2,2), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
plot(density(altura), main="altura", ylim=c(0, 0.045), xlab="")
curve(dnorm(x, mean=alt.n$est[1], sd=alt.n$est[2]), from=140, to=200, col=2, add=T)
curve(dlnorm(x, meanlog=alt.ln$est[1], sdlog=alt.ln$est[2]), from=140, to=200, col=3, add=T)
curve(dgamma(x, shape=alt.g$est[1], rate=alt.g$est[2]), from=140, to=200, col=4, add=T)
legend("topright", c("Normal","LogNormal","Gamma"), col=2:4, lty=1, cex=0.8)
##
plot(density(peso), main="peso", ylim=c(0, 0.035), xlab="")
curve(dnorm(x, mean=peso.n$est[1], sd=peso.n$est[2]), from=30, to=150, col=2, add=T)
curve(dlnorm(x, meanlog=peso.ln$est[1], sdlog=peso.ln$est[2]), from=30, to=150, col=3, add=T)
curve(dgamma(x, shape=peso.g$est[1], rate=peso.g$est[2]), from=30, to=140, col=4, add=T)
legend("topright", c("Normal","LogNormal","Gamma"), col=2:4, lty=1, cex=0.8)
##
plot(ecdf(altura))
curve(pnorm(x, mean=alt.n$est[1], sd=alt.n$est[2]), from=140, to=200, col=2, add=T)
curve(plnorm(x, meanlog=alt.ln$est[1], sdlog=alt.ln$est[2]), from=140, to=200, col=3, add=T)
curve(pgamma(x, shape=alt.g$est[1], rate=alt.g$est[2]), from=140, to=200, col=4, add=T)
legend("topleft", c("Normal","LogNormal","Gamma"), col=2:4, lty=1, cex=0.8)
plot(ecdf(peso))
curve(pnorm(x, mean=peso.n$est[1], sd=peso.n$est[2]), from=30, to=150, col=2, add=T)
curve(plnorm(x, meanlog=peso.ln$est[1], sdlog=peso.ln$est[2]), from=30, to=150, col=3, add=T)
curve(pgamma(x, shape=peso.g$est[1], rate=peso.g$est[2]), from=30, to=140, col=4, add=T)
legend("topleft", c("Normal","LogNormal","Gamma"), col=2:4, lty=1, cex=0.8)


## ----echo=TRUE,results="markup"------------------------------------------
(alt.n <- MASS:::fitdistr(altura, densfun="normal"))
(alt.ln <- MASS:::fitdistr(altura, densfun="lognormal"))
(alt.g <- MASS:::fitdistr(altura, densfun="gamma"))
(alt.bc <- MASS:::fitdistr(-((1/altura)-1), densfun="normal"))


## ----echo=TRUE,results="markup"------------------------------------------
(altlL <- c(logLik(alt.n), logLik(alt.ln), logLik(alt.g), 
            (-1-1)*sum(log(altura)) + logLik(alt.bc)))
max(altlL) - altlL 


## ----echo=TRUE,results="markup"------------------------------------------
(peso.n <- MASS:::fitdistr(peso, densfun="normal"))
(peso.ln <- MASS:::fitdistr(peso, densfun="lognormal"))
(peso.g <- MASS:::fitdistr(peso, densfun="gamma"))
(peso.bc <- MASS:::fitdistr(-((1/peso)-1), densfun="normal"))


## ----echo=TRUE,results="markup"------------------------------------------
(pesolL <- c(logLik(peso.n), logLik(peso.ln), logLik(peso.g), 
             (-1-1)*sum(log(peso))+logLik(peso.bc)))
max(pesolL) - pesolL


## ----fig.width=8, fig.height=6, out.width = "0.9\\textwidth"-------------
par(mfrow=c(2,3), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
hist(altura, prob=T, main=""); lines(density(altura))
MASS:::boxcox(altura~1, lambda=seq(-3.5, 1, 1/10))
MASS:::boxcox(altura~sexo, lambda=seq(-3.5, 1, 1/10))
hist(peso, prob=T, main=""); lines(density(peso))
MASS:::boxcox(peso~1, lambda=seq(-3, 1, 1/10))
MASS:::boxcox(peso~sexo, lambda=seq(-3, 1, 1/10))


## ----echo=TRUE,results="markup"------------------------------------------
## Ajustando distribuições
alt.n <- glm(altura ~ 1, family="gaussian", data=quest)
alt.g <- glm(altura ~ 1, family=Gamma(), data=quest)
alt.bc <- glm(-((1/altura)-1) ~ 1, family="gaussian", data=quest)
## Ajustando distribuições com médias diferentes para  masculino e feminino
alt.n1 <- glm(altura ~ sexo, family="gaussian", data=quest)
alt.g1 <- glm(altura ~ sexo, family=Gamma(), data=quest)
alt.bc1 <- glm(-((1/altura)-1) ~ sexo, family="gaussian", data=quest)
## Verossimilhanças dos ajustes
fits <- rbind(c(logLik(alt.n), logLik(alt.g), 
                (-1-1)*sum(log(altura)) + logLik(alt.bc), 2),
              c(logLik(alt.n1), logLik(alt.g1), 
                (-1-1)*sum(log(altura)) + logLik(alt.bc1), 3)
              )
dimnames(fits) <- list(c("Modelo 1","Modelo 2"), 
                       c("Normal","Gamma","Box-Cox","npar"))
fits


## ----results="hide",echo=FALSE,fig.width=8,fig.height=6,out.width="0.98\\textwidth"----
par(mfrow=c(2,2), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
#names(summary(alt.n1))
#coef(alt.n1)
alt.n1.med <- cumsum(coef(alt.n1))
alt.n1.sd <- with(summary(alt.n1), sqrt(deviance/df.residual))
#altlm1 <- lm(altura ~sexo)
#names(summary(altlm1))
#summary(altlm1)$sigma
propFM <- prop.table(table(sexo))
plot(density(altura), main="altura", ylim=c(0, 0.045), xlab="")
curve(0.6*dnorm(x, m=alt.n1.med[1], sd=alt.n1.sd), from=140, to=200, col=4, add=T, lwd=1, lty=2)    
curve(0.4*dnorm(x, m=alt.n1.med[2], sd=alt.n1.sd), from=140, to=200, col=4, add=T, lwd=1, lty=2)    
dalt.n1 <- function(x){
    propFM[1] * dnorm(x, m=alt.n1.med[1], sd=alt.n1.sd) +
        propFM[2] * dnorm(x, m=alt.n1.med[2], sd=alt.n1.sd)
}
curve(dnorm(x, mean=coef(alt.n), sd=with(summary(alt.n), sqrt(deviance/df.residual))), 
      from=140, to=200, col=2, add=T, lwd=2)
curve(dalt.n1(x), from=140, to=210, add=T, col=4, lwd=2)
legend("topright", c("Normal 1","Normal 2"), lty=1, col=c(2,4), lwd=2)

plot(ecdf(altura))
palt.n1 <- function(x){
    propFM[1] * pnorm(x, m=alt.n1.med[1], sd=alt.n1.sd) +
        propFM[2] * pnorm(x, m=alt.n1.med[2], sd=alt.n1.sd)
}
curve(pnorm(x, mean=coef(alt.n), sd=with(summary(alt.n), sqrt(deviance/df.residual))), 
      from=140, to=200, col=2, add=T, lwd=2)
curve(palt.n1(x), from=140, to=210, add=T, col=4, lwd=2)
legend("topleft", c("Normal 1","Normal 2"), lty=1, col=c(2,4), lwd=2)
#@ 
#\end{frame}
#
#
#\begin{frame}[fragile]
#  \frametitle{Modelos para \code{altura}}
#
#  \red{Comparando ajustes para o modelo com distribuição Gama}
#<<results="hide",echo=FALSE,fig.width=8,fig.height=4,out.width="0.9\\textwidth"#>>=
##
## Distribuição Y ~Gama(shape = alpha, rate = beta)
## Ajuste Gama GLM (link = "inverse")
## Y ~Gama (\mu, \phi)
## Para y ~1:
## \mu = 1/beta_0 
## alpha = 1/phi = 1/summary(fit)$dispersion
## beta = beta_0/phi = coef(fit)/summary(fit)$dispersion
## Modelo 1 ~1
alphaG <- 1/summary(alt.g)$dispersion
betaG <- coef(alt.g) * alphaG
## Modelo 2 ~ sexo
alphaG1 <- 1/summary(alt.g1)$dispersion
betaGF <- coef(alt.g1)[1] * alphaG1
betaGM <- sum(coef(alt.g1)) * alphaG1

propFM <- prop.table(table(sexo))
#par(mfrow=c(1,2), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
plot(density(altura), main="altura", ylim=c(0, 0.045), xlab="")
curve(0.6*dgamma(x, sh=alphaG1, rate=betaGF), from=140, to=200, col=4, add=T, lwd=1, lty=2)    
curve(0.4*dgamma(x, sh=alphaG1, rate=betaGM), from=140, to=200, col=4, add=T, lwd=1, lty=2)    
dalt.g1 <- function(x){
    propFM[1] * dgamma(x, sh=alphaG1, rate=betaGF) +
        propFM[2] * dgamma(x, sh=alphaG1, rate=betaGM) 
}
curve(dgamma(x,sh=alphaG, rate=betaG), from=140, to=200, col=2, add=T, lwd=2)
curve(dalt.g1(x), from=140, to=210, add=T, col=4, lwd=2)
legend("topright", c("Gama 1","Gama 2"), lty=1, col=c(2,4), lwd=2)
##
plot(ecdf(altura))
palt.g1 <- function(x){
    propFM[1] * pgamma(x, sh=alphaG1, rate=betaGF) +
        propFM[2] * pgamma(x, sh=alphaG1, rate=betaGM) 
}
curve(pgamma(x,sh=alphaG, rate=betaG), from=140, to=200, col=2, add=T, lwd=2)
curve(palt.g1(x), from=140, to=210, add=T, col=4, lwd=2)
legend("topleft", c("Gama 1","Gama 2"), lty=1, col=c(2,4), lwd=2)


## ----echo=TRUE,results="markup"------------------------------------------
## Ajustando distribuições
peso.n <- glm(peso ~ 1, family="gaussian", data=quest)
peso.g <- glm(peso ~ 1, family=Gamma(), data=quest)
peso.bc <- glm(-((1/peso)-1) ~ 1, family="gaussian", data=quest)
## Ajustando distribuições com médias diferentes para  masculino e feminino
peso.n1 <- glm(peso ~ sexo, family="gaussian", data=quest)
peso.g1 <- glm(peso ~ sexo, family=Gamma(), data=quest)
peso.bc1 <- glm(-((1/peso)-1) ~ sexo, family="gaussian", data=quest)
## Verossimilhanças dos ajustes
fits.peso <- rbind(c(logLik(peso.n), logLik(peso.g), 
                     (-1-1)*sum(log(peso)) + logLik(peso.bc), 2), 
                   c(logLik(peso.n1), logLik(peso.g1), 
                     (-1-1)*sum(log(peso)) + logLik(peso.bc1), 3))
dimnames(fits.peso) <- list(c("Modelo 1","Modelo2"), 
                            c("Normal","Gamma","Box-Cox","npar"))
fits.peso


## ----results="hide",echo=FALSE,fig.width=8,fig.height=6,out.width="0.98\\textwidth"----
par(mfrow=c(2,2), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
#names(summary(peso.n1))
#coef(peso.n1)
peso.n1.med <- cumsum(coef(peso.n1))
peso.n1.sd <- with(summary(peso.n1), sqrt(deviance/df.residual))
#altlm1 <- lm(altura ~sexo)
#names(summary(altlm1))
#summary(altlm1)$sigma
propFM <- prop.table(table(sexo))
plot(density(peso), main="peso", ylim=c(0, 0.03), xlab="")
curve(0.6*dnorm(x, m=peso.n1.med[1], sd=peso.n1.sd), from=35, to=150, col=4, add=T, lwd=1, lty=2)    
curve(0.4*dnorm(x, m=peso.n1.med[2], sd=peso.n1.sd), from=35, to=150, col=4, add=T, lwd=1, lty=2)    
dpeso.n1 <- function(x){
    propFM[1] * dnorm(x, m=peso.n1.med[1], sd=peso.n1.sd) +
        propFM[2] * dnorm(x, m=peso.n1.med[2], sd=peso.n1.sd)
}
curve(dnorm(x, mean=coef(peso.n), sd=with(summary(peso.n), sqrt(deviance/df.residual))), 
      from=35, to=150, col=2, add=T, lwd=2)
curve(dpeso.n1(x), from=35, to=150, add=T, col=4, lwd=2)
legend("topright", c("Normal 1","Normal 2"), lty=1, col=c(2,4), lwd=2)

plot(ecdf(peso))
ppeso.n1 <- function(x){
    propFM[1] * pnorm(x, m=peso.n1.med[1], sd=peso.n1.sd) +
        propFM[2] * pnorm(x, m=peso.n1.med[2], sd=peso.n1.sd)
}
curve(pnorm(x, mean=coef(peso.n), sd=with(summary(peso.n), sqrt(deviance/df.residual))), 
      from=35, to=150, col=2, add=T, lwd=2)
curve(ppeso.n1(x), from=35, to=150, add=T, col=4, lwd=2)
legend("topleft", c("Normal 1","Normal 2"), lty=1, col=c(2,4), lwd=2)
#@ 
#\end{frame}
#
#
#\begin{frame}[fragile]
#  \frametitle{Modelos para \code{altura}}
#
#  \red{Comparando ajustes para o modelo com distribuição Gama}
#<<results="hide",echo=FALSE,fig.width=8,fig.height=4,out.width="0.9\\textwidth"#>>=
##
## Distribuição Y ~Gama(shape = alpha, rate = beta)
## Ajuste Gama GLM (link = "inverse")
## Y ~Gama (\mu, \phi)
## Para y ~1:
## \mu = 1/beta_0 
## alpha = 1/phi = 1/summary(fit)$dispersion
## beta = beta_0/phi = coef(fit)/summary(fit)$dispersion
## Modelo 1 ~1
alphaG <- 1/summary(peso.g)$dispersion
betaG <- coef(peso.g) * alphaG
## Modelo 2 ~ sexo
alphaG1 <- 1/summary(peso.g1)$dispersion
betaGF <- coef(peso.g1)[1] * alphaG1
betaGM <- sum(coef(peso.g1)) * alphaG1

propFM <- prop.table(table(sexo))
#par(mfrow=c(1,2), mar=c(3,3, 1, 0), mgp=c(1.8, 0.8, 0))
plot(density(peso), main="peso", ylim=c(0, 0.03), xlab="")
curve(0.6*dgamma(x, sh=alphaG1, rate=betaGF), from=35, to=150, col=4, add=T, lwd=1, lty=2)    
curve(0.4*dgamma(x, sh=alphaG1, rate=betaGM), from=35, to=150, col=4, add=T, lwd=1, lty=2)    
dpeso.g1 <- function(x){
    propFM[1] * dgamma(x, sh=alphaG1, rate=betaGF) +
        propFM[2] * dgamma(x, sh=alphaG1, rate=betaGM) 
}
curve(dgamma(x,sh=alphaG, rate=betaG), from=35, to=150, col=2, add=T, lwd=2)
curve(dpeso.g1(x), from=35, to=150, add=T, col=4, lwd=2)
legend("topright", c("Gama 1","Gama 2"), lty=1, col=c(2,4), lwd=2)
##
plot(ecdf(peso))
ppeso.g1 <- function(x){
    propFM[1] * pgamma(x, sh=alphaG1, rate=betaGF) +
        propFM[2] * pgamma(x, sh=alphaG1, rate=betaGM) 
}
curve(pgamma(x,sh=alphaG, rate=betaG), from=35, to=150, col=2, add=T, lwd=2)
curve(ppeso.g1(x), from=35, to=140, add=T, col=4, lwd=2)
legend("topleft", c("Gama 1","Gama 2"), lty=1, col=c(2,4), lwd=2)

