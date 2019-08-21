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
source("../dados/mandible.R")
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
ll1 <- function(mu, theta, lsigma, am1, am2){
    sigma <- exp(lsigma)
    l1 <- sum(dnorm(am1, m=mu, sd=sigma, log=T))    
    l2 <- sum(dnorm(am2, m=mu+theta, sd=sigma, log=T))
    return(-(l1+l2))
}    


## ----echo=TRUE-----------------------------------------------------------
fit1 <- mle(ll1, start=list(mu=110, theta=0, lsigma=log(10)), 
           fixed=list(am1 = mandible$fem, am2 = mandible$male))


## ----echo=TRUE-----------------------------------------------------------
summary(fit1)
confint(fit1, level=0.95)
prof1 <- profile(fit1)

## ----results="hide",fig.width=7, fig.height=3, out.width = "0.99\\textwidth"----
par(mfrow=c(1,3), mar=c(2.5,2.5,0.3, 0.3), mgp=c(1.7, 0.7, 0))
plot(prof1, level=c(0.5, 0.8, 0.95, 0.99))


## variâncias diferentes
ll2 <- function(mu, theta, lsigma1, lsigma2, am1, am2){
    sigma1 <- exp(lsigma1)
    sigma2 <- exp(lsigma2)
    l1 <- sum(dnorm(am1, m=mu, sd=sigma1, log=T))    
    l2 <- sum(dnorm(am2, m=mu+theta, sd=sigma2, log=T))
    return(-(l1+l2))
}    
fit2 <- mle(ll2, start=list(mu=110, theta=0, lsigma1=log(10),
            lsigma2=log(10)), 
           fixed=list(am1 = mandible$fem, am2 = mandible$male))
summary(fit2)
confint(fit2, level=0.95)
prof2 <- profile(fit2)
logLik(fit2)
par(mfrow=c(1,4), mar=c(2.5,2.5,0.3, 0.3), mgp=c(1.7, 0.7, 0))
plot(prof2, level=c(0.5, 0.8, 0.95, 0.99))


## ----echo=TRUE-----------------------------------------------------------
fot0 <- lm(mandible~1, data=mand)
(lLs <- c(logLik(fit0), logLik(fit1), logLik(fit2)))
nps <- c(2, 3, 4)
(dlLs <-2* diff(lLs))
dnps <- diff(nps)
2*pchisq(dlLs, df=dnps, lower=FALSE)

## ----echo=TRUE-----------------------------------------------------------
-2 * logLik(fit0) + 2 * 2
-2 * logLik(fit1) + 2 * 3
-2 * logLik(fit2) + 2 * 4
c(AIC(fit0), AIC(fit1), AIC(fit2))
##
-2 * logLik(fit0) + log(20) * 2
-2 * logLik(fit1) + log(20) * 3
-2 * logLik(fit2) + log(20) * 4
c(AIC(fit0, k=log(20)), AIC(fit1, k=log(20)), AIC(fit2, k=log(20)))

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
