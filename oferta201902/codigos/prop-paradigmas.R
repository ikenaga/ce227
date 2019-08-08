##
## Comandos utilizados para gerar slides comparando paradigmas para inferência
## Utilizando a inferẽncia sobre uma proporção como exemplo 
## - frequentista (distribuição amostral)
## - verossimilhança (função de verossimilhança) 
## - bayesiano (distribuição à posteriori)
##
## Responsável: Paulo Justiniano Ribeiro Junior


##
## Inferência frequentista - baseada na distribuição amostral
##

## ----echo=F--------------------------------------------------------------
rm(list=ls())
set.seed(2018)
th <- 0.17
POP <- sample(c(0,1), 10000000, prob=c((1-th), th), rep=TRUE)
AMs <- matrix(sample(POP, 80*10000, rep=T), nrow=80)
ps <- apply(AMs, 2, mean)

## ----echo=TRUE,results="markup"------------------------------------------
summary(ps)

## ----echo=FALSE,fig.width=6, fig.height=4, out.width = "0.75\\textwidth"----
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0))
hist(ps, prob=TRUE, main="", xlab=expression(hat(theta)), ylab="densidade")
lines(density(ps, bw=0.01), lwd=2)
#abline(v=quantile(ps, prob=c(0.025, 0.975)), lty=2)


## ------------------------------------------------------------------------
LWD=2

## ----amostralestimada----------------------------------------------------
par(mgp=c(2,0.8,0))
p <- 19/80
vp <- p*(1-p)/80
curve(dnorm(x, m=p, sd=sqrt(vp)), from=0.05, to=0.40, 
      ylab=expression(f(hat(theta))), xlab=expression(hat(theta)), ylim=c(0,10), lwd=LWD)
abline(v=p, lty=2, lwd=LWD)
abline(v=prop.test(19, 80)$conf, lty=2, lwd=LWD)


## ------------------------------------------------------------------------
LWD=1
par(mgp=c(2,0.8,0))
p <- 19/80
vp <- p*(1-p)/80
curve(dnorm(x, m=p, sd=sqrt(vp)), from=0.05, to=0.40, 
      ylab=expression(f(hat(theta))), xlab=expression(hat(theta)), ylim=c(0,10), lwd=LWD)
abline(v=p, lty=2, lwd=LWD)
abline(v=prop.test(19, 80)$conf, lty=2, lwd=LWD)
p0 <- 0.20
vp0 <- p0*(1-p0)/80
curve(dnorm(x, m=p0, sd=sqrt(vp0)), from=0.05, to=0.40, 
      ylab=expression(f(hat(theta))), xlab=expression(hat(theta)), add=TRUE, col=2, lwd=2)
segments(p, 0, p, dnorm(p, m=p0, sd=sqrt(vp0)), col=2)
pseq <- seq(p, 0.4, l=100)
polygon(c(pseq, rev(pseq)), c(rep(0, 100), rev(dnorm(pseq, m=p0, sd=sqrt(vp0)))), density=5, col=2)
legend("topright", c("estimada","sob H_0"), col=c(1,2), lwd=c(1, 2))
text(0.27, 1, "p-valor", col=2)


## ----echo=TRUE,results="markup"------------------------------------------
prop.test(19, 80)$conf
prop.test(19, 80, p=0.20, alt="greater")


## ------------------------------------------------------------------------
## Simulando sob H_0
th0 <- 0.20    ## supondo theta da hipótese
##POP0 <- sample(c(0,1), 10000000, prob=c((1-th0), th0), rep=TRUE)
##mean(POP0)   ## só para conferir...
##AM0s <- matrix(sample(POP0, 80*10000, rep=T), nrow=80)
set.seed(2018)
AM0s <- matrix(rbinom(80*10000, size=1, prob=0.20), nrow=80)
#dim(AMs)
p0s <- apply(AM0s, 2, mean)

## ----echo=TRUE,results="markup"------------------------------------------
summary(p0s)
(pvalor <- mean(p0s >= 19/80))

## ----fig.width=6, fig.height=4, out.width = "0.75\\textwidth"------------
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0))
hist(p0s, prob=TRUE, main="")
lines(density(p0s, bw=0.01), lwd=2)
abline(v=19/80, lty=2, col=2, lwd=2)
text(0.25, 8, substitute(pvalor==p, list(p=pvalor)), pos=4)
curve(dnorm(x, m=th0, sd=sqrt(th0*(1-th0)/80)), from=0, to=0.40, col=2, add=TRUE)
legend("topright", c("teórica","simulada"), lty=1, col=2:1)


##
## Inferência baseada na função de verossimilhança
##

## ----fig.width=5, fig.height=4, out.width = "0.7\\textwidth"-------------
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0))
Mf <- Vectorize(function(theta, y, size, log=F) dbinom(y, prob=theta, size, log=log))
th0=0.20
n=80
y=19
est <- 19/80
LMAX <- dbinom(y, size=n, prob=est)
th <- seq(0.08, 0.45, l=101)
sX17 <- drop(outer(th, y, Mf, size=n))
plot(th, sX17/LMAX, type="l", xlab=expression(theta), ylab=expression(paste("L(", theta, "|", "y=17)")), 
     xlim=range(th), ylim=c(-0.05, 1), lwd=2)
arrows(est, 1, est, 0, length=0.1, col=2)
LINT <- 0.25*dbinom(y, size=n, prob=est)/LMAX
abline(h=LINT, lty=3)
LL <- function(x){dbinom(y, size=n, prob=x)/LMAX - LINT}
int <- rootSolve:::uniroot.all(LL, c(0,1))
arrows(int, rep(LINT, 2), int, 0, length=0.1, lty=1, col=4)
LHIP <- dbinom(y, size=n, prob=th0)/LMAX
segments(th0, 0, th0, LHIP, col="darkolivegreen") 
arrows(th0, LHIP, min(th), LHIP, length=0.1, col="darkolivegreen") 
#text(0.5, 0, "0,5", cex=0.7)
text(th0, 0, expression(theta[0]), cex=0.7, col="darkolivegreen", pos=1, offset=0.2)
text(est, 0, expression(hat(theta)), pos=1, cex=0.7, offset=0.2, col=2)
text(int, 0, c(expression(hat(theta)[I]),expression(hat(theta)[S])),
             pos=1, cex=0.7, offset=0.2, col=4)
text(0.1, LHIP, expression(L(theta[0])/L(hat(theta))), pos=3, cex=0.8, col="darkolivegreen", offset=0.2)


## ----fig.width=10, fig.height=5, out.width = "0.99\\textwidth"-----------
par(mfrow=c(1,2), mar=c(2.9,2.9,0.2,1), mgp=c(1.7, 0.7, 0))
plot(th, -2*log(sX17/LMAX), type="l", xlab=expression(theta), ylab=expression(paste("D(", theta, "|", "y=17)")), 
     xlim=range(th), lwd=2)
abline(h=qchisq(c(0.9, 0.95, 0.99), df=1), lty=2)
plot(th, sqrt(-2*log(sX17/LMAX)), type="l", xlab=expression(theta), ylab=expression(sqrt(paste("D(", theta, "|", "y=17)"))), 
     xlim=range(th), lwd=2)
abline(h=sqrt(qchisq(c(0.9, 0.95, 0.99), df=1)), lty=2)


## ----fig.width=5, fig.height=4, out.width = "0.7\\textwidth"-------------
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0))
Mf <- Vectorize(function(theta, y, size, log=F) dbinom(y, prob=theta, size, log=log))
th0=0.20
n=80
y=19
est <- 19/80
LMAX <- dbinom(y, size=n, prob=est)
th <- seq(0.08, 0.45, l=101)
sX17 <- drop(outer(th, y, Mf, size=n))
plot(th, sX17/LMAX, type="l", xlab=expression(theta), ylab=expression(paste("L(", theta, "|", "y=17)")), 
     xlim=range(th), ylim=c(-0.05, 1), lwd=2)
arrows(est, 1, est, 0, length=0.1, col=2)
LINT <- 0.25*dbinom(y, size=n, prob=est)/LMAX
abline(h=LINT, lty=3)
LL <- function(x){dbinom(y, size=n, prob=x)/LMAX - LINT}
int <- rootSolve:::uniroot.all(LL, c(0,1))
arrows(int, rep(LINT, 2), int, 0, length=0.1, lty=1, col=4)
LHIP <- dbinom(y, size=n, prob=th0)/LMAX
segments(th0, 0, th0, LHIP, col="darkolivegreen") 
arrows(th0, LHIP, min(th), LHIP, length=0.1, col="darkolivegreen") 
#text(0.5, 0, "0,5", cex=0.7)
text(th0, 0, expression(theta[0]), cex=0.7, col="darkolivegreen", pos=1, offset=0.2)
text(est, 0, expression(hat(theta)), pos=1, cex=0.7, offset=0.2, col=2)
text(int, 0, c(expression(hat(theta)[I]),expression(hat(theta)[S])),
             pos=1, cex=0.7, offset=0.2, col=4)
text(0.1, LHIP, expression(L(theta[0])/L(hat(theta))), pos=3, cex=0.8, col="darkolivegreen", offset=0.2)


## ----fig.width=10, fig.height=5, out.width = "0.99\\textwidth"-----------
par(mfrow=c(1,2), mar=c(2.9,2.9,0.2,1), mgp=c(1.7, 0.7, 0))
plot(th, -2*log(sX17/LMAX), type="l", xlab=expression(theta), ylab=expression(paste("D(", theta, "|", "y=17)")), 
     xlim=range(th), lwd=2)
abline(h=qchisq(c(0.9, 0.95, 0.99), df=1), lty=2)
plot(th, sqrt(-2*log(sX17/LMAX)), type="l", xlab=expression(theta), ylab=expression(sqrt(paste("D(", theta, "|", "y=17)"))), 
     xlim=range(th), lwd=2)
abline(h=sqrt(qchisq(c(0.9, 0.95, 0.99), df=1)), lty=2)


##
## Inferência bayesiana - baseada na distribuição a posteriori
##

## Obtenção de prioris e posterioris a partir de diferentes opiniões
## ----echo=FALSE----------------------------------------------------------
source("bayes-fun-01.R")

## ----priori03,echo=FALSE-------------------------------------------------
prI <- prioriBeta(0.4, c(0.30, 0.50), 0.70)
postI <- postBinom(19, 80, prI, plot=FALSE)


## ----echo=FALSE,results = "hide", fig.width=6, fig.height=4, out.width = "0.75\\textwidth"----
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0))
postBinom(19, 80, prI)


## ----echo=TRUE-----------------------------------------------------------
(prI <- prioriBeta(0.4, c(0.30, 0.50), 0.70)) 
postBinom(19, 80, prI, plot=FALSE)


## ----echo=FALSE----------------------------------------------------------
prII <- prioriBeta(0.08, c(0.03, 0.20), 0.90)
postII <- postBinom(19, 80, prII, plot=FALSE)


## ----results="hide",fig.width=6, fig.height=4, out.width = "0.75\\textwidth"----
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0))
prII <- prioriBeta(0.08, c(0.03, 0.20), 0.90)
postII <- postBinom(19, 80, prII)


## ----echo=TRUE-----------------------------------------------------------
(prII <- prioriBeta(0.08, c(0.03, 0.20), 0.90))
postBinom(19, 80, prII, plot=FALSE)


## ----echo=FALSE----------------------------------------------------------
prIII <- prioriBeta(0.5, c(0.05, 0.95), 0.90)
postIII <- postBinom(19, 80, prIII, plot=FALSE)


## ----results="hide",fig.width=6, fig.height=4, out.width = "0.75\\textwidth"----
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0))
postIII <- postBinom(19, 80, prIII)


## ----echo=TRUE-----------------------------------------------------------
(prIII <- prioriBeta(0.50, c(0.05, 0.95), 0.90))
postBinom(19, 80, prIII, plot=FALSE)


## ----results="hide",fig.width=6, fig.height=4, out.width = "0.99\\textwidth"----
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0), mfrow=c(2,3))
res <- postBinom(19, 80, c(3,3), cex.leg=0.6, ylim=c(0,10))
res <- postBinom(19, 80, c(1,1), cex.leg=0.6, ylim=c(0,10))
res <- postBinom(19, 80, c(0.5,5), cex.leg=0.6, ylim=c(0,10))
res <- postBinom(19, 80, c(10,7), cex.leg=0.6, ylim=c(0,10))
res <- postBinom(19, 80, c(16,11), cex.leg=0.6, ylim=c(0,10))
res <- postBinom(19, 80, c(43,29), cex.leg=0.6, ylim=c(0,10))


## ----fig.width=6, fig.height=4, out.width = "0.99\\textwidth"------------
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.7, 0), mfrow=c(2,3))
res <- postBinom(6, 27, c(16,11), cex.leg=0.6, ylim=c(0,25))
text(0.8, 18, "n=27 e y=7")
res <- postBinom(9, 40, c(16,11), cex.leg=0.6, ylim=c(0,25))
text(0.8, 18, "n=40 e y=9")
res <- postBinom(19, 80, c(16,11), cex.leg=0.6, ylim=c(0,25))
text(0.8, 18, "n=80 e y=19")
res <- postBinom(38, 160, c(16,11), cex.leg=0.6, ylim=c(0,25))
text(0.8, 18, "n=160 e y=38")
res <- postBinom(95, 400, c(16,11), cex.leg=0.6, ylim=c(0,25))
text(0.8, 18, "n=400 e y=95")
res <- postBinom(190, 800, c(16,11), cex.leg=0.6, ylim=c(0,25))
text(0.8, 18, "n=800 e y=190")

