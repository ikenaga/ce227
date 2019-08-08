## Em detalhes: informação à priori-------------------------------------------------
## Y \sim B(alpha, beta)
## E[Y] = \frac{alpha}{alpha+beta}
## Mo[Y] = \frac{alpha-1}{alpha+beta-2} ; para alpha, beta > 1
## Var[Y] = \frac{alpha beta}{(alpha+beta)^2 (alpha+beta+1)}

priori.beta <- function(est, valor, prob, tipo = c("media","moda"), curve=FALSE, ...){
    ## se tipo = "media" : valor é a ME, prob é a "confiança" de um intervalo (simétrico, assintótico)
    ## se tipo = "moda"  : valor é um limite inferior de um intervalo e prob a a prob abix deste limite
    tipo <- match.arg(tipo, c("media","moda"))
    if(tipo == "moda"){
        moda <- est
        beta.f <- function(beta){
            alpha <- ((2*moda - 1 - beta*moda)/(moda - 1))
            pbeta(valor, alpha, beta) - prob
        }
        beta.est <- uniroot(beta.f, c(1, 100))$root
        alpha.est <- ((2*moda - 1 - beta.est*moda)/(moda - 1))
    }
    else{
        beta.f <- function(beta){
            alpha <- (est*beta)/(1-est)
            z <- qnorm(1-(1-prob)/2)
            v <- (alpha*beta)/(((alpha+beta)^2)*(alpha+beta+1))
            return(valor - z * sqrt(v))
        }
        beta.est <- uniroot(beta.f, lower=1, upper=100)$root
        alpha.est <- (est*beta.est)/(1-est)
    }
    if(curve){
        distbeta <- function(x){dbeta(x, alpha.est, beta.est)}
        curve(distbeta(x), 
              ylab=expression(p(theta)), xlab=expression(theta),
              ... )
        abline(v=est)
    }
    return(c(alpha=alpha.est, beta=beta.est))
}

par(mfrow=c(1,2))
priori.beta(est=.35, valor=0.20, prob=0.05, curve=TRUE, tipo="moda")

priori.beta(est=.35, valor=0.15, prob=0.95, curve=TRUE, tipo="media")


pri[1]/(pri[1]+pri[2])
(pri[1]-1)/(pri[1]+pri[2]-2)

pbeta(0.20, pri[1], pri[2])

curve(dbeta(x, pri[1], pri[2]), from=0, to=1) 


priori.beta <- function(est, low, conf){
  moda <- est
  fd <- function(beta, conf, ic, moda){
    diff(pbeta(c(ic[1],ic[2]), ((2*moda - 1 - beta*moda)/(moda - 1)), beta))}
  op <- optim(1.00001, fd, method="BFGS", conf=0.9, ic=c(0.4, 0.7), moda=0.4, 
              control=list(fnscale=-1))
  beta <- op$par
  alpha <- ((2*moda - 1 - beta*moda)/(moda - 1))
  distbeta <- function(x){dbeta(x, alpha, beta)}
  curve(distbeta(x), ylim=c(0,10), ylab=expression(beta(theta)), xlab=expression(theta))
  return(list(alpha=alpha, beta=beta))
}

pri <- priori.beta(est=.35, ic=c(.28,.42), conf=.85); pri


#-----------------------------------------------------------------------------------

## Coletei os dados-----------------------------------------------------------------

# y = 24
# n = 80

## Gráfico dos dados coletados
#curve(dbinom(24,80,x),0,1)

## Gráfico dos dados coletados (escalonado)
f.bin <- function(x){dbeta(x, 25, 57)}
curve(f.bin, col=2, lty=2, add=T)
#-----------------------------------------------------------------------------------

## Posteriori-----------------------------------------------------------------------

post.beta <- function(post){
  post <- unname(unlist(post))
  beta <- post[2] - post[1] + post[4]
  alpha <- post[1] + post[3]
  distbeta <- function(x){dbeta(x, alpha, beta)}
  curve(distbeta(x), lty=3, col=4, add=T)
  legend("topright", legend=c(expression(paste("[",theta, "]")), 
                              expression(paste("[",y, "|", theta, "]")),
                              expression(paste("[",theta, "|", x, "]"))), 
         col=c(1, 2, 4), bty="n", lty=1:3)
  return(list(alpha=alpha, beta=beta))
}

post <- c(24, 80, pri[1], pri[2])

post.beta(post)
#-----------------------------------------------------------------------------------




## Direto, com lattice -------------------------------------------------------------

library(lattice)

post.beta <- function(est, ic, conf, n, y){
  vec <- seq(0.0001, 0.9999, l=200)
  moda <- est
  # Cálculo da priori
  fd <- function(b, conf, ic, moda){
    diff(pbeta(c(ic[1],ic[2]), ((2*moda - 1 - b*moda)/(moda - 1)), b))}
  op <- optim(1.00001, fd, method="BFGS", conf=0.9, ic=c(0.4, 0.7), moda=0.4, 
              control=list(fnscale=-1))
  b.pri <- op$par
  a.pri <- ((2*moda - 1 - b.pri*moda)/(moda - 1))
  beta.pri <- function(x, a.pri, b.pri){dbeta(x, a.pri, b.pri)}
  priori <- beta.pri(vec, a.pri, b.pri)
  # coleta de dados
  f.bin <- function(x){dbeta(x, (y+1), (n-y+1))}
  dados <- f.bin(vec)
  # Cálculo da posteriori
  b.pos <- n - y + b.pri
  a.pos <- y + a.pri
  beta.pos <- function(x, a.pos, b.pos){dbeta(x, a.pos, b.pos)}
  posteriori <- beta.pos(vec, a.pos, b.pos)
  resul <- data.frame(Priori=priori, Dados=dados, Posteriori=posteriori, Parametro=vec)
  graf <- xyplot(Priori + Dados + Posteriori ~ Parametro , data=resul, type="l", auto.key=T, 
                 ylab=expression(beta(theta)), xlab=expression(theta))
  return(list(alpha.pri=a.pri, beta.pri=b.pri, alpha.pos=a.pos, beta.pos=b.pos, graf))
}

post.beta(est=.5, ic=c(.1,.7), conf=.75, n=80, y=24)
#-----------------------------------------------------------------------------------




## Entrando com dados---------------------------------------------------------------

alpha <- 25
beta <- 57
conf <- .95
#-----------------------------------------------------------------------------------




## Cálculo da média ----------------------------------------------------------------

media.beta <- function(alpha, beta){
  esp <- function(x){x * dbeta(x, alpha, beta)}
  valor <- integrate(esp, 0,1)$value
  return(valor)
}

media.beta(alpha, beta)

media.analit <- alpha/(alpha + beta)
#-----------------------------------------------------------------------------------



## Cálculo da mediana --------------------------------------------------------------

mediana.beta <- function(a, b, chute){
  fun <- function(med, a, b){
    (pbeta(med, a, b) - pbeta(0, a, b) - 0.5)^2
  }
  med <- optim(chute, fun, a=alpha, b=beta, 
               control=list(fnscale=1), method="BFGS")$par
  #curve(fun(x, alpha, beta))
  return(med)
}

mediana.beta(chute=0.2, alpha, beta)

mediana.analit <- (alpha - 1/3)/(alpha + beta - 2/3)
#-----------------------------------------------------------------------------------



## Cálculo da moda -----------------------------------------------------------------

moda.beta <- function(a, b, chute){
  fun <- function(x, a, b){dbeta(x, a, b)}
  moda <- optim(chute, fun, a=alpha, b=beta, 
                control=list(fnscale=-1), method="BFGS")$par
  return(moda)
}

moda.beta(alpha, beta, chute=0.2)
moda.analitica <- (alpha - 1)/(alpha + beta - 2)
#-----------------------------------------------------------------------------------



## Intervalo HPD--------------------------------------------------------------------


HPD.beta <- function(conf=0.95, tol=0.00000001, alpha, beta){
  # Significância
  conf <- min(conf, 1-conf) 
  # função que será itimizada em relação a x
  f <- function(x, conf, a, b){
    diff(qbeta(c(x, 1-conf+x), a, b))
  }
  # otimização: minimizará o valor de x que irá manter f válida
  val <- optimize(f, c(0, conf), a=alpha, b=beta, conf=conf, tol=tol)
  #obtendo os valores do intervalo HPD
  ic <- c(qbeta(val$minimum, alpha, beta), qbeta(1-conf+val$minimum, alpha, beta))
  
  f.int <- dbeta(ic, alpha, beta)
  #unname(diff(pbeta(int, alpha, beta)))
  curve(dbeta(x, alpha, beta), ylab=expression(beta(theta)), xlab=expression(theta))
  lines(c(ic), c(f.int), col=2)
  arrows(ic, f.int, ic, 0, length=0.1, col=2)
  return(c(Int.I=ic[1], Int.S=ic[2]))
}

HPD.beta(alpha=alpha, beta=beta)

text(0.3, 4, "HPD")
text(0.31, 3, "95%")



## HPD mais genérico
HPD <- function(distpost, conf=0.95, tol=0.00000001,...){
  conf <- min(conf, 1-conf)
  f <- function(x,distpost,conf,...){
    diff(distpost(c(x, 1-conf+x),...))
  }
  out <- optimize(f, c(0,conf), distpost = distpost,
                  conf=conf, tol=tol, ...)
  return(c(distpost(out$minimum,...), distpost(1-conf+out$minimum,...)))
}

HPD(qgamma, shape=25, rate=1)

#-----------------------------------------------------------------------------------






#### Rascunho HPD....


## Primeira tentativa....
library(rootSolve)

HPD <- function(x, aa, bb, cc){
  c(alt = (diff(dbeta(c(x[1],x[2]), aa, bb))^2),
    larg = ((diff(pbeta(c(x[1],x[2]), aa, bb, lower.tail=F)) - cc)^2))
}

resp <- multiroot(HPD, start=c(0.21,.40), aa=alpha, bb=beta, cc=conf)

# Conferindo
dbeta(resp$root, alpha, beta)


# Altura do intervalo
dbeta(S, a, b) - dbeta(I, a, b)

# Largura do intervalo
pbeta(S, a, b) - pbeta(I, a, b) - conf



## Segunda tentativa....

a <- 25
b <- 57
c <- .95
m <- moda.beta(a, b, chute=0.2)

HPD <- function(a, b, c, m){
  f.mo <- dbeta(m, a, b)
  #ic <- qbeta(c(0.025, 0.975), a, b)
  #val <- unname(c(ic, f.mo-.1))
  
  int <- function(val, a, b, c, f.mo){
    (diff(pbeta(c(val[1],val[2]), a, b)) - c)^2
    f.mo - diff(dbeta(c(val[1],val[2]), a, b)) - val[3]
  }
  
  op <- optim(c(0.21, 0.4, 7), int, a=alpha, b=beta, c=conf, f.mo=f.mo, 
              control=list(fnscale=-1), method="BFGS") # lower=, upper=,
  op$par
  return(op$par)
}


## Terceira tentativa....

HPD <- function(a, b, c){
  
  int <- function(val, a, b, c){
    (diff(pbeta(c(val[1],val[1]+val[2]), a, b)) - c)^2
  }
  
  op <- optim(c(0.1, 0.3), int, a=alpha, b=beta, c=conf,
              control=list(fnscale=-1), method="BFGS") # lower=, upper=,
  # Conferindo
  op$par
  dbeta(c(op$par[1], op$par[1]+op$par[2]), alpha, beta)
  val <- c(op$par[1], op$par[1]+op$par[2])
  
  return(op$par)
}


val <- c(op$par[1], op$par[1]+op$par[2])





## 26/03----------------------------------------------------------------------------


## 1) Priori uniforme U(0,1)
unif <- function(x, a=a, b=b) (x^{0})*1/(b-a)

## 2) Priori beta(alpha, beta)
beta <- function(x, a=A, b=B){
  gamma(a+b)/(gamma(a)*gamma(b))*(x^{a-1})*((1-x)^{b-1})
}

## 3) Priori Jeffreys (integrando 1)
jef <- function(x){
  f <- function(x){
    sqrt(1/(x*(1-x)))
  }
  c <- integrate(f, 0, 1)$value
  sqrt(1/(x*(1-x)))/c
}

## Função de verossimilhança da dist. binomial
L <- function(theta, n=N, y=Y){
  dbinom(y, size=n, prob=theta)
}

# Deixando a fç L na escala do parâmetro theta
plot.L <- function(x, y=Y, n=N){
  dbeta(x, shape1=y+1, shape2=n-y+1)
}

## Posteriori 1)
post1 <- function(x, n=N, y=Y){
  dbeta(x, shape1=(y+1), shape2=(n-y+1))
}

## Posteriori 2)
post2 <- function(x, n=N, y=Y, a=A, b=A){
  dbeta(x, shape1=(y+a), shape2=(n+b-y))
}

## Posteriori 3)
post3 <- function(x, n=N, y=Y){
  dbeta(x, shape1=(y+0.5), shape2=(n-y+0.5))
}


## Valores iniciais

# Fç vero
N <- 15
Y <- 6
# priori unif
a <- 0
b <- 1
# priori beta
A <- 2.625
B <- 2.625


par(mfrow=c(1,3))

## Cenário 1
curve(post1, col=2, ylim=c(0, 4), lwd=3,
      xlab=expression(theta), ylab="y")
curve(unif(x, a=a, b=b), 0,1, col=4, add=T)
curve(plot.L(x), col=1, add=T)
legend("topright", 
       legend=c(expression(paste("[",theta, "]")), 
                expression(paste("[",y, "|", theta, "]")),
                expression(paste("[",theta, "|", x, "]"))), 
       col=c(4,1,2), bty="n", lty=1)

## Cenário 2
curve(post2, col=2, ylim=c(0, 4),
      xlab=expression(theta), ylab="y")
curve(beta(x), col=4, add=T)
curve(plot.L(x), col=1, add=T)
legend("topright", 
       legend=c(expression(paste("[",theta, "]")), 
                expression(paste("[",y, "|", theta, "]")),
                expression(paste("[",theta, "|", x, "]"))), 
       col=c(4,1,2), bty="n", lty=1)

## Cenário 3
curve(post3, col=2, ylim=c(0,4),
      xlab=expression(theta), ylab="y")
curve(jef(x), col=4, add=T)
curve(plot.L(x), col=1, add=T)
legend("topright", 
       legend=c(expression(paste("[",theta, "]")), 
                expression(paste("[",y, "|", theta, "]")),
                expression(paste("[",theta, "|", x, "]"))), 
       col=c(4,1,2), bty="n", lty=1)

layout(1)
##----------------------------------------------------------------------------------




## 02/04 - amostrador de Gibbs -----------------------------------------------------

## simulando dados
set.seed(15)
x <- sort(runif(20, 0, 20))
y <- 2 + 0.25*x + rnorm(20, m=0, sd=2)
## visualizando dados e reta de regressão "usual"
plot(x,y)
reg <- lm(y ~ x)
abline(reg)

## Código ("ingênuo") para amostrador de Gibbs
reg.lin.GIBBS <- function(y, x, Nsim, iniBeta){ 
  library(MASS)
  n <- length(y) # comprimento de y
  X <- cbind(1, x) # matriz do modelo de regressão
  XX <- crossprod(X) # X'X
  bhat <- solve(XX,crossprod(X,y)) # Estimativa do parâmetro beta
  ## informações para o laço
  sim <- matrix(0, nrow = Nsim, ncol=3) # Objeto nulo, receberá valores simulados
  sim[1,1:2] <- iniBeta # primeira linha de sim com o chute de beta
  sim[1, 3]  <- 1/rgamma(1, shape=n/2, scale=2/crossprod(y-X%*%sim[1,1:2]))
  ## laço
  for(i in 2:Nsim){
    sim[i, 3]  <- 1/rgamma(1, shape=n/2, scale=2/crossprod(y-X%*%sim[(i-1),1:2]))
    sim[i,1:2] <- mvrnorm(n=1, mu=bhat, Sigma = solve(XX)*sim[i, 3])
  }
  return(sim)
}
## Obtendo amostras
rlG0 <- reg.lin.GIBBS(y=y, x=x, Nsim=15000, iniBeta=c(0,0.6))
## Burn-in (3000) e thining (1/10)
rlG1 <- rlG0[-(1:3000),] # aquecimento da simulação
rlG2 <- rlG1[10*(1:1200),] # file mignon da simulação (raleamento)

## comparando estimativas "usuais" e médias das posterioris
c(coef(reg), summary(reg)$sigma^2) 
apply(rlG2, 2, mean)
# comparação de variâncias não é linear, ex: 1.678^2 e 2.378^2


## traços e posterioris (por simulação) com indicação das estimativas "usuais"
x11()
par(mfrow=c(2,3))
plot(rlG2[,1], type="l", xlab="b0", ylab="y")
plot(rlG2[,2], type="l", xlab="b1", ylab="y")
plot(rlG2[,3], type="l", xlab=expression(sigma^2), ylab="y")
plot(density(rlG2[,1]), xlab="b0", ylab="y", main=""); abline(v=coef(reg)[1])
plot(density(rlG2[,2]), xlab="b1", ylab="y", main=""); abline(v=coef(reg)[2])
plot(density(rlG2[,3]),xlab=expression(sigma^2), ylab="y", main=""); abline(v=summary(reg)$sigma^2)
layout(1)
##----------------------------------------------------------------------------------




## 04/04 - código agora incluindo predições ----------------------------------------

reg.lin.GIBBS <- function(y, x, Nsim, iniBeta, x.pred = NULL){  
  require(MASS)
  n <- length(y)
  X <- cbind(1, x)
  XX <- crossprod(X)
  bhat <- solve(XX,crossprod(X,y))
  ## preparando objetos para sair a coluna da produção em sim
  NC <- ncol(X) + 1
  if(!is.null(x.pred)){ 
    Xpred <- cbind(1, x.pred)
    NC <- NC + nrow(Xpred) 
  }
  sim <- matrix(0, nrow = Nsim, ncol=NC)
  sim[1,1:2] <- iniBeta
  sim[1, 3]  <- 1/rgamma(1, shape=n/2, scale=2/crossprod(y-X%*%sim[1,1:2]))
  ##
  for(i in 2:Nsim){
    sim[i, 3]  <- 1/rgamma(1, shape=n/2, scale=2/crossprod(y-X%*%sim[(i-1),1:2]))
    sim[i,1:2] <- mvrnorm(n=1, mu=bhat, Sigma = solve(XX)*sim[i, 3])
    if(!is.null(x.pred)) # simulo y de uma normal com u=b0hat+b1hat*x e sigma=(X'X)^-1*sigmahat, x que quero predição
      sim[i,-(1:3)] <- rnorm(nrow(Xpred), m=Xpred%*%sim[i,1:2], sd=sqrt(sim[i,3]))
  }
  return(sim)
}
## agora com predição
rlG0 <- reg.lin.GIBBS(y=y, x=x, Nsim=15000, iniBeta=c(0,0.6), x.pred=1:20)
dim(rlG0)
rlG1 <- rlG0[-(1:3000),]
dim(rlG1)
rlG2 <- rlG1[10*(1:1200),]
dim(rlG2)
rbind(dim(rlG0), dim(rlG1), dim(rlG2))

reg
summary(reg)$sigma^2

x11()
par(mfrow=c(4,5))
apply(rlG2[,-(1:3)], 2, function(x) plot(density(x),type="l", main=""))
par(mfrow=c(1,1))

x11()
#usando rlG2
y.pred <- apply(rlG2[,-(1:3)], 2, mean)
plot(y ~ x)
abline(reg)
lines(1:20, y.pred, col=2)
#usando rlG1
y.pred <- apply(rlG1[,-(1:3)], 2, mean)
plot(y ~ x)
abline(reg)
lines(1:20, y.pred, col=4)

##----------------------------------------------------------------------------------




## 04/04 Exemplos JAGS/rjags -------------------------------------------------------
# Download do Jags
# http://sourceforge.net/projects/mcmc-jags/files/latest/download

## Média e variância (precisão) da distribuição normal
n <- 20
x <- rnorm(n, 70, 5)

#write.table(x,
#            file = 'normal.data',
#            row.names = FALSE,
#            col.names = FALSE)

cat( "model {
     for (i in 1:n){
     x[i] ~ dnorm(mu, tau)
     }
     mu ~ dnorm(0, 0.0001)
     sigma ~ dunif(0, 100)     
     tau <- pow(sigma, -2)
     }", file="normal.modelo"
)


require(rjags)

## OBS: valores iniciais são dispensáveis neste exemplo
inis <- list(list(mu=10, sigma=2),
             list(mu=50, sigma=5),
             list(mu=70, sigma=10))

jags <- jags.model('normal.modelo',
                   data = list('x' = x, 'n' = n),
                   n.chains = 3,
                   inits = inis,
                   n.adapt = 100)

#update(jags, 1000)

#sam <- jags.samples(jags, c('mu', 'tau'), 1000)
sam <- coda.samples(jags, c('mu', 'tau'), n.iter=10000, thin=10)

par(mfrow=c(2,2))
plot(sam)
str(sam)
summary(sam)
HPDinterval(sam)



## regressão linear simples

cat( "model {
     for (i in 1:n){
     y[i] ~ dnorm(mu[i], tau)
     mu[i] <- b0 + b1 * x[i]
     }  		
     b0 ~ dnorm(0, .0001)
     b1 ~ dnorm(0, .0001)
     tau <- pow(sigma, -2)
     sigma ~ dunif(0, 100)
     }", file="reglin.modelo")
 
## alternativa:
## tau ~ dgamma(0.001, 0.001)
## sigma2 <- 1/tau

n <- 20
x <- sort(runif(n, 0, 20))
epsilon <- rnorm(n, 0, 2.5)
y <- 2 + 0.5*x + epsilon

#write.table(data.frame(X = x, Y = y, Epsilon = epsilon),
#            file = 'reglin.dados',
#            row.names = FALSE,
#            col.names = TRUE)

require(rjags)
# require(R2jags)
# require(dclone)

inis <- list(list(b0=0, b1=1, sigma=1),
             # list(b0=1, b1=0.5, sigma=2),
             list(b0=2, b1=0.1, sigma=5))

jags <- jags.model('reglin.modelo',
                   data = list('x' = x,'y' = y,'n' = n),
                   n.chains = 2,
                   # inits=inits,
                   n.adapt = 100)

#update(jags, 1000)
class(jags)

sam <- jags.samples(jags,c('b0', 'b1', 'sigma'),1000)
class(sam)

sam <- coda.samples(jags, c('b0', 'b1', 'sigma'),1000)
class(sam)
str(sam)
plot(sam)

sam <- coda.samples(jags, c('b0', 'b1', 'sigma', "y"),1000)
str(sam)
int <- HPDinterval(sam)
str(int)
## complementar com gráficos, resumos, inferências de interesse, etc

##----------------------------------------------------------------------------------
