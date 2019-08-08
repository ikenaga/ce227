#-------------------------------------------------------------------------------
# distribuição beta
# função densidade: f(x, a, b) = gamma(a+b)/(gamma(a)*gamma(b))*x^(a-1)*(1-x)^(b-1)
# suporte: 0 < x < 1
# parâmetros:  shape1 > 0  e shape2 > 0 (a=shape1, b=shape2)
# E(X) = a/(a+b) e  V(X) = a*b/((a+b)^2*(a+b+1))

#beta.panel <- function(panel){
#  curve(dbeta(x, shape=panel$shape1, shape2=panel$shape2), 0, 1)
#  panel
#}
 
#panel <- rp.control()
#rp.slider(panel, shape1, 0.1, 20, initval=5, showvalue=TRUE, action=beta.panel)
#rp.slider(panel, shape2, 0.1, 20, initval=5, showvalue=TRUE, action=beta.panel)

##
## Função para definir os parâmetros de uma priori beta
## a partir de:
##    - uma estimativa pontual (considerada como moda) (p)
##    - um intervalo de valores (int) e a "confiança" a ele atribuída (conf)
##
##  Lembrete: a moda da Beta é \frac{\alpha - 1}{\alpha + \beta -2} , para \alpha e \beta>1
#prioriBeta <- function(p, int, conf){
#	betapars.f <- function(beta, p, int, conf){
#		alfa <- ((beta-2)*p + 1)/(1-p)
#		val <- (diff(pbeta(int, shape1=alfa, shape2=beta)) - conf)^2
#                # print(c(alfa, beta, val))
#		return(val)
#		}
#	res <- optimize(betapars.f, interval=c(1,50), p=p, conf=conf, int=int)
#	pars <- c(alpa = ((res$minimum-2)*p + 1)/(1-p), beta=res$minimum)
#        # attr(pars, "resNum") <- res
#        return(pars)
#}

#(pr01 <- prioriBeta(0.4, c(0.1, 0.7), 0.90))
#(pr02 <- prioriBeta(0.4, c(0.38, 0.42), 0.95))
#(pr03 <- prioriBeta(0.5, c(0.45, 0.55), 0.90))
#curve(dbeta(x, pr02[1], pr02[2]))
#curve(dbeta(x, pr01[1], pr01[2]), add=T)
#curve(dbeta(x, pr03[1], pr03[2]), add=T)

## conferindo:
#diff(pbeta(c(0.15, 0.35), pr01[1], pr01[2]))
#unname((pr01[1]-1)/(pr01[1] + pr01[2] - 2))


## Em detalhes: informação à priori-------------------------------------------------
## Y \sim B(alpha, beta)
## E[Y] = \frac{alpha}{alpha+beta}
## Mo[Y] = \frac{alpha-1}{alpha+beta-2} ; para alpha, beta > 1
## Var[Y] = \frac{alpha beta}{(alpha+beta)^2 (alpha+beta+1)}

prioriBeta <- function(p, int, conf){
    priori.beta(est=p, valor=diff(int)/2, prob=conf, tipo="media")
}

priori.beta <- function(est, valor, prob, tipo = c("media","moda"), curve=FALSE, ...){
    ## se tipo = "media" : valor é a ME, prob é a "confiança" de um intervalo (simétrico, assintótico)
    ## se tipo = "moda"  : valor é um limite inferior de um intervalo e prob a a prob abaixo deste limite
    tipo <- match.arg(tipo, c("media","moda"))
    if(tipo == "moda"){
        moda <- est
        beta.f <- function(logbeta){
            beta <- exp(logbeta)
            alpha <- ((2*moda - 1 - beta*moda)/(moda - 1))
            pbeta(valor, alpha, beta) - prob
        }
        beta.est <- exp(uniroot(beta.f, c(0, 50))$root)
        alpha.est <- ((2*moda - 1 - beta.est*moda)/(moda - 1))
    }
    else{
        beta.f <- function(logbeta){
            beta <- exp(logbeta)
            alpha <- (est*beta)/(1-est)
            z <- qnorm(1-(1-prob)/2)
            v <- (alpha*beta)/(((alpha+beta)^2)*(alpha+beta+1))
            return(valor - z * sqrt(v))
        }
        beta.est <- exp(uniroot(beta.f, lower=0, upper=50)$root)
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


#(pr01 <- prioriBeta(0.4, c(0.1, 0.7), 0.90))
#priori.beta(0.4, 0.3, 0.90)
#(pr02 <- prioriBeta(0.4, c(0.38, 0.42), 0.95))
#priori.beta(0.4, 0.02, 0.95)
#(pr03 <- prioriBeta(0.5, c(0.45, 0.55), 0.90))
#priori.beta(0.5, 0.05, 0.90)
#curve(dbeta(x, pr02[1], pr02[2]))
#curve(dbeta(x, pr01[1], pr01[2]), add=T)
#curve(dbeta(x, pr03[1], pr03[2]), add=T)

## conferindo:
#diff(pbeta(c(0.1, 0.7), pr01[1], pr01[2]))
#diff(pbeta(c(0.38, 0.42), pr02[1], pr02[2]))
#diff(pbeta(c(0.45, 0.55), pr03[1], pr03[2]))

#par(mfrow=c(1,2))
#priori.beta(est=.35, valor=0.20, prob=0.05, curve=TRUE, tipo="moda")
#priori.beta(est=.35, valor=0.15, prob=0.95, curve=TRUE, tipo="media")


#pri[1]/(pri[1]+pri[2])
#(pri[1]-1)/(pri[1]+pri[2]-2)

#pbeta(0.20, pri[1], pri[2])

#curve(dbeta(x, pri[1], pri[2]), from=0, to=1) 


 
## Posteriori

postBinom <- function(y, size, prioriPars, plot=TRUE, cex.leg=1, ...){
    ## priori deve ser um vetor com os dois parametros da distribuicao Beta
    if(length(prioriPars) != 2) stop("prioriPars deve ter dois elementos")
    postPars <- c(prioriPars[1]+y, prioriPars[2] + size - y)
    mediaPrior <- prioriPars[1]/(prioriPars[1]+prioriPars[2])
    mediaPost = postPars[1]/(postPars[1]+postPars[2])
    modaPrior = ifelse(prioriPars[2] > 1,
    (prioriPars[1]-1)/(prioriPars[1]+prioriPars[2]-2), NA)
    modaPost = ifelse(postPars[2] > 1,
    (postPars[1]-1)/(postPars[1]+postPars[2]-2), NA)
    variaPrior <- prod(prioriPars)/((sum(prioriPars)^2)*(sum(prioriPars)+1))
    variaPost <- prod(postPars)/((sum(postPars)^2)*(sum(postPars)+1))
    pars <- rbind(prioriPars, postPars)
    dimnames(pars) <- list(c("priori","posteriori"), c("alpha","beta"))
    summs <- matrix(c(modaPrior, modaPost, mediaPrior, mediaPost, variaPrior, variaPost),
                    nrow=2, byrow=TRUE)
    dimnames(summs) <- list(c("priori","posteriori"), c("moda","media", "variancia"))
    if(plot){
        par.seq <- seq(0.01,0.99, length=201)
        prior.seq <- dbeta(par.seq, prioriPars[1], prioriPars[2])
        post.seq <- dbeta(par.seq, postPars[1], postPars[2])
        vero.seq   <-  dbeta(par.seq,  y+1,  size-y+1)
        plot(c(0,1), c(0, max(max(prior.seq), max(post.seq), max(vero.seq))), type="n",
                       xlab= expression(theta), ylab= expression(group("[",theta,"]")), ...)
        lines(par.seq, prior.seq, col=2, lty=2, lwd=2)
        lines(par.seq, post.seq, col=4, lwd=2)
        lines(par.seq, vero.seq, col=1, lty=1, lwd=2)
        legend("topright", c("priori","posteriori", "verossimilhança"), lty=c(2,1,1), col=c(2,4,1), cex=cex.leg)
    }
#    return(c(postPars, media=unname(mediaPost), moda=unname(modaPost), EMV=y/size))
    return(list(pars=pars, summary=summs, EMV=y/size))
}

#(post01 <- postBinom(24, 80, pr01))
#(post03 <- postBinom(24, 80, pr03))

#qbeta(c(0.025, 0.975), post01$pars[2,1], post01$pars[2,2])

hpdBeta <- function(pars, prob){
    int.f  <- function(inf, pars, prob){
        p.inf <- pbeta(inf, pars[1], pars[2])
        sup <-  qbeta(prob+p.inf, pars[1], pars[2])
        return(sup - inf) 
    }
    inf.max <- qbeta(1-prob, pars[1], pars[2])
    res.num <- optimize(int.f, c(0,inf.max), pars=pars, prob=prob)
    inf <- res.num$minimum
    sup <-  qbeta(prob+pbeta(inf, pars[1], pars[2]), pars[1], pars[2])
    return(c(inf, sup))
}

#hpdBeta(post01$pars[1,], 0.95)
                        
## rpanel da Beta
