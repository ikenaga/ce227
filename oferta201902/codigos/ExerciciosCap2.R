##
## Exercício 2.1
##
## a)
## [Y|\theta] \sim Geo(1-\theta)
## f(y|\theta) = \theta^{y-1} (1-\theta)
## Verossimilhança \propto Be(y, 2)

## Priori Be(p,q)
## f(\theta) = \frac{\theta^{p-1} (1-\theta)^{q-1}}{{\rm Be}(p,q)}

## f(\theta|y) = \theta^{y+p-2} (1-\theta)^{q}
## Posteriori Be(y+p-1, q+1)

## Cenario I:
## p=3, q=5 e y=9
## Priori Be(3,5)
## Verossimilhança \propto Be(9, 2)
## Posteriori Be(11, 6)
curve(dbeta(x, 3, 5), from=0, to=1, col=2, ylim=c(0,4))
curve(dbeta(x, 9, 2), from=0, to=1, add=TRUE)
curve(dbeta(x, 11, 6), from=0, to=1, col=4, add=TRUE)

## Cenario II:
## p=9, q=15 e y=9
## Priori Be(9,15)
## Verossimilhança \propto Be(9, 2)
## Posteriori Be(17, 16)
curve(dbeta(x, 9, 15), from=0, to=1, col=2, ylim=c(0,5))
curve(dbeta(x, 9, 2), from=0, to=1, add=TRUE)
curve(dbeta(x, 17, 16), from=0, to=1, col=4, add=TRUE)

## Vamos ver agora com uma amostra mais "rica"
## Para isto suponha agora uma verrossimilhança Binomal Negativa com k=3
##
## [Y|\theta] \sim Bin(k=3, 1-\theta)
## f(y|\theta) = \binom{y-1}{2} \theta^{y-1} (1-\theta)^3
## Verossimilhança \propto Be(y, 4)

## Priori Be(p,q)
## f(\theta) = \frac{\theta^{p-1} (1-\theta)^{q-1}}{{\rm Be}(p,q)}

## f(\theta|y) = \theta^{y+p-2} (1-\theta)^{q+2}
## Posteriori Be(y+p-1, q+3)

## Cenario III:
## p=3, q=5 e y=27
## Priori Be(3,5)
## Verossimilhança \propto Be(27, 4)
## Posteriori Be(29, 18)
curve(dbeta(x, 3, 5), from=0, to=1, col=2, ylim=c(0,7))
curve(dbeta(x, 27, 4), from=0, to=1, add=TRUE)
curve(dbeta(x, 29, 18), from=0, to=1, col=4, add=TRUE)

## Cenario IV:
## p=9, q=15 e y=9
## Priori Be(9,15)
## Verossimilhança \propto Be(27, 4)
## Posteriori Be(35, 18)
curve(dbeta(x, 9, 15), from=0, to=1, col=2, ylim=c(0,7))
curve(dbeta(x, 27, 4), from=0, to=1, add=TRUE)
curve(dbeta(x, 35, 18), from=0, to=1, col=4, add=TRUE)

##
## Repetindo todos gráficos juntos na mesma escala para comparação
##
par(mfrow=c(2,2), mar=c(3,3,0,0), mgp=c(2,1,0))
curve(dbeta(x, 3, 5), from=0, to=1, col=2, ylim=c(0,7), xlab=expression(theta), ylab="")
curve(dbeta(x, 9, 2), from=0, to=1, add=TRUE)
curve(dbeta(x, 11, 6), from=0, to=1, col=4, add=TRUE)

curve(dbeta(x, 9, 15), from=0, to=1, col=2, ylim=c(0,7), xlab=expression(theta), ylab="")
curve(dbeta(x, 9, 2), from=0, to=1, add=TRUE)
curve(dbeta(x, 17, 16), from=0, to=1, col=4, add=TRUE)

curve(dbeta(x, 3, 5), from=0, to=1, col=2, ylim=c(0,7), xlab=expression(theta), ylab="")
curve(dbeta(x, 27, 4), from=0, to=1, add=TRUE)
curve(dbeta(x, 29, 18), from=0, to=1, col=4, add=TRUE)

curve(dbeta(x, 9, 15), from=0, to=1, col=2, ylim=c(0,7), xlab=expression(theta), ylab="")
curve(dbeta(x, 27, 4), from=0, to=1, add=TRUE)
curve(dbeta(x, 35, 18), from=0, to=1, col=4, add=TRUE)
