##
## Inferência na distribuição normal
##

## Conjunta:
##f(\mu, \sigma^2|y) = (\sigma^2)^{\frac{n}{2}-1} \exp\left{-\frac{A}{2\sigma^2} \right\}
##    A = SQ + n(\mu - \overline{y})^2
##    SQ = \sum_{i=1}^{n} (y_i - \overline{y})^2
##
## Condicionais
##    [\mu|\sigma^2, y] \sim {\rm N}(\overline{y}, \sigma^2/n)
##    [\sigma^2|\mu, y] \sim {\rm IG}(\frac{n}{2}, \frac{2}{A})
##
## Marginais
##    [\mu|y] \sim {\rm t}_{n-1}(\overline{y}, SQ/n)
##    \frac{\mu - \overline{y}}{\sqrt{sigma^2/n}} \sim {\rm t}_{n-1}
##
##    [\sigma^2|y] \sim {\rm IG}(\frac{n-1}{2}, \frac{2}{SQ})
##    \frac{SQ}{\sigma^2} \sim \chi^2_{n-1}
##

set.seed(20180419)
(y <- rnorm(12, mean=50, sd=8))
dados <- list(n=length(y), m=mean(y), v = var(y), SQ = sum((y-mean(y))^2))
##
## Amostra (exata) da posteriori
##
## para amostrar de pode-se explorar a fatoração:
## [\mu, \sigma^2|y] = [\sigma^2|y] \cdot [\mu|\sigma^2,y] =
## ou, alternativamente
## [\mu, \sigma^2|y] = [\mu|y] \cdot [\sigma^2|\mu,y] =
##
## Vamos adotar aqui a primeira fatoração:
## Obtendo uma amostra
##  (i) Amostrar \sigma^2 de [\sigma^2|y]
(sigma2.sim <- with(dados, 1/rgamma(1, shape=(n-1)/2, scale=2/SQ)))
## (ii) Amostrar \mu de [\mu |\sigma^2,y]
(mu.sim <- with(dados, rnorm(1, mean=m, sd=sqrt(sigma2.sim/n))))
## Obtendo 25.000 amostras
N <- 25000
sigma2.sim <- with(dados, 1/rgamma(N, shape=(n-1)/2, scale=2/SQ))
mu.sim <- with(dados, rnorm(N, mean=m, sd=sqrt(sigma2.sim/n)))

## Gráficos das amostras (correespondem às marginais)
par(mfrow=c(1,2))
t.sim <- with(dados, (mu.sim - m)/sqrt(v/n))
curve(dt(x, df=dados$n-1), from=-4, to=4)
lines(density(t.sim), col=4)
## note a diferença para uma distribuição normal:
curve(dnorm(x), from=-4, to=4, col=2, lty=3, add=TRUE)

chi.sim <- with(dados, SQ/sigma2.sim)
curve(dchisq(x, df=dados$n-1), from=0, to=40)
lines(density(chi.sim), col=4)

##
## Amostra (Gibbs) da posteriori
##
## A estatégia de Gibbs é alternar as simulações entre **as distribuições condicionais**
## o que "parece" errado ,as provouse que a cadeia de valores assim simulados **converge** para a distribuição conjunta
##    [\mu|\sigma^2, y] \sim {\rm N}(\overline{y}, \sigma^2/n)
##    [\sigma^2|\mu, y] \sim {\rm IG}(\frac{n}{2}, \frac{2}{A})
## Obtendo uma amostra
## Como a distribuição de um parâmetro depende da distribuição do outro,
## é necessário fornecer/arbitrar um valor para inicial o algoritmo
mu0 <- 50
##  (i) Amostrar \sigma^2 de [\sigma^2|\mu, y]
A <- with(dados, SQ + n*(mu0 - m)^2)
(sigma2.simG <- with(dados, 1/rgamma(1, shape=n/2, scale=2/A)))
## (ii) Amostrar \mu de [\mu |\sigma^2,y]
(mu.simG <- with(dados, rnorm(1, mean=m, sd=sqrt(sigma2.sim/n))))

## Gerando agora 25.000 amostras
N <- 25000
mu.simG <- sigma2.simG <- numeric(N)
mu.simG[1] <- 30
sigma2.simG[1] <- 100

{for(i in 2:N){
    A <- with(dados, SQ + n*(mu.simG[i-1]-m)^2)
    sigma2.simG[i] <- with(dados, 1/rgamma(1, shape=n/2, scale=2/A))
    mu.simG[i] <- with(dados, rnorm(1, mean=m, sd=sqrt(sigma2.simG[i]/n)))
 }
}

par(mfrow=c(2,1))
plot(mu.simG, type="l")
plot(mu.simG[-(1:1000)], type="l")

plot(sigma2.simG, type="l")
plot(sigma2.simG[-(1:1000)], type="l")

plot(log(sigma2.simG), type="l")
plot(log(sigma2.simG[-(1:1000)]), type="l")

par(mfrow=c(1,2))
t.sim <- with(dados, (mu.sim - m)/sqrt(v/n))
curve(dt(x, df=dados$n-1), from=-4, to=4)
lines(density(t.sim), col=4)
##curve(dnorm(x), from=-4, to=4, col=2, add=TRUE)
t.simG <- with(dados, (mu.simG - m)/sqrt(v/n))
lines(density(t.simG), col=3, lwd=2)

chi.sim <- with(dados, SQ/sigma2.sim)
curve(dchisq(x, df=dados$n-1), from=0, to=40)
lines(density(chi.sim), col=4)
chi.simG <- with(dados, SQ/sigma2.simG)
lines(density(chi.simG), col=3, lwd=2)




## Relação entre Chi-quadrado escalonada inversa e Gamma
a1 <- dados$SQ/rchisq(10000, df=24)
a2 <- 1/rgamma(1000, shape=(dados$n-1)/2, scale=2/dados$SQ)
plot(density(a1))
lines(density(a2), col=2)
