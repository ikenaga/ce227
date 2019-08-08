##
## Inferência
## 

##
## Exemplo 1: amostra binomial
##

## Contexto: estimar a proporção $\theta$ de uma população que possui
##           determinado atributo
##
## Emulando "processo" (contexto) no computador
## sob o ponto de vista frequentista

## Geração (artificial) de uma população
th <- 0.17    ## supondo theta conhecido = 0,17

## para testar inicialmente uma pequena população de 100 indivíduos
(POP <- sample(c(0,1), 100, prob=c(1-th, th), rep=TRUE))

## agora sim uma população "grande" com 10 milhões de indivíduos 
set.seed(2019)
POP <- sample(c(0,1), 10000000, prob=c(1-th, th), rep=TRUE)

mean(POP)   ## só para conferir o valor do parâmetro na população...

##
## 1.1 Estimação e distribuição amostral
##

## tirando uma amostra de 80 indivíduos desta população ...
(AM1 <- sample(POP, 80))
## a partir desta amostra podemos estimar o parâmetro
(p1 <- mean(AM1))
## vemos que o valor obtido da amostra é diferente do populacional.

## Além disto ele vai variar de uma amostra para outra como vemos a seguir
## onde tomamos outras amostras de tamanho 80
(AM2 <- sample(POP, 80))
(p2 <- mean(AM2))

(AM3 <- sample(POP, 80))
(p3 <- mean(AM3))

(AM4 <- sample(POP, 80))
(p4 <- mean(AM4))

(AM5 <- sample(POP, 80))
(p5 <- mean(AM4))

(AM6 <- sample(POP, 80))
(p6 <- mean(AM4))

(AM7 <- sample(POP, 80))
(p7 <- mean(AM4))
### ... os seja as estimativas mudam (variam) de uma amostra para outra
### e de forma aleatória. Portanto o estimador (média dos dados o caso = proporção de 1's na amostra)
### As estimativas têm portanto uma distribuição que é chamada de "distribuição amostral" 

### Vamos então obter (nesta simulação computacional) a "distribuição amostral"
### Tomando então "um grande número" (10.000 aqui) de amostras
AMs <- matrix(sample(POP, 80*10000, rep=T), nrow=80)
dim(AMs)
### calculamos a estimativa (a proporção estimada) em cada uma das 10.000 amostras
ps <- apply(AMs, 2, mean)
### e vemos o comportamento das 10000 estimativas:
summary(ps)
hist(ps, prob=TRUE, main="")
lines(density(ps, bw=0.01), lwd=2)

## E podemos obter então várias quantidades de interesse como por exemplo os valores entre os quais estão
## uma determinada fração (95% por exemplo) das proporções obtidas 
abline(v=quantile(ps, prob=c(0.025, 0.975)), lty=2)

## Em certas categorias de problemas existe uma forma analítica (teórica) conhecida (exata ou aproximada)
## No exemplo as proporções seguem uma distribuição aproximadamente normal 
curve(dnorm(x, m=th, sd=sqrt(th*(1-th)/80)), add=T, col=2)

##
## 1.2 Intervalos de confiança
##
## "out of the box":
prop.test(sum(AM1), 80, correct=TRUE)
## mas de onde vem isto? Como interpretar?

## Intervalos de confiança obtidos pela distribuição amostral empírica (obtida por simulação) e teórica:
(ICemp <- quantile(ps, prob=c(0.025, 0.975)))
(ICteo <- qnorm(c(0.025, 0.975), m=th, sd=sqrt(th*(1-th)/80)))

## e os intervalos coincidem!
plot(density(ps, bw=0.01), lwd=2, main="", xlab=expression(hat(theta)), ylab="densidade", ylim=c(0,10))
curve(dnorm(x, m=th, sd=sqrt(th*(1-th)/80)), add=T, col=2)
abline(v=ICemp, lwd=2)
abline(v=ICteo, lty=1, col=2)

##
## OK até aqui,
## mas como saber se não temos as várias amostras da população, mas apenas uma delas?
##
## O problema é que não conhecemos $\theta$, então usando a distribuição "estimada" 
## vamos que ele pode "se desviar" da correta
##
est <- p1
curve(dnorm(x, m=est, sd=sqrt(est*(1-est)/80)), add=T, col=4)
ic <- est + c(-1, 1) * 1.96 * sqrt(est*(1-est)/80)
abline(v=ic, col=4, lty=3)
## ... e agora o intervalo é diferente dos anteriores

## mas como aqui sabemos o valor verdadeiro podemos ver se este valor está ou não no intervalo
curve(dnorm(x, m=est, sd=sqrt(est*(1-est)/80)), col=4, from=0, to=0.5)
abline(v=ic, col=4, lty=3)
arrows(th, 1, th, 0, length=0.15, lwd=2)

## Aqui entra o argumento é "sutil" que leva ao
## considere então onter o intervalo para várias amostras

## para facilitar vamos fazer uma função que obtém uma amostra, calcula o intervalo e adiciona ao gráfico
ic.f <- function(am, plot=FALSE, add=FALSE, ...){
    est <- mean(am)
    ic <- est + c(-1, 1) * 1.96 * sqrt(est*(1-est)/80)
    if(plot){
        curve(dnorm(x, m=est, sd=sqrt(est*(1-est)/80)), add=add, col=1, lty=1, from=0, to=0.5)
        abline(v=ic, ...)
        arrows(th, 1, th, 0, length=0.15, lwd=2)
    }
    return(ic)
}

## e vejamos os intervalos para algumas amostras:
ic.f(AM1, plot=TRUE, col=1, lty=3)
ic.f(AM2, plot=TRUE, col=1, lty=3)
ic.f(AM3, plot=TRUE, col=1, lty=3)
ic.f(AM4, plot=TRUE, col=1, lty=3)
ic.f(AM5, plot=TRUE, col=1, lty=3)
ic.f(AM6, plot=TRUE, col=1, lty=3)
ic.f(AM7, plot=TRUE, col=1, lty=3)

## agora para as milhares de amostras
ICs <- apply(AMs, 2, ic.f)
dim(ICs)
## vamos ver de 20 delas
ICs[, 1:20]
## e note que o intervalo de algumas  delas não contém o valor verdadeiro !!!
## na minha rodasa o 4o, 9o, 13o e 16o
ic.f(ICs[,4], plot=TRUE, col=1, lty=3)


## os intervalos "em geral" contém o valor verdadeiro, mas pode acontecer
## que para algumas amostras não contenha.
## Mais precisamente, a proporção que contém é o "nivel de confiança"
## que foi adotado para construir o intervalo (95%)

## Vamos verificar isto, calculando primeiro se o intervalo de cada amostra
## contém ou não o valor verdadeiro ($\theta=0,20$)

## segue uma função para verificar se o intervalo de cada amostra contém o valor verdadeiro ...
in.f <- function(x, val)
    return(val > x[1] & val < x[2])
## vamos ver para as 10 primeiras amostras
(apply(ICs[,1:20], 2, in.f, val=th))
mean(apply(ICs[,1:20], 2, in.f, val=th))
## e agora calculando a proporção de todas as milhares de amostras que contém o valor verdadeiro
mean(apply(ICs, 2, in.f, val=th))
## deveria ser 95\% ! (diferenças podem ser devido à aproximação e/ou número de simulações)

##
## 1.3 Teste de Hipótese
##

## Vamos considerar agora que há um interesse prático em saber se a proporção está
## acima de um valor de referência 0,20
## Formalmente temos um teste estatístico de hipóteses no qual
##     H_0 : p \leq 0,20 (p0)
##     H_a : p > 0,20  

## usando um procedimento "out of de box" implementado em uma função
prop.test(sum(AM1), 80, p=0.20, alt="greater", correct=FALSE)$conf.int
## mas de onde vem isto? Como interpretar?

p0 <- 0.20
## Inicialmente vamos ver o que ocorre nas milhares de amostras obtidas da população
## das quais extraímos a distribuição amostral empírica
plot(density(ps, bw=0.01), lwd=2, main="", xlab=expression(hat(theta)), ylab="densidade", ylim=c(0, 10))
mean(ps > p0)   ## proporção de amostras com p<0.20 mesmo geradas com theta=0,17 !!
abline(v=p0, lty=3)   

## se conhecessemos $\theta$ a distribuição amostral teórica seria:
curve(dnorm(x, m=th, sd=sqrt(th*(1-th)/80)), add=T, col=2)
1-pnorm(p0, m=th, sd=sqrt(th*(1-th)/80))

## ... mas o problema é que não conhecemos $\theta$,
## e tipícamente, só temos uma amostra...
## No contexto de teste de hipótese usamos a distribuição amostral "sob H_0"
## ou seja para $\theta$ dado pelo valor da hipótese nula
curve(dnorm(x, m=p0, sd=sqrt(p0*(1-p0)/80)), add=T, col=4)

## e vamos onde o valor obtido se posiciona nesta distribuição,
## ou seja, se ele é ou não compatível com a distribuição supondo H_0
(est <- p1) ## aqui p1 = 19/80 = 0.2375
abline(v=est, lty=3, col=4)
## a "compatibilidade" com a distribuição amostral sob $H_0$ é dada pela probabilidade
## (ou proporção) de valores mais extremos que observado,
## e este é o conceito de p-valor!
pnorm(est, m=p0, sd=sqrt(p0*(1-p0)/80), lower=F)

## Note a distinção da distribuição amostral e da distribuição amostal sob H_0 !

## Teste para proporção
sum(AM1)  ## supondo 19
prop.test(sum(AM1), 80, p=0.20, alt="greater", correct=FALSE)
## compare o p-valor da saída e função com a probabilidade calculada anteriormente!

##
## discussão a seguir é opcional
##

## Uma alternativa para teste de hipótese baseado na distribuição amostrar teórica:
## o teste aleatorizado : pode substitui a proporção amostral quando esta
## não é conhecida (ou adequada)
## Não é o caso deste exaemplo mas faremos ainda assim como ilustração

## O primeiro passo é obter a distribuição amostral empírica sob H_0
## simulando da mesmo forma que obtivemos sob o valor verdadeiro
## Simulando a população sob H_0
th0 <- p0             ## supondo $\theta$ definindo em H_0
POP0 <- sample(c(0,1), 10000000, prob=c((1-th0), th0), rep=TRUE)
mean(POP0)   ## só para conferir...

## Obtendo várias amostras (de mesmo tamanho dos dados) desta população 
AM0s <- matrix(sample(POP0, 80*100000, rep=T), nrow=80)
dim(AMs)
p0s <- apply(AM0s, 2, mean)
summary(p0s)
## histograma e densidade empirica das estmativas na amostra
## forneces a distribuição amostral empírica (por simulação)!!!
hist(p0s, prob=TRUE, main="")
lines(density(p0s, bw=0.01), lwd=2)

abline(v=sum(AM1)/80, lty=2, col=2)
(pvalor <- mean(p0s > sum(AM1)/80))  
text(0.25, 8, substitute(pvalor==p, list(p=pvalor)))
## a diferença para o p-valor obtido anteriormente pela distribuição teórica
## pode ser devida a:
## (i) os fato do resultado teórico ser aproximado
## (ii) variação devido à simulação (erro Monte Carlo)

