############### bibliotecas #################
library(readr)
require(tidyverse)
require(dplyr)
library(car)
library(zoo)
library(forecast)
library(glarma)
library("Hmisc")



#################### Leitura e organizacao dos dados##################


dados <- read_csv("dados_aplicacao (2).csv")
View(dados)

#criando tempo
Tempo <- ts(dados[,1],start=2007,frequency=12) 
n <- length(Tempo)   # Tamanho da serie
print(Tempo)

tempo1 = dados[,1] 

y <- dados$bronq_aguda #variável resposta
x <- dados[,2:9]  #apenas as variáveis atmosféricas
dados['intercept'] <- 1

dados <- dados %>% select(bronq_aguda,intercept,everything())

CO <- dados$CO
PM10 <- dados$PM10
NO <- dados$NO
O3 <- dados$O3
NO2 <- dados$NO2
Temp_min <- dados$Temp_min
RH <- dados$RH
sen12 <- dados$sen12
cos12 <- dados$cos12
sen6 <- dados$sen6
cos6 <- dados$cos6
trend <- dados$trend
intercept <- dados$intercept
##################### modelo mlg ##############################

ajuste_glm <- glm(y ~ NO + O3 + RH + sen12 + cos12 + sen6 + cos6 + trend, family = poisson(link="log"))
summary(ajuste_glm)
acf(ajuste_glm$residuals)
pacf(ajuste_glm$residuals)
plot(ajuste_glm$residuals)

vif(ajuste_glm)



###################### modelo glarma #################
X1 <- data.frame(intercept,NO,RH,sen12,cos12,sen6,cos6,trend) #com NO, RH
X1 <- as.matrix(X1)


ajuste_glarmaX1 <- glarma(y,X1,type="Poi",residuals = "Pearson",phiLags = c(1))
summary(ajuste_glarmaX1)
par(mfrow=c(3,1))
acf(ajuste_glarmaX1$residuals)
pacf(ajuste_glarmaX1$residuals)
plot(ajuste_glarmaX1$residuals)

###################### ajuste previsao ###########
H <- 6


P.ajuste_glarmaX1 <- glarma(Tempo[1:(n-H)], X1[1:(n-H),], phiLags = c(1), thetaLags = NULL, 
                            type = "Poi", residuals = "Pearson")
P.ajuste_glarmaX1
summary(P.ajuste_glarmaX1)
par(mfrow=c(3,1))

acf(P.ajuste_glarmaX1$residuals)
pacf(P.ajuste_glarmaX1$residuals)
plot(P.ajuste_glarmaX1$residuals)




ncov1 <- ncol(X1)




# Previsoes
nP <- length(P.ajuste_glarmaX1$delta) 
nAR <- nP-ncov1	
XP <- matrix(X1[(n-H+1):n,], nrow = H, ncol=ncov1)	
Real <- Tempo[(n-H+1):n]	
rp=residuals(P.ajuste_glarmaX1)

# Usando o modelo com resIduos de Pearson
XL <- X1%*%P.ajuste_glarmaX1$delta[1:ncov1]	

Z0 <- log(P.ajuste_glarmaX1$mu)[n-H]-XL[n-H] 
Z <- NULL
Z[1] <- P.ajuste_glarmaX1$delta[nP]*(Z0+rp[n-H])
for(i in 2:H)
  Z[i] <- P.ajuste_glarmaX1$delta[nP]*Z[i-1]
Prev_MP <- exp(XP%*%P.ajuste_glarmaX1$delta[1:ncov1]+Z)
EQMPp <- sum((Real-Prev_MP)^2)/H
EQMPp


require(forecast)
Ajuste_AR=P.ajuste_glarmaX1$mu
Previsto_AR <- rep(NA,n)

for(i in 1:H)
  Previsto_AR[n-H+i] <- Prev_MP[i]
par(mfrow=c(1,1))
plot(unlist(dados[,1]),type='l',xlab='tempo',ylab='Bronquite')
lines(Ajuste_AR, col='red')
lines(Previsto_AR,type='p', col='red')



############## ajuste_correto################

init = 78
fim = length(y)
Y.prev = NULL
for (s  in init:(fim-1)) {
  y.p = y[1:s]
  x.p = X1[1:s,]
  ajusteglarma = glarma( y.p,x.p, phiLags = c(1),type="Poi",method="FS",residuals = "Pearson")
  x.prev <- matrix(X1[s+1,],nrow = 1)
  Y.prev[s-(init-1)] = forecast(ajusteglarma,1,x.prev)$mu
  
  
  
}


Previsto_AR <- rep(NA,n)


for(i in 1:H)
  Previsto_AR[n-H+i] <- Y.prev[i]


par(mfrow=c(1,1))
plot(unlist(dados[,1]),type='l',xlab='tempo',ylab='Bronquite')
lines(Ajuste_AR, col='red')
lines(Previsto_AR, col='red')

