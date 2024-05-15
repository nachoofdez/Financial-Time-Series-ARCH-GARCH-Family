#install.packages("tidyverse") 
library(tidyverse) 
library(lubridate) 
library(car) 
library(urca) 
library(tseries) 
#install.packages("astsa") 
library(astsa) 
library(forecast) 
library(foreign) 
library(timsac) 
#install.packages("vars") 
library(vars) 
library(lmtest) 
#install.packages("mFilter") 
library(mFilter) 
#install.packages("dynlm") 
library(dynlm) 
library(nlme) 
library(quantmod) 
library(xts) 
library(broom) 
#install.packages("rugarch") 
library(rugarch) 
#install.packages("rmgarch") 
library(rmgarch) 
#install.packages("kableExtra")
library(kableExtra) 
library(knitr) 
library(MASS) 
#install.packages("FinTS") 
library(FinTS) 
#install.packages("fGarch") 
library(fGarch) 
library(parallel)
library(readxl)

BBDD_Trabajo_ARCH_GARCH <- read_excel("C:/Users/nacho/OneDrive/Desktop/Master/Series Temporales/Practicas/BBDD Trabajo ARCH.GARCH.xlsx")

Datos_ITALIA <- BBDD_Trabajo_ARCH_GARCH$`31`
Datos_ITALIA<- Datos_ITALIA[2:length(Datos_ITALIA)]

# Convierte los datos a formato numérico
Datos_ITALIA <- as.numeric(Datos_ITALIA)

# Crea una secuencia de fechas desde 1960 a 2022
fechas <- seq(as.Date("1960-01-01"), as.Date("2022-12-31"), by = "year")

# Crea la serie temporal utilizando la biblioteca xts
IPC_Italia <- xts(Datos_ITALIA, order.by = fechas)

# Imprime la serie temporal
plot(IPC_Italia, main = "Serie Temporal de Datos Italia", ylab = "Datos", xlab = "Año")
chartSeries(IPC_Italia)

#Miramos la FAS y la FAP para estimar parametros ARIMA
acf_ipc<-acf(IPC_Italia, lag.max=20)
acp_ipc <- pacf(IPC_Italia,lag.max=20)
#Parece que es un AR(1)

# Revisamos si hay efectos ARCH (cálculo de modelo ARMA/ARIMA, chequeando 
#coeficientes y residuos) 
auto.arima(IPC_Italia, allowdrift=F,trace=T) 
##Falla más que una escopeta de feria!! En realidad tenemos un ARMA (2,0,2) 
fitARIMA <- arima(IPC_Italia, order=c(1,0,0), seasonal = list(order = c(0,0,0), period 
                                                       = 12),method="ML")
coeftest(fitARIMA) 
summary(fitARIMA) 
checkresiduals(fitARIMA)

#Como son significativos, probamos a añadir la media movil
fitARIMA2 <- arima(IPC_Italia, order=c(1,0,1), seasonal = list(order = c(0,0,0), period 
                                                              = 12),method="ML")
coeftest(fitARIMA2) 
summary(fitARIMA2) 
checkresiduals(fitARIMA2)
#Siguen siendo significativos, probamos a meter otra componente 
#autoregresiva

fitARIMA3 <- arima(IPC_Italia, order=c(2,0,1), seasonal = list(order = c(0,0,0), period 
                                                               = 12),method="ML")
coeftest(fitARIMA3) 
summary(fitARIMA3) 
checkresiduals(fitARIMA3)
#Ya no son significativos, asique tenemos un ARMA(1,1)

fitARIMA4 <- arima(IPC_Italia, order=c(1,0,2), seasonal = list(order = c(0,0,0), period 
                                                               = 12),method="ML")
coeftest(fitARIMA4) 
summary(fitARIMA4) 
checkresiduals(fitARIMA4)
#Ya no son significativos, asique tenemos un ARMA(1,1)

fitARIMA5 <- arima(IPC_Italia, order=c(0,0,1), seasonal = list(order = c(0,0,0), period 
                                                               = 12),method="ML")
coeftest(fitARIMA5) 
summary(fitARIMA5) 
checkresiduals(fitARIMA5)


##Calculamos los residuos al cuadrado 
rescuad <- resid(fitARIMA2)^2 
rescuad 
chartSeries(rescuad) 
#Se ve que los residuos al cuadrado son heterocedasticos

#Verificar que tienes algun efecto ARCH. Tiempo vs Res^2 -> Regresión.
fitARIMA.arch <- dynlm(rescuad ~ L(rescuad), data = IPC_Italia) 
summary(fitARIMA.arch)
#p-valor = 0.00041. Tenemos efecto ARCH parece.

## FAS y FAP: si hay barras significativas = no ruido blanco = heterocedasticidad
#Calculamos la FAS 
acfres2=acf(rescuad,lag.max = 200)
#Calculamos la FAP 
pacfres2=pacf(rescuad,lag.max = 200)
#Hay barras significativas, por lo tanto heterocedasticidad
#Sirve para corroborar los ordenes que vienen del ARMA.


## Probamos ARCHTEST
IPCArchTest <- ArchTest(IPC_Italia, lags=1, demean=TRUE)
IPCArchTest

# Tenemos ARCH(2) muy signifativo. También significativo.
IPCArchTest2 <- ArchTest(IPC_Italia, lags=2, demean=TRUE)
IPCArchTest2

# Tenemos ARCH(3) muy signifativo. También significativo.
IPCArchTest3 <- ArchTest(IPC_Italia, lags=3, demean=TRUE)
IPCArchTest3

# Tenemos ARCH(21)  signifativo. También significativo.
IPCArchTest21 <- ArchTest(IPC_Italia, lags=21, demean=TRUE)
IPCArchTest21



##Sabemos que sí hay efectos ARCH, heterocedasticidad 
##¿La varianza depende sólo de los residuos al cuadrado o también de la varianza = GARCH (ARMA)?
#En teoría deberiamos tener un GARCH(1,1), ya que tenemos componente MA:
#Veamos las pruebas:


##Calibramos por defecto GARCH(1,1), equivale a un ARMA (1,1), en RStudio y ver coeficientes
ug_spec=ugarchspec(mean.model = list(armaOrder=c(1,1)))
ug_spec

#Estimamos los coeficientes del modelo (alpha1 = coeficiente residuos^2; beta1 = coeficiente varianza^2)
#Probamos con el modelo por defecto
ugfit=ugarchfit(spec=ug_spec,data=IPC_Italia)
ugfit

## Obtenemos coeficientes del modelo que queremos
ugfit@fit$coef

## Analizamos la varianza 
ug_var=ugfit@fit$var 
ug_var 
## Analizamos los residuos (residuos + línea (varianza)) 
ug_res=(ugfit@fit$residuals)^2 
plot(ug_res,type="l") 
lines(ug_var,col="green")  #Parece un GARCH

# Pronósticos del modelo (ARMA (1,1) + GARCH (1,1)) 
ugfore=ugarchforecast(ugfit,n.ahead=10) 
ugfore 
plot(ugfore)

#Los coeficientes no son significativos, por lo que optamos
# por un modelo ARCH. Tenemos que ver cuál es el orden en el ARCH


#Aunque no tengamos valores significativos en el GARCH.
#Tiene que ser así al tener componente MA.
#El modelo no va a ser todo lo bien que nos gustaría. Si no
#seríamos ricos.



