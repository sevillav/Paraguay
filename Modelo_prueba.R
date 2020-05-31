#Mapeo digital de carbono en Paraguay
#Cogido preparado por Mario Guevara
#Victor Sevilla
#Arnulfo Encina

#librerias
library(raster)
library(SuperLearner)
#importar datos
dat <- read.csv('Paraguay_Matriz_Regresion.csv')
#datos de validacion
val <- read.csv('Paraguay_Matriz_Validacion.csv')	
#columnas de interes (predictores)
x <- dat[8:106]
#variable de respuesta
y <- dat$COS_0.30
#datos para validar
VAL <- val[8:106]
#ESTAS VARIABLES SON INCONSISTENTES 
#REVISAR, POR EL MOMENTO LAS ESTOY QUITANDO
VAL$ln1dms3a <- NULL
x$lasmod3a <- NULL
#ensamble de modelos basados en arboles
model <- SuperLearner(		y,
                                x,
                          	SL.library=list("SL.ranger", "SL.extraTrees", 
                                          "SL.xgboost")	)
#prediccion a nuevos datos
VAL$predSL <- predict.SuperLearner(model,  onlySL=TRUE, newdata=VAL)
#validacion de modelo con datos de carbono independientes
VAL$COS_0.30 <- val$COS_0.30

##END ...

cv.model.tune <- CV.SuperLearner(   y,
                                    x,
                                     V=5, SL.library=list("SL.ranger", "SL.extraTrees", 
                                          "SL.xgboost")	)
