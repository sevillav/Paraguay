#==================================================================
# PROCEDIMIENTO PARA ESTIMAR EL CARBONO ORGANICO EN LOS SUELOS=====
# ============= DE PARAGUAY (KG/M2) Y (TN/HA) ======================
#==================================================================

# MODELO EMPLEADO REGRESION - KRIGING.
# CANTIDAD DE PERFILES DE SUELOS PARA CALIBRACION:2939.
# CANTIDAD DE PERFILES DE SUELOS DEJADOS PARA VALIDACION: 328.

# Establecemos el directorio de trabajo.

setwd("C:/Sevilla/Mario/Paraguay/Corrida")

load("C:/Sevilla/Mario/Paraguay/Corrida/Pry_Anexo_6_Regresion_Kriging.Rdata")

# Verificamos que estamos en el espacio de trabajo requerido:

getwd()

# Cargamos las librerias o paquetes requeridos.

library(raster)
library(car)
library(rgdal)
library(gstat)
library(caret)
library(reshape)
library(sp)
library(lattice)
library(ggplot2)
library(automap)
library(Metrics)

# Cargamos las funciones requeridas.

load("DSM_supportfunctions.RData")

#----- Codigo para seleccionar los perfiles para la validacion---------------------

#Cargamos los datos del splines.

#dat <- read.csv("ParaguayRegMatrix2.csv")

#dat.indices <- sample(1:nrow(dat),size=2939)
#dat.entrenamiento <- dat[dat.indices,]
#dat.test <- dat[-dat.indices,]

#write.csv(dat.entrenamiento, "entrenamiento.csv")
#write.csv(dat.test, "validacion.csv")
#----------------------------------------------------------------------------------

# Cargamos los datos entrenemiento:

dat <- read.csv("entrenamiento.csv")

# Observamos los nombres de los campos o columnas.

summary(dat)
names(dat)

dat$g08esa3a <- NULL
dat$g09esa3a <- NULL
dat$g10esa3a <- NULL
dat$g12esa3a <- NULL
dat$g14esa3a <- NULL
dat$g15esa3a <- NULL
dat$g16esa3a <- NULL
dat$g17esa3a <- NULL
dat$g19esa3a <- NULL
dat$g20esa3a <- NULL
dat$g22esa3a <- NULL
dat$gabhws3a <- NULL
dat$galhws3a <- NULL
dat$ganhws3a <- NULL
dat$gathws3a <- NULL
dat$gclhws3a <- NULL
dat$gcrhws3a <- NULL
dat$ggyhws3a <- NULL
dat$ghshws3a <- NULL
dat$gkshws3a <- NULL
dat$gphhws3a <- NULL
dat$gpthws3a <- NULL
dat$gpzhws3a <- NULL
dat$gumhws3a <- NULL
dat$l01igb3a <- NULL
dat$l15igb3a <- NULL
dat$l16igb3a <- NULL
dat$lmbgsh3a <- NULL
dat$lmtgsh3a <- NULL
dat$cntgad3a <- NULL
dat$g01esa3a <- NULL
summary(dat)

names(dat)

# Realizamos un histogramos de los datos de carbono organico de los perfiles de suelos.

summary(dat$COSkgm2)
hist(dat$COSkgm2, breaks = 150)

# Disenamos un histogramos de los datos de carbono organico de los perfiles de suelos.

hist(dat$COSkgm2, breaks = 150, main = "Histograma de Valores de COS (kg m-2)", ylab = 'Frecuencia', xlab = 'COS (kg m2)',
     sub='Panel (a): Histograma sobre datos de COS de perfiles de suelos.', cex.main=1.2)

# Removemos valores atipicos, segun Bonferroni.

dat <- dat[-c(392, 510, 2531, 1637, 1938, 566, 1880, 413, 2441, 1820, 2395, 711, 1074, 156, 1988, 553, 2360, 1975, 2031, 891,
              1847),]

# Vemos la estructura de los datos.

str(dat)

# Transfomamos a log y realizamos un histogramos de los datos de COS de los perfiles de suelos.

hist(log(dat$COSkgm2), breaks=150)

## Recreamos el objeto con la ubicacion de los puntos

dat_sp <- dat
coordinates(dat_sp) <- ~ LONGITUD + LATITUD

### Analisis de correlacion

names(dat_sp@data)
COR <- cor(as.matrix(dat_sp@data[,6]), as.matrix(dat_sp@data[,-c(1:6)]))
COR
x <- subset(melt(COR), value != 1 | value != NA)
x <- x[with(x, order(-abs(x$value))),]
#as.character(x$X2[1:10])

# Vemos las primeras 10 covariables de mayor correlacion con el COS.

x[1:25,]

idx <- as.character(x$X2[1:25])
idx

# Creamos el archivo de datos para emplear en la regresion lineal multiple.

dat2 <- dat[c('COSkgm2', idx, 'LATITUD', 'LONGITUD')]

# Observamos los nombres de los campos o columnas.

names(dat2)

dat2[dat$COSkgm2 == 0, 1] <- NA

## Modelo de Regresion lineal multiple.

modelo.MLR <- lm(log(COSkgm2) ~ . -LATITUD-LONGITUD, data = dat2)  

# Vemos un resumen de los resultados del modelo de regresion lineal multiple.

summary(modelo.MLR)

# Analisis de varianza.

anova(modelo.MLR)

## Hacemos seleccion de variables por stepwise

modelo.MLR.step <- step(modelo.MLR, direction="both")

summary(modelo.MLR.step)

# Analisis de varianza.

anova(modelo.MLR.step)

# Dividimos el area de graficos en 2 filas y 2 columnas.

par(mfrow=c(2,2)) 

# Hacemos los graficos del modelo de regresion lineal multiple.

plot(modelo.MLR.step)

# Dividimos el area de graficos en 1 filas y 1 columnas.

par(mfrow=c(1,1))

#Falta de multicolinealidad en las variables x: podemos comprobar esto mediante
#el calculo de los Factores de Inflacion de la Varianza (FIVs)

vif(modelo.MLR.step)

#Variables problematicas tienen sqrt(FIV) > 2

sqrt(vif(modelo.MLR.step)) 
 
# Eliminamos del MRL multiple las covariables con Multicolinealidad.

modelo.MLR.step <- update(modelo.MLR.step, . ~ . -px2wcl3a)

# Revisamos de nuevo la multicolinealidad.

sqrt(vif(modelo.MLR.step)) 

#Vamos usar la prueba de Bonferroni para valores atipicos:

outlierTest(modelo.MLR.step)

# Incorporamos las covariables requeridas, segun MRL Multiple.

topo <- stack("PRYtopo.tif")
namesTopo <- readRDS('namesTOPO.rds')
names(topo)
names(topo) <- namesTopo 
names(topo)

# Covariables ambinetales del wordgrid:

cov <- stack('PRY_worldgridsCOVS.tif')
namesCov <- readRDS('worldgridsCOVS_names.rds')
names(cov)
names(cov) <- namesCov
names(cov)

# Apilamos todas las covariales ambientales en un objeto stack (apilado):

COV <- stack(topo, cov) 

# Revisamos los nombres de los campos del objeto COV:

names(COV)

# Seleccionamos solo las primeras 25 covariables de mayor correlacion.

COV <- COV[[idx]]

# Observamos los nombres de los campos o columnas.

names(COV)

# Cambiamos resolucion espacial de las covariables solo para verlos. 
# En la corrida final se debe dejar en la resolucion original de 1 km.

#COV <- aggregate(COV, 10)

# Adecuamos proyeccion cartograficas.
# Project point data 

dat_sp@proj4string <- COV@crs
dat_sp <- spTransform(dat_sp, CRS("+init=epsg:32721"))
COV <- projectRaster(COV, crs = CRS("+init=epsg:32721"), method='ngb')

# Convertimos las covariabes a tabla de datos espaciales.

COV.sp <- as(COV, "SpatialGridDataFrame")

## Eliminamos Datos duplicados.

zerodist(dat_sp)

dat_sp <- dat_sp[dat_sp$COSkgm2 != 0,]

# Ejecutamos estimacion del COS segun ecuacion de RLM Multiple y el kriging de los residuos para la parte
# Continental del Venezuela.

start <- Sys.time()
OCS.krige <- autoKrige(formula = as.formula(modelo.MLR.step$call$formula), 
                       input_data = dat_sp, 
                       new_data = COV.sp,
                       verbose = TRUE,
                       block = c(1000, 1000),
                       model = c("Sph", "Exp"))
print(Sys.time() - start)

# Devolvemos los valores de COS a su condcion original.

RKprediction <- exp(raster(OCS.krige$krige_output[1]))
RKpredsd <- exp(raster(OCS.krige$krige_output[3]))

# Vemos el resumen estadistico de los resultados en kg/m2.

summary(RKprediction)
summary(RKpredsd)

# Si existen valores atipico se pueden eliminar aqui.

#values(RKprediction )[values(RKprediction ) < 0]  <- NA
#values(RKprediction )[values(RKprediction ) > 100]  <- NA
#values(RKpredsd)[values(RKpredsd ) > 10]  <- NA

# Vemos el resumen estadistico de los resultados en kg/m2.

#summary(RKprediction)
#summary(RKpredsd)

# Graficamos los resultados. 

plot(RKprediction)
plot(RKpredsd)

# Reproyectamos la prediccion a geografica WGS84.

RKprediction_geo <- projectRaster(RKprediction, crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), method='ngb')

# Guardamos los resultados en archivos tiff.

writeRaster(RKprediction, filename = "PRY_OCS_RK_kgm2.tif")
#writeRaster(RKprediction, filename = "PRY_OCS_RK_kgm2a.asc")
#writeRaster(RKprediction_geo, filename = "PRY_OCS_RK_kgm2_geo.asc")
writeRaster(RKprediction_geo, filename = "PRY_OCS_RK_kgm2_geot.tif")
writeRaster(RKpredsd, filename = "PRY_OCS_RKpredsd_kgm2.tif")

# Convertimos los resultados de kg/m2 a Tn/ha.
# Importamos el raster resultados

r1 <- raster ('PRY_OCS_RK_kgm2.tif')
r2 <- r1 *10

r3 <- raster ('PRY_OCS_RKpredsd_kgm2.tif')
r4 <- r3 *10

# Graficamos los resultados en Tn/ha.

plot(r2)
plot(r4)

# Resumen del mapa de COS en tn/ha.

summary(r2)
hist(r1, breaks = 100, main  = "Histograma de frecuencia de COS en mapa de RK (kg m-2)", xlab= 'COS (kg m-2)', ylab='Frecuencia', sub='Panel (b): Histograma sobre datos de COS producto de Regresion-Kriging', cex.main=1.2)


hist(r2, breaks = 100, main  = "Histograma de frecuencia de COS en mapa de RK (tn/ha)", xlab= 'COS (tn/ha)', ylab='Frecuencia', sub='Histograma sobre datos de COS producto de Regresion-Kriging')

# Reproyectamos la prediccion a geografica WGS84.

r2_geo <- projectRaster(r2, crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), method='ngb')
r4_geo <- projectRaster(r4, crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), method='ngb')

# Se guarda en formato tif.
writeRaster(r2, 'PRY_Mapa_COS_tnha.tif')
writeRaster(r4, 'PRY_Mapa_COS_Res_tnha.tif')
writeRaster(r4_geo, 'PRY_Mapa_COS_Res_tnha_geo.tif')
writeRaster(r2_geo, 'PRY_Mapa_COS_tnha_geo.tif')

# Estimacion de la incertidumbre segun validacion cruzada.
# Eliminamos datos duplicados.

dat_sp = dat_sp[which(!duplicated(dat_sp@coords)), ]

# Corremos la validacion cruzada.

OCS.krige.cv <- autoKrige.cv(formula = as.formula(modelo.MLR.step$call$formula), 
                             input_data = dat_sp, nfold = 5)

# Vemos un resumen estadistico de la validacion cruzada.

summary(OCS.krige.cv)

#===================================================================================
# ============= SCRIPT PARA VALIDACION EXTERNA =====================================
# ============================== PARAGUAY ==========================================
#===================================================================================

# Para esta validacion se emplearon los 33 puntos dejados fuera de la calibracion
# del modelo de Regresion - Kriging.

# Cargamos los datos de los perfiles de validacion.

datv <- read.csv("validacion1.csv", header = TRUE, sep = ",")

datv <- datv[-c(51),]

# Observamos los nombres de los campos o columnas.

names(datv)

# Vemos un resumen de los datos de COS (Kg/m2) de los perfiles de validacion.

summary(datv$COSkgm2)

## Recreamos el objeto con la ubicacion de los puntos.

coordinates(datv) <- ~ LONGITUD + LATITUD

# Adecuamos proyeccion cartograficas.
# Project point data 

datv@proj4string <- CRS(projargs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
OCSKGM_RK <- raster("PRY_OCS_RK_kgm2_geot.tif")

# Extraemos los datos de COS en kg/m2 de la capa estimada para los puntos de validacion.

datv <- extract(x = OCSKGM_RK, y = datv, sp = TRUE)

# Calculamos la diferencia entre los valores de COS medidos y los COS estims=ados.

datv$PE_RK <- datv$PRY_OCS_RK_kgm2_geot - datv$COSkgm2

# Guardamos los resultados de esta validacion.

#write.csv(datv, "Pry_validacion10.csv", row.names = F)

# Exponemos un resumen de los errores absolutos de prediccion.

summary(res_rk <- abs(datv$PRY_OCS_RK_kgm2_geot - datv$COSkgm2))

# Estimacion de las medidas de calidad del mapa.

# Calculamos el cuartil 75%.

s <- quantile(res_rk,.75, na.rm=TRUE)

# Calculamos e imprimimos el error medio cuadrado entre el valor predicho
# y el valor medido.

a <-(rmse(datv$PRY_OCS_RK_kgm2_geot, datv$COSkgm2))

# Calculamos el R2 entre los valores estimados o predichos y los medidos u observados.

g <- (cor(datv$PRY_OCS_RK_kgm2_geot, datv$COSkgm2)^2)

# Calculamos el Error medio de todos los puntos de validacion.

ME_RK <- mean(datv$PE_RK, na.rm=TRUE)

# Calculamos el error promedio absoluto (MAE).

MAE_RK <- mean(abs(datv$PE_RK), na.rm=TRUE)

# Calculamos el cuadrado del error promedio (MSE).

MSE_RK <- mean(datv$PE_RK^2, na.rm=TRUE)

# Calculamos la raiz cuadrada del error promedio cuadrado (RMSE).

RMSE_RK <- sqrt(sum(datv$PE_RK^2, na.rm=TRUE) / length(datv$PE_RK))

# Estimamos la varianza explicada (Amount of Variance Explained (AVE)).

AVE_RK <- 1 - sum(datv$PE_RK^2, na.rm=TRUE) / 
  sum( (datv$COSkgm2 - mean(datv$COSkgm2, na.rm = TRUE))^2, 
       na.rm = TRUE)

# Impresion de los errores.

metodo <- factor("Regresion-Kriging")
metodo <- data.frame(metodo)

resultados <- cbind(metodo, ME_RK, MAE_RK, MSE_RK, RMSE_RK, AVE_RK, s, g)
etiquetas <- c("Metodos", "ME", "MAE", "MSE", "RMSE", "AVE", "Err Q75", "R2")
names(resultados) <- etiquetas

print(resultados)

# Graficamos las medidas de calidada del mapa.

# Graficamos el Scatter.

par(mfrow=c(1,1))
plot(datv$PRY_OCS_RK_kgm2_geot, datv$COSkgm2, main="Comparacion entre valores COS predichos por Regresion-Kriging y valores reales observados (kg/m2)", xlab=" Valor COS predicho", 
     ylab='Valor COS observado', text(15,0.5, "TEXT"))


# Dibujamos una linea con relacion 1:1 color negro.

abline(0,1, lty=2, col='red')

# Establecemos una linea de regresion entre los valores estimados y los medidos color azul.

abline(lm(datv$COSkgm2 ~ datv$PRY_OCS_RK_kgm2_geot), col = 'blue', lty=2)


legend(x = 11, y = 38,  legend = c("Relacion 1:1", "Regresion"), col = c("Red", "Blue"), 
       title = "Leyenda", lty = 2, cex = 0.75)

text(20,4,"R2 : 0.41", cex = 1.5)

# Graficamos las borbujas de prediccion espacial (Spatial bubbles for prediction errors).

Paraguay <- shapefile("gadm36_PRY_0.shp")

Paraguay@proj4string <- CRS(projargs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

Paraguay.sp <- as(Paraguay, "SpatialPolygonsDataFrame", filled.contour())

bubble(datv[!is.na(datv$PE_RK),], "PE_RK", pch = 21, 
       col=c('red', 'green'), main = "Errores espaciales de prediccion", sp.layout=list("sp.polygons", Paraguay.sp, fill="0") )

# Grabamos el esapcio de trabajo.

save.image("C:/Sevilla/Mario/Paraguay/Corrida/Pry_Anexo_6_Regresion_Kriging.Rdata")
load("C:/Sevilla/Mario/Paraguay/Corrida/Pry_Anexo_6_Regresion_Kriging.Rdata")
