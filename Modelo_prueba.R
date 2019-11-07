#Codigo base para ensamble de modelos usando Super Learner

###10cm
library(raster)
library(SuperLearner)
DAT <- readRDS('DAT_2001_2005_SOC.rds')
preds <- readRDS('prediction_domain_2001_2005.rds')
e <- read.csv('prediction_domain_example.csv')

PREDICTED <- stack()
#i = 2001
for (i in 2001:2005){
DAT1 <- DAT[DAT$year==i,]
y <- DAT1$SOC
x <- DAT1[1:41]
model <- SuperLearner(
				y,
                                x,
                          	SL.library=list("SL.ranger", "SL.extraTrees", 
                                          "SL.xgboost"))
predsA <- preds[preds$year==i,]
pred<- predict.SuperLearner(model,  onlySL=TRUE, newdata=predsA[names(x)])
predAll <- data.frame(e , pred[[1]])
coordinates(predAll) <- ~ x+y
gridded(predAll) <- TRUE
predAll <- raster(predAll)
names(predAll) <- paste0('predicted', i)
PREDICTED <- stack(PREDICTED, predAll)
print(i)
}

