#### LIBRARIES ####
library(spatstat); library(GISTools); library(rgdal)
library(fMultivar); library(gstat); library(sp)
library(dplyr); library(tidyr);library(raster)
library(rgeos); library(rgbif); library(viridis)
library(gridExtra); library(rasterVis); library(cleangeo)
library(spdep); library(spgwr); library(tictoc)
library(sf); library(GWmodel); library(ggplot2)
library(scatterplot3d); library(PerformanceAnalytics); library(factoextra)


# ******** ANALISIS ******** 
load("Data_for_analysis_final.RData")

# Neighboring matrix
hex.grid.nb.2 <- poly2nb(hex2)
hex.grid.lw.2 <- nb2listw(hex.grid.nb.2)

# Principal Components Analysis for fragmentation information
pc.frag <- princomp(hex2@data[,c(12:19, 26:29)],cor = TRUE) # Excluding forest data
hex2$pc.frag1 <- pc.frag$scores[,1]
hex2$pc.frag2 <- pc.frag$scores[,2]
loadings(pc.frag)

# Alternative PCA using prcomp instead of princomp
pc.frag4 <- prcomp(hex2@data[,c(12:19, 24:27)], center = TRUE, scale. = TRUE)
loadings(pc.frag4)

# Detrending
detrend.1 <- lm(hex2$rab.out.dens ~ hex2$Y * hex2$X)
summary(detrend.1)
hex2$residuals <- residuals(detrend.1)


##### SPATIAL MODELING: SAR #####
#### Spatial Autorregressive Model

# Null Model
mod.sar.null <- spautolm(residuals~1, data=hex2, hex.grid.lw.2, family = "SAR")

### Foward Selection
## STEP 1
pred.var <- c("pc.frag1", "pc.frag2", "liv.biom", "c.fnp", "c.fpl", "temp_jul")
mod = 1
aic.values1 <- NA 
models1 <- list(NA)
best.models <- list(NA)
for(i in 1:length(pred.var)){
  models1[[i]] <- spautolm(as.formula(paste("residuals ~ ", pred.var[i])),
                           data=hex2, hex.grid.lw.2, family = "SAR")
  aic.values1[i] <- AIC(models1[[i]])
}
best.var <- which.min(aic.values1)
best.models[[mod]] <- models1[[best.var]]
prev.var <- pred.var[best.var]

## STEP 2
mod = mod+1
models2 <- list(NA)
aic.values2 <- NA
pred.var <- pred.var[-best.var]
for(i in 1:length(pred.var)){
  models2[[i]] <- spautolm(as.formula(paste0("residuals ~ ", paste(prev.var, pred.var[i], sep= " * "))),
                           data=hex2, hex.grid.lw.2, family = "SAR")
  aic.values2[i] <- AIC(models2[[i]])
}
best.var <- ifelse(min(aic.values2) < AIC(best.models[[1]]) - 2,
                   which.min(aic.values2), NA)
best.models[[mod]] <- models2[[best.var]]
prev.var <- c(prev.var, pred.var[best.var])

## STEP 3
mod = mod+1
models3 <- list(NA)
aic.values3 <- NA
pred.var <- pred.var[-best.var]
for(i in 1:length(pred.var)){
  f1 <- as.formula(paste0("residuals ~ ", 
                          paste(
                            paste(prev.var[1], prev.var[2], pred.var[i], sep = " + "),
                            paste(prev.var[1], prev.var[2],sep = ":"),
                            paste(prev.var[1], pred.var[i],sep = ":"),
                            paste(prev.var[2], pred.var[i],sep = ":"), 
                            sep = " + ")))
  models3[[i]] <- spautolm(f1, data=hex2, hex.grid.lw.2, family = "SAR")
  aic.values3[i] <- AIC(models3[[i]])
}
best.var <- ifelse(min(aic.values3) < AIC(best.models[[2]]) - 2,
                   which.min(aic.values3), NA)
best.models[[mod]] <- models3[[best.var]]
prev.var <- c(prev.var, pred.var[best.var])

#### AIC and AIC Weights *including null models*
aic1 <- c(14.44, 8.22, 0.00)
# weights
exp(-aic1/2) / sum(exp(-aic1/2))
rm(aic1)

## Residual's Autocorrelation -Best Model-
geary.mc(best.models[[2]]$fit$residuals, hex.grid.lw.2, 1000, 
         alternative = "greater") # The residuals show no significant spatial autocorrelation

## Incorporate the data in the main object
hex2$pred.SAR <- best.models[[2]]$fit$fitted.values
hex2$res.SAR <- best.models[[2]]$fit$residuals
hex2$pred.LM.SAR <- hex2$pred.SAR + hex2$fitted_linear

##### SPATIAL MODELING: GWR #####
##### Geographically Weighted Regression

mod.gwr.null <- gwr.basic(residuals~1, hex2, bw = 6, adaptive = TRUE, kernel = "boxcar")

## STEP 1
pred.var.gwr <- c("pc.frag1", "pc.frag2", "liv.biom", "c.fnp", "c.fpl", "temp_jul")
mod.gwr1 <- list(NA)
aic.gwr1 <- NA
best.gwr <- list(NA)
for(i in 1:length(pred.var.gwr)){
  mod.gwr1[[i]] <- gwr.basic(as.formula(paste("residuals ~ ", pred.var.gwr[i])), 
                             hex2, bw = 6, adaptive = TRUE, kernel = "boxcar")
  aic.gwr1[i] <- AIC(mod.gwr1[[i]]$lm)
}
best.var.gwr <- which.min(aic.gwr1)
best.gwr[[1]] <- mod.gwr1[[best.var.gwr]]
prev.var.gwr <- pred.var.gwr[best.var.gwr]

## STEP 2
mod.gwr2 <- list(NA)
aic.gwr2 <- NA
pred.var.gwr <- pred.var.gwr[-best.var.gwr]
for(i in 1:length(pred.var.gwr)){
  mod.gwr2[[i]] <- gwr.basic(as.formula(paste0("residuals ~ ", paste(prev.var.gwr, pred.var.gwr[i], sep= " + "))), 
                             hex2, bw = 6, adaptive = TRUE, kernel = "boxcar")
  aic.gwr2[i] <- AIC(mod.gwr2[[i]]$lm)
}
best.var.gwr <- ifelse(min(aic.gwr2) < AIC(best.gwr[[1]]$lm) - 2,
                       which.min(aic.gwr2), NA)
best.gwr[[2]] <- mod.gwr2[[best.var.gwr]]
prev.var.gwr <- paste(prev.var.gwr, pred.var.gwr[best.var.gwr], sep= " + ")

##STEP 3
mod.gwr3 <- list(NA)
aic.gwr3 <- NA
pred.var.gwr <- pred.var.gwr[-best.var.gwr]
for(i in 1:length(pred.var.gwr)){
  mod.gwr3[[i]] <- gwr.basic(as.formula(paste0("residuals ~ ", paste(prev.var.gwr, pred.var.gwr[i], sep= " + "))), 
                             hex2, bw = 6, adaptive = TRUE, kernel = "boxcar")
  aic.gwr3[i] <- AIC(mod.gwr3[[i]]$lm)
}
best.var.gwr <- ifelse(min(aic.gwr3) < AIC(best.gwr[[2]]$lm) - 2,
                       which.min(aic.gwr3), NA)
best.gwr[[3]] <- mod.gwr3[[best.var.gwr]]
prev.var.gwr <- paste(prev.var.gwr, pred.var.gwr[best.var.gwr], sep= " + ")

## STEP 4
mod.gwr4 <- list(NA)
aic.gwr4 <- NA
pred.var.gwr <- pred.var.gwr[-best.var.gwr]
for(i in 1:length(pred.var.gwr)){
  mod.gwr4[[i]] <- gwr.basic(as.formula(paste0("residuals ~ ", paste(prev.var.gwr, pred.var.gwr[i], sep= " + "))), 
                             hex2, bw = 6, adaptive = TRUE, kernel = "boxcar")
  aic.gwr4[i] <- AIC(mod.gwr4[[i]]$lm)
}
best.var.gwr <- ifelse(min(aic.gwr4) < AIC(best.gwr[[3]]$lm) - 2,
                       which.min(aic.gwr4), NA)
best.gwr[[4]] <- mod.gwr4[[best.var.gwr]]
prev.var.gwr <- paste(prev.var.gwr, pred.var.gwr[best.var.gwr], sep= " + ")

## STEP 5
mod.gwr5 <- list(NA)
aic.gwr5 <- NA
pred.var.gwr <- pred.var.gwr[-best.var.gwr]
for(i in 1:length(pred.var.gwr)){
  mod.gwr5[[i]] <- gwr.basic(as.formula(paste0("residuals ~ ", paste(prev.var.gwr, pred.var.gwr[i], sep= " + "))), 
                             hex2, bw = 6, adaptive = TRUE, kernel = "boxcar")
  aic.gwr5[i] <- AIC(mod.gwr5[[i]]$lm)
}
best.var.gwr <- ifelse(min(aic.gwr5) < AIC(best.gwr[[4]]$lm) - 2,
                       which.min(aic.gwr4), NA)
best.gwr[[5]] <- mod.gwr5[[best.var.gwr]]
prev.var.gwr <- paste(prev.var.gwr, pred.var.gwr[best.var.gwr], sep= " + ")

#### AIC and AIC Weights *including null models*
aic1 <- c(35.41, 15.09, 7.14, 2.75, 0.00)
# weights
exp(-aic1/2) / sum(exp(-aic1/2))
rm(aic1)

aic1 <- c(162.46, 282.37, 255.74, 228.05, 0.00)
# weights
exp(-aic1/2) / sum(exp(-aic1/2))
rm(aic1)

## Residual's Autocorrelation -Best Model-
geary.mc(best.gwr[[4]]$SDF$residual, hex.grid.lw.2, 1000, 
         alternative = "greater") # The residuals show no significant spatial autocorrelation

## Variable coefficients' Autocorrelation -Best Model-
geary.mc(best.gwr[[4]]$SDF$pc.frag1, hex.grid.lw.2, 1000, 
         alternative = "greater")
geary.mc(best.gwr[[4]]$SDF$c.fnp, hex.grid.lw.2, 1000, 
         alternative = "greater")
geary.mc(best.gwr[[4]]$SDF$c.fpl, hex.grid.lw.2, 1000, 
         alternative = "greater")
geary.mc(best.gwr[[4]]$SDF$temp_jul, hex.grid.lw.2, 1000, 
         alternative = "greater")

## Incorporate the data in the main object
hex2$pred.GWR <- best.gwr[[4]]$SDF$yhat
hex2$res.GWR <- best.gwr[[4]]$SDF$residual
hex2$pred.LM.GWR <- hex2$pred.GWR + hex2$fitted_linear

# save.image("data_models_SAR_GWR_result.RData")

