#### LIBRARIES ####
library(spatstat); library(GISTools); library(rgdal)
library(fMultivar); library(gstat); library(sp)
library(dplyr); library(tidyr);library(raster)
library(rgeos); library(rgbif); library(viridis)
library(gridExtra); library(rasterVis); library(cleangeo)
library(spdep); library(spgwr); library(tictoc)
library(sf); library(GWmodel); library(ggplot2)
library(scatterplot3d); library(PerformanceAnalytics); library(factoextra)


# ****** FIGURES ******
load("data_models_SAR_GWR_result.RData")

hex2.sf <- st_as_sf(hex2) # transformation into a SF object
names(hex2.sf)

#### FIG 1 #### 
# Map of the country and distribution of outbreaks
# created in QGIS from the spatial raw data


#### FIG 2 ####
# predictor variables used in the spatial modeling 
# panels compiled into one figure in Inkscape
plot(hex2.sf[31], main = NA) # pc.frag1: current fragmentation 
plot(hex2.sf[32], main = NA) # pc.frag2: recent change in fragmentation 
plot(hex2.sf[28], main = NA) # c.fnp: recent change in number of forestry patches  
plot(hex2.sf[29], main = NA) # c.fpl: recent change in proportion of landscape occupied by forestry 
plot(hex2.sf[5], main = NA) # liv.biom: livestock biomass 
plot(hex2.sf[30], main = NA) # temp_jul: average minimum temperature in the coldest month 


#### FIG 3 ####
# Observed density of outbreaks and predictions from linear and spatial models
# panels edited in Inkscape
plot(hex2.sf[c(10,34,37,40)], breaks = seq(-0.5, 6.5, 0.5))


#### FIG S1 ####
# Time series of rabies outbreaks
load("outbreaks_dates.RData")
hist(outbreaks$date, "weeks", freq = TRUE, main=NA, xlab = NA)


#### FIG S2 ####
# Fragmentation variables correlation plot
chart.Correlation(hex2@data[,c(12:19, 24:27)], histogram=TRUE, pch=19)


#### FIG S3 ####
# PCA biplot
plot(pc.frag$scores[,1], pc.frag$scores[,2], 
     type = "p", pch = 20, xlab= "pc.frag1 (46.74%)", 
     ylab = "pc.frag2 (20.94%)", 
     xlim = c(-10,10), ylim = c(-10,10), cex = 0.6)
arrows(0,0,pc.frag$loadings[,1]*15, 
       pc.frag$loadings[,2]*15, length = 0.05, 
       col = "blue")
text(pc.frag$loadings[,1]*15, 
     pc.frag$loadings[,2]*15,
     dimnames(pc.frag$loadings)[[1]], 
     cex=0.8, col = "blue")


#### FIG S4 ####
# comparison between both PCA methods
# panels edited in Inkscape
# summary(pc.frag4) # variances: 46.7%; 20.9%
fviz_pca_biplot(pc.frag4, geom.ind = "point", title = "A: prcomp", 
                xlab = "pc.frag1 (46.7%)", ylab = "pc.frag2 (20.9%)")
# summary(pc.frag) # variances: 46.7%; 20.9%
fviz_pca_biplot(pc.frag, geom.ind = "point", title = "B: princomp", 
                xlab = "pc.frag1 (46.7%)", ylab = "pc.frag2 (20.9%)")


#### FIG S5 ####
# First order effect (linear trend)
plot(hex2.sf[34], main = NA)


#### FIG S6 ####
# Moran Plot for the residuals of the first order model (linear trend)
# modified from spdep::moran.plot
plot(hex2$residuals, lag.listw(hex.grid.lw.2, hex2$residuals), 
     xlab = 'Detrended rabies density', 
     ylab = 'Spatially lagged density', pch = 20)
abline(lm(lag.listw(hex.grid.lw.2, hex2$residuals) ~ hex2$residuals))
abline(h=0, v=0, lty = 3)


#### FIG S7 ####
# Geary's C Correlogram for the residuals of the first order model (linear trend)
plot(sp.correlogram(hex.grid.nb.2, hex2$residuals, order = 10, method = "C", randomisation = TRUE), main = NA)


#### FIG S8 ####
# residuals from the SAR
plot(hex2.sf[36], main = NA)


#### FIG S9 ####
# distribution of fitted coefficients for the best ranking GWR model
# panels compiled into one figure in Inkscape
gwr.4sf <- st_as_sf(best.gwr[[4]]$SDF) 
names(gwr.4sf)
plot(gwr.4sf[2], main = NA) # pc.frag1: current fragmentation
plot(gwr.4sf[3], main = NA) # c.fnp: recent change in number of forestry patches  
plot(gwr.4sf[4], main = NA) # c.fpl: recent change in proportion of landscape occupied by forestry
plot(gwr.4sf[5], main = NA) # temp_jul: average minimum temperature in the coldest month 


#### FIG S10 ####
# scatterplots for fitted coefficients for the best ranking GWR model
# panels compiled into one figure in Inkscape
par(mfrow=c(2,2))
plot(hex2.sf$pc.frag1, gwr.4sf$pc.frag1, xlab=NA, ylab="fitted coefficients", main=NA, pch=20, cex=0.5)
plot(hex2.sf$c.fnp, gwr.4sf$c.fnp, xlab=NA, ylab=NA, main=NA, pch=20, cex=0.5)
plot(hex2.sf$c.fpl, gwr.4sf$c.fpl, xlab="values", ylab="fitted coefficients", main=NA, pch=20, cex=0.5)
plot(hex2.sf$temp_jul, gwr.4sf$temp_jul, xlab="values", ylab=NA, main=NA, pch=20, cex=0.5)
par(mfrow=c(1,1))


#### FIG S11 ####
# residuals from the GWR
plot(hex2.sf[39], main = NA)

#### FIG S12 ####
# Performance of spatial models: predicted values from SAR and GWR vs. observed density of rabies outbreaks
par(mfrow=c(1,2))
plot(hex2$pred.LM.SAR, hex2$rab.out.dens, xlab = "Predicted values (LM+SAR)", ylab = "Rabies outbreaks density", pch = 20,
     main="A", xlim = c(0,6))
abline(a=-0.03881, b= 1.16305)
plot(hex2$pred.LM.GWR, hex2$rab.out.dens, xlab = "Predicted values (LM+GWR)", ylab = NA, pch = 20, 
     main="B", xlim = c(0,6))
abline(a=-0.00495, b= 1.00319)
par(mfrow=c(1,1))
