#### LIBRARIES ####
library(spatstat); library(GISTools); library(rgdal)
library(fMultivar); library(gstat); library(sp)
library(dplyr); library(tidyr);library(raster)
library(rgeos); library(rgbif); library(viridis)
library(gridExtra); library(rasterVis); library(cleangeo)
library(spdep); library(spgwr); library(tictoc)
library(sf); library(GWmodel); library(ggplot2)


# RData containing:
# Country's polygon (border)
# rabies quarantined ranches per police precincts (rabies)
# data from agricultural census per census tracts (agr.cens)
load("basic_data.RData") 


###### GRID PREPARATION ######
# Function to build the hexagon-cell regular grid
# modified from:
# http://strimas.com/spatial/hexagonal-grids/ Matt Strimas-Mackey, 2016

set.seed(1000)
make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate array of hexagon centers
  g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                offset = c(0.5, 0.5))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
} 


## GRID URUGUAY
hex_grid_500 <- make_grid(border, cell_area = 500000000, clip = TRUE)
hex_grid_500$X <- coordinates(hex_grid_500)[,1]
hex_grid_500$Y <- coordinates(hex_grid_500)[,2]
hex_grid_500$ID <- rownames(hex_grid_500@data)
hex_grid_500$ID.num <- as.numeric(as.character(hex_grid_500$ID))

# Livestock Data Processing
# Already done in the RData
agr.cens$area <- gArea(agr.cens, byid = TRUE)
agr.cens$cattle <- as.integer(as.character(agr.cens$Vacunos)) 
agr.cens$horses <- as.integer(as.character(agr.cens$Equinos)) 
agr.cens$sheep <- as.integer(as.character(agr.cens$Ovinos))
agr.cens$catt.dens <- agr.cens$cattle / (agr.cens$area / 1000000)
agr.cens$hors.dens <- agr.cens$horses / (agr.cens$area / 1000000)
agr.cens$sheep.dens <- agr.cens$sheep / (agr.cens$area / 1000000)
agr.cens$tot.liv.dens <- (agr.cens$cattle + agr.cens$horses + agr.cens$sheep) / (agr.cens$area / 1000000)
agr.cens$catt.biom <- agr.cens$catt.dens * 618
agr.cens$hors.biom <- agr.cens$hors.dens * 403
agr.cens$sheep.biom <- agr.cens$sheep.dens * 39
agr.cens$tot.liv.biom <- agr.cens$catt.biom + agr.cens$hors.biom + agr.cens$sheep.biom

##### RASTERIZATION ##### 
# Create Empty rasters
rst.uy <- raster(ncols = 49, nrows = 54, xmn = rabies@bbox[1,1], 
                 xmx = rabies@bbox[1,2], ymn = rabies@bbox[2,1], 
                 ymx = rabies@bbox[2,2])
rst.uy@crs <- agr.cens@proj4string

# Land cover data
# Land cover layers were intersected with the hexagon/cell grid in ArcGIS 
# to overcome geometric errors
# The Geodatabase from ArcGIS was transformed to an RData file
load("intersect_lccs_hex.RData")
lccs.2000@data <- lccs.2000@data[,-c(1,5:11)]
lccs.2011@data <- lccs.2011@data[,-c(1,13:18)]
lccs.2000$ID.num <- as.numeric(as.character(lccs.2000$ID))
lccs.2011$ID.num <- as.numeric(as.character(lccs.2011$ID))


# Livestock
rst.catt.dens <- rasterize(agr.cens, rst.uy, "catt.dens", background = 0)
rst.hors.dens <- rasterize(agr.cens, rst.uy, "hors.dens", background = 0)
rst.sheep.dens <- rasterize(agr.cens, rst.uy, "sheep.dens", background = 0)
rst.tot.liv.dens <- rasterize(agr.cens, rst.uy, "tot.liv.dens", background = 0)
rst.tot.liv.biom <- rasterize(agr.cens, rst.uy, "tot.liv.biom", background = 0)
# Rabies
rst.rabies.out <- rasterize(rabies, rst.uy, "conteo", background = 0)


#### VECTORIZATION ####
hex_grid_500$liv.biom <- as.vector(raster::extract(rst.tot.liv.biom, 
                                                   hex_grid_500, fun = mean, 
                                                   small = FALSE, na.rm = TRUE, 
                                                   weights = TRUE))
hex_grid_500$liv.dens <- as.vector(raster::extract(rst.tot.liv.dens, 
                                                   hex_grid_500, fun = mean, 
                                                   small = FALSE, na.rm = TRUE, 
                                                   weights = TRUE))
hex_grid_500$cat.dens <- as.vector(raster::extract(rst.catt.dens,
                                                   hex_grid_500, fun = mean, 
                                                   small = FALSE, na.rm = TRUE, 
                                                   weights = TRUE))
hex_grid_500$hor.dens <- as.vector(raster::extract(rst.hors.dens,
                                                   hex_grid_500, fun = mean, 
                                                   small = FALSE, na.rm = TRUE, 
                                                   weights = TRUE))
hex_grid_500$she.dens <- as.vector(raster::extract(rst.sheep.dens, 
                                                   hex_grid_500, fun = mean, 
                                                   small = FALSE, na.rm = TRUE, 
                                                   weights = TRUE))
hex_grid_500$rab.out.dens <- as.vector(raster::extract(rst.rabies.out, 
                                                       hex_grid_500, fun = mean, 
                                                       small = FALSE, na.rm = TRUE, 
                                                       weights = TRUE))
#### LAND COVER ANALYSIS ####
#Livestock suitable areas
lccs.2000$liv.area <- ifelse(lccs.2000$USLB_00 == "Ar" | 
                               lccs.2000$USLB_00 == "ANi" |
                               lccs.2000$USLB_00 == "He" |
                               lccs.2000$USLB_00 == "Pa", 1, 0 )
lccs.2011$liv.area <- ifelse(lccs.2011$USLB_11 == "Ar" | 
                               lccs.2011$USLB_11 == "ANi" |
                               lccs.2011$USLB_11 == "He" |
                               lccs.2011$USLB_11 == "Pa", 1, 0 )
#Commercial afforestation areas
lccs.2000$forest <- ifelse(lccs.2000$USLB_00 == "PF", 1, 0 )
lccs.2011$forest <- ifelse(lccs.2011$USLB_11 == "PF", 1, 0 )

#fragmentation indexes
area.hex <- data.frame("area" = gArea(hex_grid_500, byid = TRUE))
area.hex <- cbind(area.hex, hex_grid_500$ID.num)
names(area.hex)[2] <- "ID.num"

# Livestock areas fragmentation indexes
areas.2000 <- list(NA)
areas.2011 <- list(NA)
for(i in 1:397){
  areas.2000[[i]] <- st_union(st_as_sf(subset(lccs.2000, 
                                              lccs.2000@data$ID.num == i & 
                                                lccs.2000@data$liv.area == 1)))
  areas.2011[[i]] <- st_union(st_as_sf(subset(lccs.2011,
                                              lccs.2011@data$ID.num == i & 
                                                lccs.2011@data$liv.area == 1)))
}

areas.2000.1 <- list(NA)
areas.2011.1 <- list(NA)
for(i in 1:397){
  if (class(areas.2000[[i]])[1] != "sfc_GEOMETRYCOLLECTION") areas.2000.1 <- st_area(st_cast(areas.2000[[i]], "POLYGON")) else areas.2000.1 <- NA
  area.hex$num.patches.2000 [i] <- length(areas.2000.1)
  area.hex$prop.land.2000 [i] <- as.numeric(sum(areas.2000.1))/area.hex$area[i]
  area.hex$med.patch.size.2000 [i] <- as.numeric(mean(areas.2000.1))
  area.hex$eff.mesh.size.2000 [i] <- as.numeric(sum(areas.2000.1^2))/area.hex$area[i]
  if (class(areas.2011[[i]])[1] != "sfc_GEOMETRYCOLLECTION") areas.2011.1 <- st_area(st_cast(areas.2011[[i]], "POLYGON")) else areas.2011.1 <- NA
  area.hex$num.patches.2011 [i] <- length(areas.2011.1)
  area.hex$prop.land.2011 [i] <- as.numeric(sum(areas.2011.1))/area.hex$area[i]
  area.hex$med.patch.size.2011 [i] <- as.numeric(mean(areas.2011.1))
  area.hex$eff.mesh.size.2011 [i] <- as.numeric(sum(areas.2011.1^2))/area.hex$area[i]
}    


#Forest areas indexes
areas.2000f <- list(NA)
areas.2011f <- list(NA)
for(i in 1:397){
  areas.2000f[[i]] <- st_union(st_as_sf(subset(lccs.2000, 
                                               lccs.2000@data$ID.num == i & 
                                                 lccs.2000@data$forest == 1)))
  areas.2011f[[i]] <- st_union(st_as_sf(subset(lccs.2011,
                                               lccs.2011@data$ID.num == i & 
                                                 lccs.2011@data$forest == 1)))
}

areas.2000f.1 <- list(NA)
areas.2011f.1 <- list(NA)
for(i in 1:397){
  if (class(areas.2000f[[i]])[1] != "sfc_GEOMETRYCOLLECTION") areas.2000f.1 <- st_area(st_cast(areas.2000f[[i]], "POLYGON")) else areas.2000f.1 <- NA
  area.hex$num.patches.f.2000 [i] <- length(areas.2000f.1)
  area.hex$prop.land.f.2000 [i] <- as.numeric(sum(areas.2000f.1))/area.hex$area[i]
  if (class(areas.2011f[[i]])[1] != "sfc_GEOMETRYCOLLECTION") areas.2011f.1 <- st_area(st_cast(areas.2011f[[i]], "POLYGON")) else areas.2011f.1 <- NA
  area.hex$num.patches.f.2011 [i] <- length(areas.2011f.1)
  area.hex$prop.land.f.2011 [i] <- as.numeric(sum(areas.2011f.1))/area.hex$area[i]
}    

# Change in fragmentation
hex_grid_500$c.num <- hex_grid_500$num.patches.2011 - hex_grid_500$num.patches.2000
hex_grid_500$c.pla <- hex_grid_500$prop.land.2011 - hex_grid_500$prop.land.2000
hex_grid_500$c.mps <- hex_grid_500$med.patch.size.2011 - hex_grid_500$med.patch.size.2000
hex_grid_500$c.ems <- hex_grid_500$eff.mesh.size.2011 - hex_grid_500$eff.mesh.size.2000

#Change in forest coverage
hex_grid_500$c.fnp <- hex_grid_500$num.patches.f.2011 - hex_grid_500$num.patches.f.2000
hex_grid_500$c.fpl <- hex_grid_500$prop.land.f.2011 - hex_grid_500$prop.land.f.2000


hex_grid_500@data <- merge(hex_grid_500@data, area.hex, by = "ID.num", all.x = TRUE)
hex_grid_500 <- subset(hex_grid_500,
                       is.na(hex_grid_500$rab.out.dens)==FALSE) # Remove one polygon with no data

miss.vals <- unique(unlist (lapply (hex_grid_500@data, function (x) which (is.na (x)))))
hex2 <- hex_grid_500[-miss.vals,] # Remove 7 small border polygons containing missing values

# After data preparation were saved as a RData file
# save.image(file = "Data_for_analysis_final.RData")