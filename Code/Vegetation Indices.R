library(raster)

bands <- brick("../Data/Workflow/LC/LC_small_topo.tif")

#plotRGB(bands, 5, 3, 2, stretch="lin")

##getting each band out of the brick for easier calculation

coastal <- bands[[1]] #goastal band of WV 2/3
blue <- bands[[2]]    #glue band of WV 2/3
green <- bands[[3]]   #green band of WV 2/3
yellow <- bands[[4]]  #yellow band of WV 2/3
red <- bands[[5]]     #red band of WV 2/3

rededge <- bands[[6]] #rededge band of WV 2/3
NIR1 <- bands[[7]]    #NIR1 band of WV 2/3
NIR2 <- bands[[8]]    #NIR2 band of WV 2/3


##Calculate the NDVIs of WV (Nouri, Becham et alt.(2014))

NDVI <- (NIR1 - red) / (NIR1 + red)
plot(NDVI)
writeRaster(NDVI, filename = "../Data/Indizes/NDVI", format = "GTiff", overwrite = T)

NDVI2 <- (NIR2 - red) / (NIR2 + red)
plot(NDVI2)
writeRaster(NDVI2, filename = "../Data/Indizes/NDVI2", format = "GTiff", overwrite = T)

NDVI3 <- (NIR2 - yellow) / (NIR2 + yellow)
plot(NDVI3)
writeRaster(NDVI3, filename = "../Data/Indizes/NDVI3", format = "GTiff", overwrite = T)

NDVI5 <- (rededge - red) / (rededge + red)
plot(NDVI5)
writeRaster(NDVI5, filename = "../Data/Indizes/NDVI5", format = "GTiff", overwrite = T)




##Calculate the EVI
##Calculation Parameters for EVI
C1 <- 6     #coefficients of the aerosol resistance term
C2 <- 7.5   #coefficients of the aerosol resistance term
L <- 0.5    #canopy background adjustment
G <- 2.5    #gain factor

EVI <- G * ((NIR2 - red) / (NIR2 + C1 * red + C2 * blue + L))
plot(EVI)
writeRaster(EVI, filename = "../Data/Indizes/EVI", format = "GTiff", overwrite = T)


##SAVI <- Soil Adjusted Vegetation Index (Later, Huete) 
##    (like NDVI, but reducing the influences of soil brightness)

L <- 0.5    #Adjusting the canopy background adjustment for the SAVI

SAVI <- ((1 + L) * (NIR2 - red)) / (NIR2 + red + L)
#plot(SAVI)
writeRaster(SAVI, filename = "../Data/Indizes/SAVI", format = "GTiff", overwrite = T)


##Calculate the NDWI (Method without SWIR)

NDWI <- (green - NIR2) / (green + NIR2)
#plot(NDWI)
writeRaster(NDWI, filename = "../Data/Indizes/NDWI", format = "GTiff", overwrite = T)


##NDBSI <- Normalised Difference Bare Soil Index (Xiaocheng et al.) - blue, coastal

NDBSI <- (blue - coastal) / (blue + coastal)
#plot(NDBSI)
writeRaster(NDBSI, filename = "../Data/Indizes/NDBSI", format = "GTiff", overwrite = T)

##FCI <- Forest and Crop Index - rededge, NIR1

FCI <- (NIR1 - rededge) / (NIR1 + rededge)
#plot(FCI)
writeRaster(FCI, filename = "../Data/Indizes/FCI", format = "GTiff", overwrite = T)

##NDSI <- Normalised Difference Soil Index <- green, yellow

NDSI <- (green - yellow) / (green + yellow)
#plot(NDSI)
writeRaster(NDSI, filename = "../Data/Indizes/NDSI", format = "GTiff", overwrite = T)

##NHFD <- Non-homogenous Feature Difference <- red-edge, coastal

NHFD <- (rededge - coastal) / (rededge + coastal)
#plot(NHFD)
writeRaster(NHFD, filename = "../Data/Indizes/NHFD", format = "GTiff", overwrite = T)

##NDRE <- Normalized Difference Red Edge <- red-edge, NIR (1 / 2)

NDRE <- (NIR1 - rededge) / (NIR1 + rededge)
#plot(NDRE)
writeRaster(NDRE, filename = "../Data/Indizes/NDRE", format = "GTiff", overwrite = T)

NDRE2 <- (NIR2 - rededge) / (NIR2 + rededge)
#plot(NDRE2)
writeRaster(NDRE2, filename = "../Data/Indizes/NDRE2", format = "GTiff", overwrite = T)

##RI <- Ratio NIR/Green Index, good for green accumulations

RI <- NIR1 / green
#plot(RI)
writeRaster(RI, filename = "../Data/Indizes/RI", format = "GTiff", overwrite = T)

##RDVI <- Renormalised Difference Vegetation Index (Roujean, Breon)

RDVI <- (NIR1 - red) / sqrt((NIR1 + red))
plot(RDVI)
writeRaster(RDVI, filename = "../Data/Indizes/RDVI", format = "GTiff", overwrite = T)
