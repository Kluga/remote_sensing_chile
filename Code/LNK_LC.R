#----------------------------Library------------------------------------------#

#install.packages("sf")
#install.packages("spData")
#install.packages("raster")
#install.packages('rgeos')
#install.packages('rgdal')
#install.packages("plyr")
#install.packages("caret")
#install.packages("ranger")
#install.packages("pillar")
#install.packages("ggplot2")
#install.packages("landsat")

library(sf)
library(spData)
library(raster)
library(rgdal)
library(plyr)
library(caret)
library(ranger)
library(pillar)
library(ggplot2)
library(landsat)
library(RStoolbox)

#Workingdirectory

wd <- "C:/Users/Alexander/Desktop/HiWi/Code"
setwd(wd)


#-------------------------------- data ---------------------------------------#


wv2_szene <- brick("../Data/Products/ortho-images/LC_ortho_0.5m.tif") #has to be an ortho-image
DTM <- raster("../Data/Products/dtm/DTM_LC_50cm.tif")                 #has to have the same resolution as the ortho-image

#mask szene and dtm to same extend and shape

DTM[DTM == 300] <- NA 
szene_crop <- crop(wv2_szene, DTM)
DTM_crop <- crop(DTM, szene_crop)
szene_mask <- mask(szene_crop, DTM_crop)
dtm_mask <- mask(DTM_crop, szene_mask[[1]])

writeRaster(szene_mask, filename = "../Data/Workflow/LC/LC_sight.tif", format = "GTiff", overwrite = TRUE)
writeRaster(dtm_mask, filename = "../Data/Workflow/LC/LC_sight_dtm.tif", format = "GTiff", overwrite = TRUE)
#szene_mask <- brick("../Data/Workflow/LC/LC_sight.tif")
#dtm_mask <- raster("../Data/Workflow/LC/LC_sight_dtm.tif")

#------------------ topocorrection to remove shadows -------------------------#

szene_res <- resample(szene_mask, dtm_mask)

crs(dtm_mask) <- CRS('+init=EPSG:32719')      #insert fitting crs
crs(szene_res) <- CRS('+init=EPSG:32719')

dtm_slope <- terrain(dtm_mask, opt="slope", unit="radians")
dtm_aspect <- terrain(dtm_mask, opt="aspect", unit="radians")

szene_topo <- topCor(szene_res, dtm_mask, solarAngles = c(0.9232792, 0.6335545), method = "C")
plotRGB(szene_topo, r=5, g=3, b=2, stretch='lin')

#writeRaster(dtm_slope, "../Data/Indizes/LC/LC_slope.tif", format = "GTiff", overwrite = T)
#writeRaster(dtm_aspect, "../Data/Indizes/LC/LC_aspect.tif", format = "GTiff", overwrite = T)
#writeRaster(szene_topo, "../Data/Workflow/LC/LC_topo.tif", format = "GTiff", overwrite = T)
#szene_topo <- brick('../Data/Workflow/LC/LC_topo.tif')

#-------------------- create random forest model -----------------------------#

ext <- extent(304100, 305200, 6353000, 6354200) #crop to a smaller extent, where the polygons are located

szene_small <- crop(szene_topo, ext)
dtm_small <- crop(dtm_mask, ext)
#dtm_slope <- terrain(dtm_small, opt="slope", unit="radians")
#dtm_aspect <- terrain(dtm_small, opt="aspect", unit="radians")

#writeRaster(szene_small, "../Data/Workflow/LC/LC_small_topo.tif", format = "GTiff", overwrite = T)
#writeRaster(dtm_small, "../Data/Workflow/LC/LC_small_dtm.tif", format = "GTiff", overwrite = T)
#writeRaster(dtm_slope, "../Data/Workflow/LC/LC_small_slope.tif", format = "GTiff", overwrite = T)
#writeRaster(dtm_aspect, "../Data/Workflow/LC/LC_small_aspect.tif", format = "GTiff", overwrite = T)

#Polygone einlesen

poly <- st_read("../Data/Land_cover_classes/LC_Polygons.shp")
poly <- st_transform(poly, crs(szene_small))
names(poly)[names(poly) == "id"] <- "Klasse"
poly$Klasse <- as.character(poly$Klasse)
is.character(poly$Klasse)

plotRGB(szene_small, r=5, g=3, b=2, stretch='lin') #check position of polygons
plot(poly, add = T)

####
# Use "Vegetetaion Indices.R" to create fitting indices for model
####

#load indices

all_files_in_distribution     <- list.files(path = "../Data/Indizes/small/LC/", recursive = T) # List all files
tiff_paths                     <- grep(".tif$", all_files_in_distribution, value=TRUE) # Select tiff-files

number_of_indices <- length(tiff_paths) # All indices

tiff_list <- list()
szene_small_stack <- szene_small
names(szene_small_stack) <- c("coastal_blue", "blue", "green", "yellow", "red", "rededge", "NIR1", "NIR2")

for(i in 1:number_of_indices){ 
  tiff_list[[i]]                 <- raster(paste0("../Data/Indizes/small/LC/", tiff_paths[i]))
  szene_small_stack              <- stack(szene_small_stack, tiff_list[[i]])
}

szene_small_brick <- brick(szene_small_stack)

#writeRaster(szene_small_brick, filename = '../Data/Workflow/LC/LC_small_brick.tif', overwrite = T)
#szene_small_brick <- brick('../Data/Workflow/LC/LC_small_brick.tif')

#extract polygons
ex <- extract(szene_small_brick, poly, df = T)
poly_class <- data.frame(ID = seq(nrow(poly)), Klasse = poly$Klasse)

#create database (link data)

polygons <- join(ex, poly_class, by = 'ID')

#saveRDS(polygons, file = "../Data/Workflow/LC/LC_polygons.rds")
#polygons <- readRDS("../Data/Workflow/LC/LC_polygons.rds")

#create model

set.seed(707)   #(optional) set Seed to make it reproducable (as in Minecraft)

# split Data in 70% trainingdata and 30% Testdata

training <- createDataPartition(y = polygons$Klasse,  #used row
                                times = 1,            #number of splits
                                p = 0.7,              #percentage used for train-data
                                list = FALSE)         #create matrix

df_train <- polygons[training, ]
df_test <- polygons[-training, ]

#saveRDS(df_train, file = "../Data/Workflow/LC/df_train.rds")
#saveRDS(df_test, file = "../Data/Workflow/LC/df_test.rds")

#Cross-Validation
fitControl <- trainControl(
  method = "cv",
  savePredictions = TRUE,
  returnResamp = "final",
  number = 5
  #repeats = 10
  )

tuneControl <- expand.grid(mtry = c(2,3,4,seq(6,length(names(szene_small_brick)),2)))

rf_model <- caret::train (x = df_train[,2:24], 
                          y = df_train$Klasse,
                          method = "rf", metric = "Kappa",
                          importance = TRUE, tuneGrid = tuneControl,
                          trControl = fitControl #to use cross-validation - needs a lot of time
)

rf_model$results
rf_model$finalModel

#saveRDS(rf_model, file = "../Data/Workflow/LC/LC_rf_model.rds")
#rf_model <- readRDS("../Data/Workflow/LC/LC_rf_model.rds")

#feature Selection

FS <- varImp(rf_model)
plot(FS)

# ----------------------- R U N - small --------------------------------------#

names(szene_small_brick)
colnames(df_train)

LNK_small <- raster::predict(object = szene_small_brick, rf_model)

#Normal plotting
plot(LNK_small)

#saveRDS(LNK_small, file = "../Data/Final/Small/LC_small_LNK.rds")
#writeRaster(LNK_small, filename = "../Data/Final/Small/LC_small_LNK.tif", format = "GTiff", overwrite = TRUE)

#LNK_small <- readRDS("../Data/Final/Small/LC_Small_LNK.rds")

LNK_small_ext <- extract(LNK_small, poly, df = T)

LNK_small_test <- LNK_small_ext[-training, ] #seperate predicted polygons

Klasse <- c("11", "12", "21", "22", "31", "32")
layer <- c( 11, 12, 21, 22, 31, 32)

class_layer <- data.frame(Klasse, layer)
LNK_small_Ref <- merge(LNK_small_test, class_layer)

conMat <- confusionMatrix(table(LNK_small_Ref$Klasse, df_test$Klasse))
print(conMat)

#-------------------------- Full R U N ---------------------------------------#

#Use "Vegetetaion Indices.R" to create fitting indices for szene

all_files_in_distribution     <- list.files(path = "../Data/Indizes/LC/", recursive = T) # List all files
tiff_paths                     <- grep(".tif$", all_files_in_distribution, value=TRUE) # Select tiff-files

number_of_indices <- length(tiff_paths) # All indices

tiff_list <- list()
szene_stack <- szene_topo
names(szene_stack) <- c("coastal_blue", "blue", "green", "yellow", "red", "rededge", "NIR1", "NIR2")

for(i in 1:number_of_indices){ 
  tiff_list[[i]]                 <- raster(paste0("../Data/Indizes/LC/", tiff_paths[i]))
  szene_stack                    <- stack(szene_stack, tiff_list[[i]])
}

szene_brick <- brick(szene_stack)

#plot(szene_brick[[1,24]])

#writeRaster(szene_brick, filename = '../Data/Workflow/LC/LC_brick.tif', overwrite = T)

#----------------------------- R U N -----------------------------------------#

names(szene_brick)
colnames(df_train)

LNK <- raster::predict(object = szene_brick, rf_model)

#Normal plotting
plot(LNK)

#saveRDS(LNK, file = "../Data/Final/LC_LNK.rds")
#writeRaster(LNK, filename = "../Data/Final/LC_LNK.tif", format = "GTiff", overwrite = TRUE)

#--------------------------- Validation --------------------------------------#

LNK_extr <- extract(LNK, poly, df = T)

#Testpixel aus den Landnutzungspolygonen auschneiden
LNK_test <- LNK_extr[-training, ]
LNK_Ref <- merge(LNK_test, class_layer)

conMat <- confusionMatrix(table(LNK_Ref$Klasse, df_test$Klasse))
print(conMat)
