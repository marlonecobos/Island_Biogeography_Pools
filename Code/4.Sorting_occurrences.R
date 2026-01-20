################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Eliminating the occurrences from the Caribbean Islands
################

setwd("WORKING DIRECTORY")

# loading packages
#install.packages("terra")
library(terra)

papi <- list.files(path = "Thinned/Papilionidae",
                   pattern = ".csv$", full.names = TRUE, recursive = T)
papi_names<- list.files(path = "Thinned/Papilionidae",
                        pattern = ".csv$", full.names = F, recursive = T)
papi_names <- paste0("Final/", papi_names)

carib <- vect("Data/Shapefiles/country.shp")

for (i in 1:length(papi_names)){
  occ <- read.csv(papi[i])
  occ1 <- occ[, c("acceptedScientificName", "decimalLongitude", "decimalLatitude")]
  
  occ1 <- vect(occ1, geom = c("decimalLongitude", "decimalLatitude"), 
               crs = crs(carib), keepgeom = TRUE)
  
  occ1 <- occ1[carib, ]
  
  occ1 <- as.data.frame(occ1)
  
  if (nrow(occ1) > 0){
    write.csv(occ1, papi_names[i], row.names = F)
  }
}

###-------------------Pieridae-------------------
pier <- list.files(path = "Thinned/Pieridae",
                   pattern = ".csv$", full.names = TRUE, recursive = T)
pier_names<- list.files(path = "Thinned/Pieridae",
                        pattern = ".csv$", full.names = F, recursive = T)
pier_names <- paste0("Final/", pier_names)

for (i in 1:length(pier_names)){
  occ <- read.csv(pier[i])
  occ1 <- occ[, c("scientificName", "decimalLongitude", "decimalLatitude")]
  
  occ1 <- vect(occ1, geom = c("decimalLongitude", "decimalLatitude"), 
               crs = crs(carib), keepgeom = TRUE)
  
  occ1 <- occ1[carib, ]
  
  occ1 <- as.data.frame(occ1)
  
  if (nrow(occ1) > 0){
    write.csv(occ1, pier_names[i], row.names = F)
  }
}
