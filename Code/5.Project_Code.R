################
# Project: Enhancing Island Biogeography: Improving Identification of Potential 
#          Species Pools via Environmental Filtering 
# Authors: Claudia Nuñez-Penichet, Jorge Soberón, Marlon E. Cobos, 
#          Fernando Machado-Stredel, A. Townsend Peterson

# Process: Ecological Niche Modeling with Ellipsoids
################

# Installing and loading packages
#remotes::install_github("marlonecobos/ellipsenm2")
#install.packages("terra")
library(terra)
library(ellipsenm)

# Establishing working directory
#setwd("WORKING DIRECTORY")

# Calling functions (it have to be on the same working directory)
source("Code/Functions_ellip.R")

# Variables for models
vars <- list.files("Data/Rasters", pattern = ".tif", full.names = TRUE)
var_stack <- terra::rast(vars)

# Reading shapefile to crop the variables to area of interest
shp <- vect("Data/Shapefiles/country.shp")

# Cropping the variables
var <- terra::crop(var_stack, shp, mask = T)

# Selecting variables of interest
var_use <- var[[c(3, 4, 1, 2)]]


######
# Reading occurrences
sp <- list.files(path = "Data/Occurrences", pattern = ".csv$", full.names = TRUE)
sp1 <- list.files(path = "Data/Occurrences", pattern = ".csv$", full.names = FALSE)
sp1 <- paste0("Models/", sp1)

errors <- vector()
sensi <- list()

dir.create("Models")

for (i in 1:length(sp1)) {
  sp2 <- read.csv(file = sp[i])
  ext <- extract(var_use, sp2[, c(2, 3)], ID = FALSE)
  
  non_na <- which(is.na(ext), arr.ind = T)[, 1]
  if (length(non_na) > 0) {
    sp2 <- sp2[-non_na, ]
    ext <- ext[-non_na, ]
  }
  
  if (nrow(sp2) > 5) {  # only speices that have more than 5 occurrences
    # initial test
    covm <- cov(ext) # ext is data with which ellipsoids will be created
    
    # testing if positive definite
    test <- is_pos_def(covm)
    
    # testing for non_singularity
    test1 <- non_sing(covm)
    
    # conditional running
    if (test == TRUE & test1 == TRUE) {
      
      sens <- try(ellipsoid_sensitivity(data = ext, variable_columns = 1:4, 
                                        level = 0.95, iterations = 10, 
                                        train_proportion = 0.75))
      
      if (class(sens) != "try-error") {
        sensi[[i]] <- data.frame(Species = sp2[1, 1], sens$summary) 
      } else {
        sensi[[i]] <- NULL
      }
      
      # creating the model with no replicates
      err <- try(ellipsoid_model(data = sp2, species = "Species",
                                 longitude = "Longitude", 
                                 latitude = "Latitude",
                                 raster_layers = var_use, method = "covmat", 
                                 level = 95, replicates = 10, 
                                 percentage = 75,
                                 prediction = "suitability",
                                 return_numeric = FALSE, format = "GTiff",
                                 overwrite = TRUE, 
                                 output_directory = sp1[i]), silent = TRUE)
      if (class(err) == "try-error") {
        errors[i] <- i 
      }
    } else {
      message("\nNon positive definite or singular cov. matrix for species ", 
              as.character(sp2[1, 1]))
    }
  } else {
    message("\nSpecies ", sp1[i], "less than 5 records.............")
  }
}


# Combine sensitivity results excluding those with no results
sensi <- do.call(rbind, sensi[lengths(sensi) > 0])

# write sensitivity results
write.csv(sensi, "sensitivity_results.csv", row.names = FALSE)
