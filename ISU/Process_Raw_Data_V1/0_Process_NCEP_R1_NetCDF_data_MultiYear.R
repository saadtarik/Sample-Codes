# Script to process NCEP NCAR Reanalysis and NCEP/DOE Reanalysis (R1 and R2) data to finer grid
#
# Author: Saad Tarik
# Created: Oct 7, 2017

# Read data
library(raster)
library(ggplot2)
library(ncdf4)
library(rgdal)
library(rgeos)
library(dplyr)

rm(list = ls())

processing_year <- seq(from = 2008, to = 2017, by = 1)

CDD_C_mean_R1 <- stack()

# Start loop through files
for (i in 1:length(processing_year)) {
  fname <- paste0("0_raw/NCEP_NNR/RAW_DATA_Daily/air.2m.gauss.", processing_year[i], ".nc")
  r1temp_stack <- stack(fname) # R1 data
  
  #### Debug Plot
  plot(r1temp_stack[[180]])
  #### Debug Plot 
  
  # Convert to regular grid
  r1temp_stack_reggrid <- rotate(r1temp_stack)
  
  #### Debug Plot
  plot(r1temp_stack_reggrid[[180]])
  #### Debug Plot 
  
  # Crop to CONUS
  CONUS_extent <- extent(-130, -61, 23, 52)
  
  CONUS_r1temp_reggrid <- crop(r1temp_stack_reggrid, CONUS_extent)
  
  #### Debug Plot
  plot(CONUS_r1temp_reggrid[[180]])
  #### Debug Plot 
  
  # Calculate CDDs
  CONUS_r1temp_reggrid_C <- CONUS_r1temp_reggrid - 273.15 # Convert to C
  CONUS_CDD_C_R1 <- CONUS_r1temp_reggrid_C - 18
  CONUS_CDD_C_R1[CONUS_CDD_C_R1 < 0] <- 0
  CONUS_CDD_C_R1_year <- calc(CONUS_CDD_C_R1, fun = sum)
  
  #### Debug Plot
  plot(CONUS_CDD_C_R1_year)
  #### Debug Plot 
  
  # Resample to finer grid
  finer_raster <- raster(nrow = 291, ncol = 691, xmn = -130, xmx = -61, ymn = 23, ymx = 52)
  CONUS_CDD_C_R1_year_finer_res <- raster::resample(CONUS_CDD_C_R1_year, finer_raster, method = "bilinear")
  
  #### Debug Plot
  plot(CONUS_CDD_C_R1_year_finer_res)
  #### Debug Plot
  
  # # CDD stacks for multiple year
  # CDD_C_mean_R1[[i]] <- CONUS_CDD_C_R1_year_finer_res
  # CDD_C_mean_R2[[i]] <- CONUS_CDD_C_R2_year_finer_res
  
  # Write to Raster
  out_fname <- paste0("0_raw/NCEP_NNR/Processed_Data/Long-Term_CDDs/Yearly/CDD_C_R1_2m_", processing_year[i], "_0p1deg.tif")
  writeRaster(CONUS_CDD_C_R1_year_finer_res, filename = out_fname, overwrite = TRUE)
}


r1_CDD_files <- list.files(path = "0_raw/NCEP_NNR/Processed_Data/Long-Term_CDDs/Yearly/", pattern = "*0p1deg.tif$")
r1_CDD_stack <- stack(paste0("0_raw/NCEP_NNR/Processed_Data/Long-Term_CDDs/Yearly/", r1_CDD_files)) # Select the files that you need using indices if you processed more years than you actually want to use

# Calculate mean CDDs for the time period
r1_mean_CDD_C <- calc(r1_CDD_stack, fun = mean)
r1_mean_CDD_C[r1_mean_CDD_C < 0] <- 0
# Visualize
plot(r1_mean_CDD_C)

r1_mean_CDD_F <- r1_mean_CDD_C * 1.8

out_fname <- paste0("1_input/NNR/Mean_CDD_C_", processing_year[1], "_to_", processing_year[length(processing_year)], "_R1.tif")
writeRaster(r1_mean_CDD_C, filename = out_fname, overwrite = TRUE)
