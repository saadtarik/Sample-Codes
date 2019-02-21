# New code for Energyuse estimation from RECS and CBECS data with CDDs
#
# Author: Saad Tarik
# Created: Aug 9, 2018

## V2.0 includes the feedback after thesis defense - uses other covariates slopes
#                                                  - uses UK predictions instead of previously-used RK predicitons

rm(list = ls())
tic <- proc.time()

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Bias
bias <- function(pred,obs){
  n <- length(obs)
  res <- pred - obs
  b <- sum(res, na.rm <- T)/n
  return(b)
}

# RMSE
rmse <- function(pred,obs){
  n <- length(obs)
  res <- pred - obs
  sum_res_sq <- sum(res*res, na.rm = T)
  b <- res/n 
  r <- sqrt(sum_res_sq/n)
  return(r)
}

## Specify user-defined variables ----
font_size <- 18
start_year <- 2008
finish_year <- 2017

library(raster)
library(ggplot2)
library(sensitivity)
library(boot)

## Specify user-defined parameters
domain_of_analysis <- "CensusDiv"
# domain_of_analysis <- "RepDom"

## To shift or not
# shift <- "yes"
shift <- "no"
PerCapita_FlAreaType <- "mean"
# PerCapita_FlAreaType <- "median"
# PerCapita_FlAreaType <- "ratio"
deltax <- 0.1 # Amount of perturbation

NNR_data_source <- "NNR_RAW"
# NNR_data_source <- "NNR_Corr"

## Use lapse rate?
# use_lapse_rate_correction <- "yes"
use_lapse_rate_correction <- "no"

if (use_lapse_rate_correction == "yes") {
  # LapseRate_type <- "Monthly_LapseRate"
  LapseRate_type <- "DryMoist_LapseRate"
  
  NNR_data <- paste0(NNR_data_source, "/", LapseRate_type)
} else {
  NNR_data <- "CDD_RK_10yr_Avg" # NNR_data_source
}

library(xlsx)
if (domain_of_analysis ==  "RepDom") {
  domain_details <- read.csv("1_input/Domain_Details/Results_Summary_V3p1.csv", header = TRUE, sep = ",")
} else if (domain_of_analysis == "CensusDiv") {
  domain_details <- read.csv("1_input/Domain_Details/Results_Summary_V3p2.csv", header = TRUE, sep = ",")
}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## Read impervious surface data
fname <- paste0("1_input/NLCD/CONUS_1km_NLCD.tif")
impsurface <- raster(fname)

impsurface_Full <- impsurface

str(impsurface)

## Read elevation raster data
fname <- paste0("1_input/DEM/SRTM_1km_CONUS_final.tif")
elevation <- raster(fname)

elevation_Full <- elevation

# Create mask
impsurface_temp <- impsurface
impsurface_temp[impsurface_temp < 15] <- NA

library(rgdal)
library(RColorBrewer)
# fname <- domain

## Read outline shapefile and Visualize
library(leaflet)
library(sp)

if (domain_of_analysis == "RepDom") {
  fname <- "1_input/CONUS_outline/CONUS_states_w2009RECSdomainDissolve_prj.shp"
  outline <- shapefile(fname)
  # qtm(outline, fill = "RECS09_Dom")
  nregions <- seq(1, 27, by = 1)
} else if (domain_of_analysis == "CensusDiv") {
  fname <- "1_input/CONUS_outline/CONUS_Subregions_CensusDiv_prj.shp"
  outline <- shapefile(fname)
  # qtm(outline, fill = "SUB_REGION")
  nregions <- seq(1, 9, by = 1)
}

outline_full <- outline

## Read gridded population data 
# fname <- paste0("Z:/Saad/Data_Repository/US_Census/usgrid_data_2000/Pop_CensusDiv/", domain, "_Pop00.tif") # 2000
fname <- paste0("1_input/GRUMP_Population/gpw-v4-population-count_2010_CONUS_Prj.tif") # 2010
popraster <- raster(fname)
popraster_full <- popraster

# popraster <- crop(popraster, impsurface)
# popraster <- mask(popraster, outline)
# popraster[popraster == 0] <- NA
gridded_tot_pop_FULL <- sum(as.vector(popraster_full), na.rm = T)

## Read station data 
library(foreign) # To read .dbf files
fname <- paste0("1_input/GHCN/GHCN_CONUS_LongTerm_MeanCDD_", start_year, "to", finish_year, "w_MeanRatioCDD.csv")
station_data <- read.csv(fname, header = TRUE)

df_all_station_data_orig <- station_data

library(plyr)
library(dplyr)

# fname <- paste0("2_output/RK_Script/",NNR_data,"/CONUS_meanCDD_URB_NTL_", start_year, "_to_", finish_year, "_air2m_F_R1_UHI.tif")
fname <- paste0("2_output/UK_Script/CONUS_meanCDD_", start_year, "_to_", finish_year, "_air2m_F_R1_UHI.tif")
# CONUS_meanCDD_2003_to_2012_air2m_F_R1_UHI
UHI_CDD_Full <- raster(fname)

# fname <-  paste0("2_output/RK_Script/",NNR_data,"/CONUS_meanCDD_URB_NTL_", start_year, "_to_", finish_year, "_air2m_F_R1_non_UHI.tif")
fname <-  paste0("2_output/UK_Script/CONUS_meanCDD_", start_year, "_to_", finish_year, "_air2m_F_R1_non_UHI.tif")
NonUHI_CDD_Full <- raster(fname)

# Nightlight data
fname <- "1_input/NightLight/NTL_Urb_Ext_prj_CONUS.tif"
URB_NTL_Full <- raster(fname)

CDD_Ratio_Full <- UHI_CDD_Full / NonUHI_CDD_Full
CDD_Ratio_Full[CDD_Ratio_Full < 1] <- NA

DiffCDD_Full_percent <- ((UHI_CDD_Full - NonUHI_CDD_Full) / NonUHI_CDD_Full) * 100

DiffCDD_Full <- UHI_CDD_Full - NonUHI_CDD_Full
plot(DiffCDD_Full)

## Read Quantile Regression Stats files
cbecs_qreg_slopes <- read.csv("2_output/EIA_Scripts/CBECS_Cooling_Qreg_Slopes_CenDiv_quantreg.csv")
cbecs_qreg_Covariates_slopes <- read.csv("2_output/EIA_Scripts/CBECS_Cooling_Qreg_Covariates_Slopes_CenDiv_quantreg.csv")
cbecs_qreg_intercepts <- read.csv("2_output/EIA_Scripts/CBECS_Cooling_Qreg_Intercepts_CenDiv_quantreg.csv")

recs_qreg_slopes <- read.csv("2_output/EIA_Scripts/RECS2015_Cooling_Qreg_Slopes_CenDiv_quantreg.csv")
recs_qreg_Covariates_slopes <- read.csv("2_output/EIA_Scripts/RECS_Cooling_Qreg_Covariates_Slopes_CenDiv_quantreg.csv")
recs_qreg_intercepts <- read.csv("2_output/EIA_Scripts/RECS2015_Cooling_Qreg_Intercepts_CenDiv_quantreg.csv")

# Initialize arrays
mean_cbecs_diff_elcl <- rep(NA, length(nregions))
sd_cbecs_diff_elcl <- rep(NA, length(nregions))
mean_recs_diff_elcl <- rep(NA, length(nregions))
sd_recs_diff_elcl <- rep(NA, length(nregions))

mean_cbecs_percent_diff_elcl <- rep(NA, length(nregions))
sd_cbecs_percent_diff_elcl <- rep(NA, length(nregions))
mean_recs_percent_diff_elcl <- rep(NA, length(nregions))
sd_recs_percent_diff_elcl <- rep(NA, length(nregions))

percent_cbecs_diff_elcl <- rep(NA, length(nregions))
percent_recs_diff_elcl <- rep(NA, length(nregions))

total_cbecs_elcl_uhi <- rep(NA, length(nregions))
total_recs_elcl_uhi <- rep(NA, length(nregions))
total_cbecs_elcl_nonuhi <- rep(NA, length(nregions))
total_recs_elcl_nonuhi <- rep(NA, length(nregions))

# Floor area data
floor_area <- read.csv("1_input/Domain_Details/Floor_Area_Cooling_CBECS_RECS_CenDiv.csv")

for (k in 1:length(nregions)) {
  cat("\n Processing region: ", k)
  titletext <- paste0(domain_details$Region_No[k], ": ", domain_details$Region[k])
  
  ## Clip outline and RECS data
  if (domain_of_analysis == "RepDom") {
    outline <- outline_full[which(outline_full@data$RECS09_Dom == k), ]
  } else if (domain_of_analysis == "CensusDiv") {
    outline <- outline_full[which(outline_full@data$CenDiv == k), ]
  }
  
  # Get floor area for the region
  floor_area_region <- floor_area %>% filter(Region_No == k)
  
  # Get data for Census Divisions
  cbecs_qreg_slopes_region <- cbecs_qreg_slopes %>% filter(Region_No == k)
  cbecs_qreg_Covariates_slopes_region <- cbecs_qreg_Covariates_slopes %>% filter(Region_No == k)
  cbecs_qreg_intercepts_region <- cbecs_qreg_intercepts %>% filter(Region_No == k)
  recs_qreg_slopes_region <- recs_qreg_slopes %>% filter(Region_No == k)
  recs_qreg_intercepts_region <- recs_qreg_intercepts %>% filter(Region_No == k)
  recs_qreg_Covariates_slopes_region <- recs_qreg_Covariates_slopes %>% filter(Region_No == k)
  
  cbecs_qreg_Covariates_slopes_region_bldgType <- cbecs_qreg_Covariates_slopes_region %>% filter(Covariate == "as.factor(PBA)02")
  cbecs_qreg_Covariates_slopes_region_people <- cbecs_qreg_Covariates_slopes_region %>% filter(Covariate == "NWKER")
  cbecs_qreg_Covariates_slopes_region_bldgYrCons <- cbecs_qreg_Covariates_slopes_region %>% filter(Covariate == "as.factor(YRCONC)06")
  
  recs_qreg_Covariates_slopes_region_bldgType <- recs_qreg_Covariates_slopes_region %>% filter(Covariate == "as.factor(TYPEHUQ)2")
  recs_qreg_Covariates_slopes_region_people <- recs_qreg_Covariates_slopes_region %>% filter(Covariate == "NHSLDMEM")
  recs_qreg_Covariates_slopes_region_bldgYrCons <- recs_qreg_Covariates_slopes_region %>% filter(Covariate == "as.factor(YEARMADERANGE)4")  
  
  UHI_CDD <- raster::crop(UHI_CDD_Full, outline)
  common_extent <- extent(UHI_CDD)
  # UHI_CDD <- raster::mask(UHI_CDD, outline)
  
  NonUHI_CDD <- raster::crop(NonUHI_CDD_Full, outline)
  # NonUHI_CDD <- raster::mask(NonUHI_CDD, outline)
  # extent(NonUHI_CDD) <- common_extent
  NonUHI_CDD <- resample(NonUHI_CDD, UHI_CDD)
  
  impsurface_domain <- raster::crop(impsurface_Full, outline)
  # impsurface_domain <- raster::mask(impsurface_domain, outline)
  
  elevation_domain <- raster::crop(elevation_Full, outline)
  # elevation_domain <- raster::mask(elevation_domain, outline)
  
  URB_NTL_domain <- raster::crop(URB_NTL_Full, outline)
  # URB_NTL_domain <- raster::mask(URB_NTL_domain, outline)
  # extent(URB_NTL_domain) <- common_extent
  URB_NTL_domain <- resample(URB_NTL_domain, UHI_CDD)
  
  popraster_domain <- raster::crop(popraster_full, outline)
  # popraster_domain <- raster::mask(popraster_domain, outline)
  
  total_pop_domain <- sum(popraster_domain[], na.rm = T)
  per_cap_FloorArea_CBECS_domain <- floor_area_region$CBECS / total_pop_domain
  per_cap_FloorArea_RECS_domain <- floor_area_region$RECS / total_pop_domain
  
  FloorArea_CBECS_ras_domain <- per_cap_FloorArea_CBECS_domain * popraster_domain
  FloorArea_RECS_ras_domain <- per_cap_FloorArea_RECS_domain * popraster_domain
  
  diffCDD <- UHI_CDD - NonUHI_CDD
  diffCDD[diffCDD == 0] <- NA
  
  ratioCDD = UHI_CDD / NonUHI_CDD
  ratioCDD[ratioCDD < 1] = NA
  
  coef_medintensity_cbecs <- c(cbecs_qreg_intercepts_region$Q50, cbecs_qreg_slopes_region$Q50,
                               cbecs_qreg_Covariates_slopes_region_bldgType$Q50, cbecs_qreg_Covariates_slopes_region_people$Q50,
                               cbecs_qreg_Covariates_slopes_region_bldgYrCons$Q50)
  coef_medintensity_recs <- c(recs_qreg_intercepts_region$Q50, recs_qreg_slopes_region$Q50,
                              recs_qreg_Covariates_slopes_region_bldgType$Q50, recs_qreg_Covariates_slopes_region_people$Q50,
                              recs_qreg_Covariates_slopes_region_bldgYrCons$Q50)
  
  median_energyint_cbecs_ras_uhi <- coef_medintensity_cbecs[1] + coef_medintensity_cbecs[2] * UHI_CDD +
    coef_medintensity_cbecs[3] * 1 + coef_medintensity_cbecs[4] * 19 + coef_medintensity_cbecs[5] * 1 # Calculate BTUELCOL intenstity based on CDD for UHI case
  median_energyint_cbecs_ras_nonuhi <- coef_medintensity_cbecs[1] + coef_medintensity_cbecs[2] * NonUHI_CDD +
    coef_medintensity_cbecs[3] * 1 + coef_medintensity_cbecs[4] * 19 + coef_medintensity_cbecs[5] * 1 # Calculate BTUELCOL intenstity based on CDD for Non-UHI case
  
  median_energyint_cbecs_ras_uhi[median_energyint_cbecs_ras_uhi < 0] <- NA
  median_energyint_cbecs_ras_nonuhi[median_energyint_cbecs_ras_nonuhi < 0] <- NA
  
  median_energyint_recs_ras_uhi <- coef_medintensity_recs[1] + coef_medintensity_recs[2] * UHI_CDD +
    coef_medintensity_recs[3] * 1 + coef_medintensity_recs[4] * 2 + coef_medintensity_recs[5] * 1 # Calculate BTUELCOL intenstity based on CDD for UHI case
  median_energyint_recs_ras_nonuhi <- coef_medintensity_recs[1] + coef_medintensity_recs[2] * NonUHI_CDD +
    coef_medintensity_recs[3] * 1 + coef_medintensity_recs[4] * 2 + coef_medintensity_recs[5] * 1 # Calculate BTUELCOL intenstity based on CDD for Non-UHI case
  
  median_energyint_recs_ras_uhi[median_energyint_recs_ras_uhi < 0] <- NA
  median_energyint_recs_ras_nonuhi[median_energyint_recs_ras_nonuhi < 0] <- NA
  
  diff_median_energyuint_cbecs <- median_energyint_cbecs_ras_uhi - median_energyint_cbecs_ras_nonuhi
  diff_median_energyuint_recs <- median_energyint_recs_ras_uhi - median_energyint_recs_ras_nonuhi
  
  # fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_EnergyuseInt/Cooling/CBECS_CenDiv_", k, "_CoolingElecIntDiff_UHI_URB_NTL.tif")
  # writeRaster(diff_median_energyuint_cbecs, filename = fname_ras, overwrite = TRUE)
  # fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_EnergyuseInt/Cooling/RECS_CenDiv_", k, "_CoolingElecIntDiff_UHI_URB_NTL.tif")
  # writeRaster(diff_median_energyuint_recs, filename = fname_ras, overwrite = TRUE)
  
  ## Total electricity from intensity
  median_energy_cbecs_ras_uhi <- median_energyint_cbecs_ras_uhi * FloorArea_CBECS_ras_domain
  median_energy_recs_ras_uhi <- median_energyint_recs_ras_uhi * FloorArea_RECS_ras_domain
  median_energy_cbecs_ras_nonuhi <- median_energyint_cbecs_ras_nonuhi * FloorArea_CBECS_ras_domain
  median_energy_recs_ras_nonuhi <- median_energyint_recs_ras_nonuhi * FloorArea_RECS_ras_domain
  
  ## Mask with URB_NTL data
  median_energy_cbecs_ras_uhi <- raster::mask(median_energy_cbecs_ras_uhi, URB_NTL_domain)
  median_energy_cbecs_ras_nonuhi <- raster::mask(median_energy_cbecs_ras_nonuhi, URB_NTL_domain)
  median_energy_recs_ras_uhi <- raster::mask(median_energy_recs_ras_uhi, URB_NTL_domain)
  median_energy_recs_ras_nonuhi <- raster::mask(median_energy_recs_ras_nonuhi, URB_NTL_domain)
  
  fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/CBECS_CenDiv_", k, "_CoolingElecUse_UHI_Full.tif")
  writeRaster(median_energy_cbecs_ras_uhi, filename = fname_ras, overwrite = TRUE)
  fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/RECS_CenDiv_", k, "_CoolingElecUse_UHI_Full.tif")
  writeRaster(median_energy_recs_ras_uhi, filename = fname_ras, overwrite = TRUE)
  fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/CBECS_CenDiv_", k, "_CoolingElecUse_nonUHI_Full.tif")
  writeRaster(median_energy_cbecs_ras_nonuhi, filename = fname_ras, overwrite = TRUE)
  fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/RECS_CenDiv_", k, "_CoolingElecUse_nonUHI_Full.tif")
  writeRaster(median_energy_recs_ras_nonuhi, filename = fname_ras, overwrite = TRUE)

  total_cbecs_elcl_uhi[k] <- sum(median_energy_cbecs_ras_uhi[], na.rm = T)
  total_recs_elcl_uhi[k] <- sum(median_energy_recs_ras_uhi[], na.rm = T)
  total_cbecs_elcl_nonuhi[k] <- sum(median_energy_cbecs_ras_nonuhi[], na.rm = T)
  total_recs_elcl_nonuhi[k] <- sum(median_energy_recs_ras_nonuhi[], na.rm = T)
  
  diff_median_energy_cbecs <- median_energy_cbecs_ras_uhi - median_energy_cbecs_ras_nonuhi
  diff_median_energy_recs <- median_energy_recs_ras_uhi - median_energy_recs_ras_nonuhi
  
  percent_cbecs_diff_elcl <- (diff_median_energy_cbecs / median_energy_cbecs_ras_uhi) * 100
  percent_recs_diff_elcl <- (diff_median_energy_recs / median_energy_recs_ras_uhi) * 100
  
  mean_cbecs_diff_elcl[k] <- mean(diff_median_energy_cbecs[], na.rm = T)
  mean_recs_diff_elcl[k] <- mean(diff_median_energy_recs[], na.rm = T)
  sd_cbecs_diff_elcl[k] <- sd(diff_median_energy_cbecs[], na.rm = T)
  sd_recs_diff_elcl[k] <- sd(diff_median_energy_recs[], na.rm = T)
  
  mean_cbecs_percent_diff_elcl[k] <- mean(percent_cbecs_diff_elcl[], na.rm = T)
  mean_recs_percent_diff_elcl[k] <- mean(percent_recs_diff_elcl[], na.rm = T)
  sd_cbecs_percent_diff_elcl[k] <- sd(percent_cbecs_diff_elcl[], na.rm = T)
  sd_recs_percent_diff_elcl[k] <- sd(percent_recs_diff_elcl[], na.rm = T)
  
  # fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse_diff_Percent/Cooling/CBECS_CenDiv_", k, "_CoolingEnergy_Percent_Diff_UHI_URB.tif")
  # writeRaster(percent_cbecs_diff_elcl, filename = fname_ras, overwrite = TRUE)
  # fname_ras <- paste0("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse_diff_Percent/Cooling/RECS_CenDiv_", k, "_CoolingEnergy_Percent_Diff_UHI_URB.tif")
  # writeRaster(percent_recs_diff_elcl, filename = fname_ras, overwrite = TRUE)
}

## Mosaic all rasters from the loop
raster_fnames <- list.files("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations", 
                            pattern = "*CoolingElecUse_UHI_Full.tif", 
                            full.names = TRUE, recursive = TRUE)

CBECS_UHI_rasters <- raster_fnames[1:9]
RECS_UHI_rasters <- raster_fnames[10:18]

raster_fnames <- list.files("2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations", 
                            pattern = "*CoolingElecUse_nonUHI_Full.tif", 
                            full.names = TRUE, recursive = TRUE)

CBECS_nonUHI_rasters <- raster_fnames[1:9]
RECS_nonUHI_rasters <- raster_fnames[10:18]

# CBECS UHI
rast.list <- list()
for(i in 1:length(CBECS_UHI_rasters)) { rast.list[i] <- raster(CBECS_UHI_rasters[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
CBECS_Energyuse_UHI_raster <- do.call(mosaic, rast.list)
plot(CBECS_Energyuse_UHI_raster)
fname_ras <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/CBECS_CenDiv_CONUS_CoolingElecUse_UHI_Full.tif"
writeRaster(CBECS_Energyuse_UHI_raster, filename = fname_ras, overwrite = TRUE)

# CBECS non-UHI
rast.list <- list()
for(i in 1:length(CBECS_nonUHI_rasters)) { rast.list[i] <- raster(CBECS_nonUHI_rasters[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
CBECS_Energyuse_nonUHI_raster <- do.call(mosaic, rast.list)
plot(CBECS_Energyuse_nonUHI_raster)
fname_ras <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/CBECS_CenDiv_CONUS_CoolingElecUse_nonUHI_Full.tif"
writeRaster(CBECS_Energyuse_nonUHI_raster, filename = fname_ras, overwrite = TRUE)

# RECS UHI
rast.list <- list()
for(i in 1:length(RECS_UHI_rasters)) { rast.list[i] <- raster(RECS_UHI_rasters[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
RECS_Energyuse_UHI_raster <- do.call(mosaic, rast.list)
plot(RECS_Energyuse_UHI_raster)
fname_ras <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/RECS_CenDiv_CONUS_CoolingElecUse_UHI_Full.tif"
writeRaster(RECS_Energyuse_UHI_raster, filename = fname_ras, overwrite = TRUE)

# RECS non-UHI
rast.list <- list()
for(i in 1:length(RECS_nonUHI_rasters)) { rast.list[i] <- raster(RECS_nonUHI_rasters[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
RECS_Energyuse_nonUHI_raster <- do.call(mosaic, rast.list)
plot(RECS_Energyuse_nonUHI_raster)
fname_ras <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/RECS_CenDiv_CONUS_CoolingElecUse_nonUHI_Full.tif"
writeRaster(RECS_Energyuse_nonUHI_raster, filename = fname_ras, overwrite = TRUE)
