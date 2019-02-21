## Script for UK in R and energyuse comparison
# 
# Author: Saad Tarik
# Created: Nov 9, 2018

## V1.0
# Unviersal Kriging of CDDs based on RK script V15.0 (RK_and_Energyuse_Script_with_CDD_V15.R)

# setwd("Z:/Saad/Data_Repository/UHI_EnergyUse") # Change Z: to the drivel letter assigned in your computer

rm(list = ls())
tic = proc.time()

## Specify user-defined variables
start_year = 2008
finish_year = 2017

## Statistics (User-defined functions)

# Bias
bias = function(pred,obs){
  n = length(obs)
  res = pred - obs
  b = sum(res, na.rm = T)/n
  return(b)
}

# RMSE
rmse = function(pred,obs){
  n = length(obs)
  res = pred - obs
  sum_res_sq = sum(res*res, na.rm = T)
  b = res/n 
  r = sqrt(sum_res_sq/n)
  return(r)
}

rep_domain = seq(1, 27, by = 1) # represents the reportable domain from the RECS 2009

## Reanalysis Data
reanalysis = 'R1' ### DO NOT CHANGE TO 'R2'- KEEP 'R1'
# reanalysis = 'R2'

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## Read .tif files (Raster files - impervious surface and elevation)
library(raster)
# Impervious surface
fname = paste0("1_input/NLCD/CONUS_1km_NLCD.tif")
impsurface = raster(fname)
str(impsurface)
impsurface_mat = as.matrix(impsurface) # Convert to matrix
impsurface_df = as.data.frame(impsurface) # Convert to data frame
ind_impsurface = which(!is.na(impsurface_mat))
sp.studygrid = as(impsurface, "SpatialGridDataFrame") # Convert to spatial grid
# sp.studygrid2 = as(impsurface, "SpatialPixelDataFrame")
class(sp.studygrid)

# Elevation
# NOTE: The path after Data_Repository should be unchanged. Change the path before Data_Repository with your drive letter.
fname = paste0("1_input/DEM/SRTM_1km_CONUS_final.tif")
elevation = raster(fname)
elevation = crop(elevation, impsurface)


## Nightlight Urban data
fname = paste0("1_input/NightLight/NTL_Urb_Ext_prj_CONUS.tif")
urb_extent_NTL = raster(fname)
urb_extent_NTL = crop(urb_extent_NTL, impsurface)

## Read outline shapefile
library(rgdal)
library(sp)
fname = "1_input/CONUS_outline/CONUS_outline_Dissolve_prj.shp"
outline = shapefile(fname)
#### Debug ####
# out2 = readOGR("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Stations","QC_GHCN_2009_Stations_CONUS_prj_ImpElevNCEP_State_Join_v4")
# plot(out2)
#### Debug ####

outline_full = outline

## CDD from NCEP
# NOTE: The path after Data_Repository should be unchanged. Change the path before Data_Repository with your drive letter.
# fname = "C:/Users/starik/Documents/Data_Repository/NCEP_NNR/air2m/CDD_2009_0p1deg_mask_prj.tif"
if (reanalysis == 'R1') {
  fname = paste0("1_input/NNR/Mean_CDD_C_", start_year,"_to_", finish_year, "_R1.tif")
  cdd_full = raster(fname)
} else if (reanalysis == 'R2') {
  fname = "1_input/NCEP-DOE_R2/Mean_CDD_C_2000_to_2009_R2_prj_mask.tif" # R2 air temp
  cdd_full = raster(fname)
}
cdd_full_prj = raster::projectRaster(cdd_full, 
                                     crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
cdd_full_prj_resample = raster::resample(cdd_full_prj, impsurface, method = "bilinear")
## Crop CDD from NCEP to study area
cdd_cropped = crop(cdd_full_prj_resample, extent(impsurface)) # cdd_full = cdd_nonuhi
cdd_nonuhi = cdd_cropped
cdd_nonuhi = mask(cdd_nonuhi, outline)

library(dplyr)
library(ggmap)
library(ggplot2)
library(maptools)
library(RColorBrewer)
library(gstat)
library(automap)

## Read shapefiles file
# Read station data 

library(foreign) # To read .dbf files
fname = paste0("1_input/GHCN/GHCN_CONUS_LongTerm_MeanCDD_", start_year, "to", finish_year, "w_MeanRatioCDD.csv")
station_data = read.csv(fname, header = TRUE)

station_data = station_data %>% filter(CDD_LT_stC > 0)
head(station_data)
station_data_original = station_data

library(xlsx)

impsurface_orig = impsurface
elevation_orig = elevation

# Filter stations in the region
# reg_station_data = station_data %>% filter(RASTERVALU > -9999)
reg_station_data = station_data


# CONUS outline
outline = outline_full #[which(outline_full@data$RECS09_Dom == i), ]

## Mask impervious area and elevation based on raster
impsurface = crop(impsurface_orig, outline)
# impsurface = mask(impsurface, outline)

elevation = crop(elevation_orig, outline)
# elevation = mask(elevation, outline)


## RasterStack with elevation and impsurface
impsurface_elevation = raster::stack(impsurface, elevation)
names(impsurface_elevation) = c("Imperv1km", "Elev_new") # Imperv1km + Elev_new

impsurface_mat = as.matrix(impsurface) # Convert to matrix
impsurface_df = as.data.frame(impsurface) # Convert to data frame
ind_impsurface = which(!is.na(impsurface_mat))
sp.studygrid = as(impsurface_elevation, "SpatialGridDataFrame") # Convert to spatial grid

class(sp.studygrid)

## Convert to spatial points
sp.station_data = SpatialPointsDataFrame(data = station_data, coords = cbind(station_data$lon_st, station_data$lat_st),
                                         proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ")))
sp.station_data_prj = spTransform(sp.station_data, CRS(as.character("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")))

sp.station_data_prj$Imperv1km = raster::extract(impsurface, sp.station_data_prj)
sp.station_data_prj$Elev_new = raster::extract(elevation, sp.station_data_prj)
sp.station_data_prj$CDD_LT_NNR_C = raster::extract(cdd_full, sp.station_data_prj)

sp.station_data_prj$NTL_Urb = raster::extract(urb_extent_NTL, sp.station_data_prj)

df.station_data_prj = as.data.frame(sp.station_data_prj)
df.station_data_prj = mutate(df.station_data_prj, diffCDD_LTC = CDD_LT_stC - CDD_LT_NNR_C)

df.station_data_prj$compare_Diffs_LTC = ifelse(df.station_data_prj$diffCDD_LTC > 0, "Yes", "No")

df.station_data_prj = df.station_data_prj %>% filter(!is.na(compare_Diffs_LTC))

# Urban stations
df.station_data_prj_urb = df.station_data_prj %>% filter(!is.na(NTL_Urb))

ggplot(data = df.station_data_prj_urb, mapping = aes(x = compare_Diffs_LTC)) + geom_bar(fill = "steelblue") + ylim(0, 1200) +
  labs(x = "GHCN > NNR ?", y = "Count", title = "Difference between Long-Term Mean CDDs",
       subtitle = "GHCN minus NNR in Urban Areas") + theme_bw(base_size = 18)


## Split stations
set.seed(1234)
ind_sample = sample(1:nrow(df.station_data_prj), nrow(df.station_data_prj)*0.8, replace=F)
cal_data = df.station_data_prj[ind_sample, ]
val_data = df.station_data_prj[-ind_sample, ]

cal_data = cal_data[complete.cases(cal_data[,14]), ]

## Filter NAs in cal_data
cal_data = cal_data %>% filter(!is.na(Elev_new))

sp.cal_data = SpatialPointsDataFrame(data = cal_data, coords = cbind(cal_data$lon_st, cal_data$lat_st),
                                     proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ")))
sp.cal_data = spTransform(sp.cal_data, CRS(as.character("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")))

sp.val_data = SpatialPointsDataFrame(data = val_data, coords = cbind(val_data$lon_st, val_data$lat_st),
                                     proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ")))
sp.val_data = spTransform(sp.val_data, CRS(as.character("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")))


## Convert the impervious layer to spatial points
impsurface_matrix = rasterToPoints(impsurface)
df_impsurface_matrix = data.frame(impsurface_matrix)
colnames(df_impsurface_matrix)[3] = "Imperv"
sub_df_impusrface_matrix = df_impsurface_matrix[df_impsurface_matrix$Imperv >= 3, ] # Create subset

## Calculate variogram
var.resid = variogram(CDD_LT_stC ~ Imperv1km + Elev_new, data = sp.cal_data) # Additional arguments: cutoff = , width = 
plot(var.resid, plot.numbers = T, pch = "+")


## Using "gstat" package And local kriging
Modeltype = "Gau" #"Exc" "Sph" "Exp" "Gau" "Pen" "Lin"
mod.resid = vgm(Modeltype) #Modeltype = "Exc" #"Exc" "Sph" "Exp" "Gau" "Pen" "Lin"
#vgm(0.022,Modeltype,70000,0.035) #vgm(0.08,"Sph",800,0.03)---> psill model, model, range model, psill nugget
plot(var.resid, pl = T, model = mod.resid)

## Fit model parameters
modelfit = fit.variogram(var.resid, mod.resid)
plot(var.resid, pl = T, model = modelfit, main = Modeltype)

## RUN UK
krig.dd_all = gstat(id = "UK_fit", formula = CDD_LT_stC ~ Imperv1km + Elev_new,
                    data = sp.cal_data, nmax = 200) #, model = modelfit
summary(krig.dd_all)

UK_predictions = predict(krig.dd_all, sp.studygrid, debug.level = 0)
ras_UK_predictions = raster(UK_predictions) # Convert to raster
ras_UK_predictions[ras_UK_predictions < 0] = 0

#Save Kriging results
# NOTE: The path after Data_Repository should be unchanged. Change the path before Data_Repository with your drive letter.
if (reanalysis == 'R1') {
  fpath = paste0("2_output/UK_Script/UK_predictions_air2m_R1_meanCDD_C_", start_year,"_to_", finish_year, ".Rda") # "_", Modeltype,
} else if (reanalysis == 'R2') {
  fpath = paste0("2_output/UK_Script/UK_predictions_air2m_R2_meanCDD_C_", start_year,"_to_", finish_year, ".Rda") # "_", Modeltype,
}

save(krig.dd_all, file = fpath) # Need to do once
# load(fpath)
 
pred_CDD = raster::extract(ras_UK_predictions, sp.val_data) # Predicted CDD at the validation points

obs_CDD = val_data$CDD_st_C
 
df.pred_val = data.frame(obs_CDD, pred_CDD)

## Save as tif file
if (reanalysis == 'R1') {
  fpath = paste0("2_output/UK_Script/CONUS_UK_pred_air2m_R1_meanCDD_pred_", start_year,"_to_", finish_year, "_degC.tif")
  writeRaster(ras_UK_predictions, filename = fpath, format = "GTiff", overwrite = TRUE)
  
  fpath = paste0("2_output/UK_Script/CONUS_UK_air2m_R1_meanCDD_pred_", start_year,"_to_", finish_year, "_degC.Rda")
  save(ras_UK_predictions, file = fpath) # Need to do once
  load(fpath)
} else if (reanalysis == 'R2') {
  fpath = paste0("2_output/UK_Script/CONUS_RK_pred_air2m_R2_meanCDD_Diff_2000_2009_degC.tif")
  writeRaster(ras_UK_predictions, filename = fpath, format = "GTiff", overwrite = TRUE)
  
  fpath = paste0("2_output/UK_Script/CONUS_RK_air2m_R2_meanCDD_Diff_2000_2009_degC.Rda")
  save(ras_UK_predictions, file = fpath) # Need to do once
  load(fpath)
}


## in Deg F
df.pred_val$obs_CDD_F <- df.pred_val$obs_CDD * 1.8
df.pred_val$pred_CDD_F <- df.pred_val$pred_CDD * 1.8

bias_RK_F = bias(df.pred_val$pred_CDD_F, df.pred_val$obs_CDD_F)
rmse_RK_F = rmse(df.pred_val$pred_CDD_F, df.pred_val$obs_CDD_F)
sesy_RK_F = rmse_RK_F / sd(obs_CDD * 1.8)
val_stat_F = summary(lm(pred_CDD_F~obs_CDD_F, data = df.pred_val))
qplot(obs_CDD_F, pred_CDD_F, data = df.pred_val, geom = "point", alpha=I(0.8), size = I(2), color = I("grey38")) + # color = I(df.pred_val$colPalette)
  # scale_color_gradient(low = "blue", high = "red", limits = c(0,100)) +#, guide = guide_legend(title.position = "top", nrow = 1)) +
  theme_light(base_size = 19) + xlim(0, 5000) + ylim(0, 5000) + 
  # ggtitle('CONUS') + #ggtitle(paste0(domain, " (", Modeltype, ")")) +
  # theme(legend.text = element_text(size = 20), legend.position = "bottom", legend.direction = "horizontal", legend.key.width=unit(4, "line")) +
  labs(x = 'Observed CDD', y = 'Predicted CDD from UK', color = "Impverviousness [%]    ") +
  annotate("text", x = 750, y = 4000, label = paste0(#"n = ", length(df.pred_val$obs_CDD_F), 
                                                     "\nBias = ", round(bias_RK_F, 1), 
                                                     "\nRMSE = ", round(rmse_RK_F, 1),
                                                     "\n RMSE_r = ", round(sesy_RK_F, 1),
                                                     "\nR2 = ", round(val_stat_F$r.squared, 1)), size = 5) + 
  geom_abline() + 
  theme(legend.direction = "vertical", legend.position = "right", legend.box = "vertical")
# image_fname = paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/RK/Images/CensusDivFull/24Nov11/", domain, "_val_RKinR_", Modeltype, ".tif")
imagefname = paste0("2_output/UK_Script/CONUS_val_CDD_UKinR_", start_year,"_to_", finish_year, "_", Modeltype, ".jpg")
ggsave(imagefname, plot = last_plot(), width = 6, height = 6, units = 'in', dpi = 320, device = "jpeg")


cdd_uhi = ras_UK_predictions

## Convert the CDDs to F
cdd_uhi = 1.8 * cdd_uhi # Convert C to F for CDD
cdd_nonuhi = 1.8 * cdd_nonuhi # Convert C to F for CDD


## Extract Predicted Value
pred_CDD_ALL = raster::extract(cdd_uhi, sp.station_data)

sp.station_data$CDD_F_UK = pred_CDD_ALL
sp.station_data$rCDD_F = sp.station_data$CDD_LT_stF 

diffCDD = cdd_uhi - cdd_nonuhi
ratioCDD = cdd_uhi / cdd_nonuhi

if (reanalysis == 'R1') {
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_", start_year,"_to_", finish_year, "_air2m_F_R1_UHI.Rda")
  save(cdd_uhi, file = fpath) # Need to do once
  
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_", start_year,"_to_", finish_year, "_air2m_F_R1_non_UHI.Rda")
  save(cdd_nonuhi, file = fpath) # Need to do once
  
  fpath = paste0("2_output/UK_Script/CONUS_mean_diffCDD_", start_year,"_to_", finish_year, "_air2m_F_R1_UHI.Rda")
  save(diffCDD, file = fpath) # Need to do once
  
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_", start_year,"_to_", finish_year, "_air2m_F_R1_non_UHI.tif")
  writeRaster(cdd_nonuhi, filename = fpath, format = "GTiff", overwrite = TRUE)
  
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_", start_year,"_to_", finish_year, "_air2m_F_R1_UHI.tif")
  writeRaster(cdd_uhi, filename = fpath, format = "GTiff", overwrite = TRUE)
  
  fpath = paste0("2_output/UK_Script/CONUS_mean_diffCDD_", start_year,"_to_", finish_year, "_air2m_F_R1_UHI.tif")
  writeRaster(diffCDD, filename = fpath, format = "GTiff", overwrite = TRUE)

} else if (reanalysis == 'R2') {
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_2000_2009_air2m_F_R2_UHI.Rda")
  save(cdd_uhi, file = fpath) # Need to do once
  
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_2000_2009_air2m_F_R2_non_UHI.Rda")
  save(cdd_nonuhi, file = fpath) # Need to do once
  
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_2000_2009_air2m_F_R2_non_UHI.tif")
  writeRaster(cdd_nonuhi, filename = fpath, format = "GTiff", overwrite = TRUE)
  
  fpath = paste0("2_output/UK_Script/CONUS_meanCDD_2000_2009_air2m_F_R2_UHI.tif")
  writeRaster(cdd_uhi, filename = fpath, format = "GTiff", overwrite = TRUE)
  
  fpath = paste0("2_output/UK_Script/CONUS_mean_diffCDD_2000_2009_air2m_F_R2_UHI.tif")
  writeRaster(diffCDD, filename = fpath, format = "GTiff", overwrite = TRUE)
  
  fpath = paste0("2_output/UK_Script/CONUS_mean_ratioCDD_2000_2009_air2m_F_R2_UHI.tif")
  writeRaster(ratioCDD, filename = fpath, format = "GTiff", overwrite = TRUE)
}

## Save all variables
if (reanalysis == 'R1') {
  save.image(file = paste0("2_output/UK_Script/CONUS_R1_meanCDD_", start_year,"_to_", finish_year, "_air2m_allVars_RK_estimation.RData"))
} else if (reanalysis == 'R2') {
  save.image(file = paste0("2_output/UK_Script/CONUS_R2_meanCDD_2000_2009_air2m_allVars_RK_estimation.RData"))
}

plot(diffCDD > 0)

## Display computational run-time
toc = proc.time()
elapsed_time = toc - tic
cat("\n Total elapsed time: ", elapsed_time[3]," seconds/", elapsed_time[3]/60," minutes/", elapsed_time[3]/3600, "hours")

