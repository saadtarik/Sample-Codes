# Script to investigate the urban cluster data and CDD difference
#
# Author: Saad Tarik
# Created: Sep 21, 2018

## Run this script after running EnergyuseCooling(Heating)_Estimation_From_CBECS_RECS.R

#This version is used for cooling only - Based on EnergyuseTotal_from_HeatingCooling_CBECS_RECS_w_urb_mask_v2.R 

## Read Urban nightlight shapefile
font_size <- 18

library(raster)
fname <- "1_input/NightLight/US_urban_cluster/urbanCluster_plg_prj_gt500km2.shp"
urb_500km2 <- shapefile(fname)

fname <- "1_input/NightLight/US_urban_cluster/urbanCluster_plg_prj_CONUS_Outline_intersect_SpJoin_CenDiv.shp"
urb_all <- shapefile(fname)

fname <- "1_input/NightLight/US_urban_cluster/urb_conus_lessthan_elev_500m_w_CenDiv.shp"
urb_all_500m <- shapefile(fname)

## Read rasters for Cooling energy use
fname <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/CBECS_CenDiv_CONUS_CoolingElecUse_nonUHI_Full.tif"
CBECS_Cooling_nonUHI <- raster::raster(fname)

fname <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/CBECS_CenDiv_CONUS_CoolingElecUse_UHI_Full.tif"
CBECS_Cooling_UHI <- raster::raster(fname)

fname <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/RECS_CenDiv_CONUS_CoolingElecUse_nonUHI_Full.tif"
RECS_Cooling_nonUHI <- raster::raster(fname)

fname <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Cooling/UK_Estimations/RECS_CenDiv_CONUS_CoolingElecUse_UHI_Full.tif"
RECS_Cooling_UHI <- raster::raster(fname)

# fname <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse/Heating/RECS_CenDiv_HeatingElecUse_nonUHI_Full.tif"
# RECS_Heating_nonUHI <- raster::raster(fname)
# RECS_Heating_nonUHI[RECS_Heating_nonUHI] <- 0 ## dummy raster

# Calculate difference in energy use
CBECS_diff_Cooling_Energy <- CBECS_Cooling_UHI - CBECS_Cooling_nonUHI
RECS_diff_Cooling_Energy <- RECS_Cooling_UHI - RECS_Cooling_nonUHI

## Calculate percent difference rasters
CBECS_percent_diff_Cooling_Energy <- (CBECS_diff_Cooling_Energy / CBECS_Cooling_UHI) * 100
RECS_percent_diff_Cooling_Energy <- (RECS_diff_Cooling_Energy / RECS_Cooling_UHI) * 100

## Mask
CBECS_percent_diff_Cooling_Energy <- raster::mask(CBECS_percent_diff_Cooling_Energy, urb_all_500m)
RECS_percent_diff_Cooling_Energy <- raster::mask(RECS_percent_diff_Cooling_Energy, urb_all_500m)

# Save Rasters
fname_ras <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse_diff_Percent/Cooling/UK_Estimations/CBECS_CONUS_percent_diff_Cooling_Energy.tif"
writeRaster(CBECS_percent_diff_Cooling_Energy, filename = fname_ras, overwrite = TRUE)
fname_ras <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse_diff_Percent/Cooling/UK_Estimations/RECS_CONUS_percent_diff_Cooling_Energy.tif"
writeRaster(RECS_percent_diff_Cooling_Energy, filename = fname_ras, overwrite = TRUE)


## Read rasters for cooling electricity use 
# fname <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse_diff_Percent/CBECS_CenDiv_CoolingElec_Percent_Diff_UHI_URB.tif"
# CBECS_cooling_elec_diff_full <- raster(fname)
# CBECS_cooling_elec_diff_urb_mask <- mask(CBECS_cooling_elec_diff_full, urb_500km2)
# 
# fname <- "2_output/Energyuse_Script/NNR_RAW/Rasters_CBEC_RECS_Energyuse_diff_Percent/RECS_CenDiv_CoolingElec_Percent_Diff_UHI_URB.tif"
# RECS_cooling_elec_diff_full <- raster(fname)
# RECS_cooling_elec_diff_urb_mask <- mask(RECS_cooling_elec_diff_full, urb_500km2)

fname <- "1_input/GRUMP_Population/gpw-v4-population-count_2010_CONUS_Prj.tif"
population_full <- raster(fname)
population_mask <- mask(population_full, urb_500km2)
population_mask_urb_all <- mask(population_full, urb_all)

fname <- "1_input/DEM/SRTM_1km_CONUS_final.tif"
elev_full <- raster(fname)
elev_mask <- mask(elev_full, urb_500km2)
elev_mask_urb_all <- mask(elev_full, urb_all)

## Extract df from shapefile
df.urban500km2 <- urb_500km2@data
df.urban_all <- urb_all@data
df.urban_all_500m <- urb_all_500m@data

library(dplyr)
df.urban500km2 <- df.urban500km2 %>% mutate(ID = seq(1, nrow(df.urban500km2), by = 1))
df.urban_all <- df.urban_all %>% mutate(ID = seq(1, nrow(df.urban_all), by = 1))
df.urban_all_500m <- df.urban_all_500m %>% mutate(ID = seq(1, nrow(df.urban_all_500m), by = 1))

## Zonal Statistics (NOTE: here *_mean_* indicates sum in each urban area) for urban areas lower than 500 m

CBECS_mean_Cooling_UHI <- raster::extract(CBECS_Cooling_UHI, urb_all_500m, fun = sum, na.rm = TRUE, df = TRUE)
CBECS_mean_Cooling_nonUHI <- raster::extract(CBECS_Cooling_nonUHI, urb_all_500m, fun = sum, na.rm = TRUE, df = TRUE)
RECS_mean_Cooling_UHI <- raster::extract(RECS_Cooling_UHI, urb_all_500m, fun = sum, na.rm = TRUE, df = TRUE)
RECS_mean_Cooling_nonUHI <- raster::extract(RECS_Cooling_nonUHI, urb_all_500m, fun = sum, na.rm = TRUE, df = TRUE)


# CBECS_diff_Total_Energy_mean <- raster::extract(CBECS_diff_Total_Energy, urb_all_500m, fun = mean, na.rm = TRUE, df = TRUE)
# RECS_diff_Total_Energy_mean <- raster::extract(RECS_diff_Total_Energy, urb_all_500m, fun = mean, na.rm = TRUE, df = TRUE)

#*** (For all urban areas)

#### Total Energy

df.cbecs_cooling_elec <- merge(CBECS_mean_Cooling_UHI, CBECS_mean_Cooling_nonUHI, by = "ID") #
colnames(df.cbecs_cooling_elec)[2:3] <- c("UHI", "Non-UHI")
# df.cbecs_cooling_elec <- merge(df.cbecs_cooling_elec, df.urban_all_500m, by = "ID")


df.recs_cooling_elec <- merge(RECS_mean_Cooling_UHI, RECS_mean_Cooling_nonUHI, by = "ID") #
colnames(df.recs_cooling_elec)[2:3] <- c("UHI", "Non-UHI")
# df.cbecs_cooling_elec <- merge(df.cbecs_cooling_elec, df.urban_all_500m, by = "ID")

merge_CBECS_sum_Cooling_Energy <- merge(df.cbecs_cooling_elec, urb_all_500m, by = "ID")

merge_RECS_sum_Cooling_Energy <- merge(df.recs_cooling_elec, urb_all_500m, by = "ID")


## Final dataframes for summary

merge_CBECS_sum_Cooling_Energy <- merge_CBECS_sum_Cooling_Energy %>% mutate(Difference = UHI - `Non-UHI`)

merge_RECS_sum_Cooling_Energy <- merge_RECS_sum_Cooling_Energy %>% mutate(Difference = UHI - `Non-UHI`)

## Summarize dataframes for visualization

summary_CBECS_sum_Cooling_Energy <- merge_CBECS_sum_Cooling_Energy %>% group_by(SUB_REGION) %>% 
  dplyr::summarize(UHI = sum(UHI, na.rm = T), 
            NonUHI = sum(`Non-UHI`, na.rm = T),
            Difference = sum(Difference, na.rm = T))

summary_RECS_sum_Cooling_Energy <- merge_RECS_sum_Cooling_Energy %>% group_by(SUB_REGION) %>% 
  dplyr::summarize(UHI = sum(UHI, na.rm = T), 
            NonUHI = sum(`Non-UHI`, na.rm = T),
            Difference = sum(Difference, na.rm = T))

## Visualize
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)
library(gridExtra)

# Total Cooling Consumption
p1 <- summary_CBECS_sum_Cooling_Energy %>% reshape2::melt(id = "SUB_REGION") %>% ggplot(aes(x = SUB_REGION, y = value, fill = variable))+ 
  geom_bar(stat = "identity", position = position_dodge())+
  labs(x = "", y = "Cooling Energy Use \n[thou. BTU]", title = "Commercial Buildings") +
  scale_x_discrete(breaks = unique(summary_CBECS_sum_Cooling_Energy$SUB_REGION), 
                   labels = addline_format(summary_CBECS_sum_Cooling_Energy$SUB_REGION)) +
  guides(fill = guide_legend(title = "")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = font_size - 2),
                     axis.text.y = element_text(size = font_size),
                     axis.title = element_text(size = font_size),
                     plot.title = element_text(size = font_size),
                     legend.position = "bottom", legend.direction = "horizontal",
                     legend.text = element_text(size = font_size))

p2 <- summary_RECS_sum_Cooling_Energy %>% reshape2::melt(id = "SUB_REGION") %>% ggplot(aes(x = SUB_REGION, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = position_dodge())+
  labs(x = "", y = "Cooling Energy Use \n[thou. BTU]", title = "Residential Buildings") +
  scale_x_discrete(breaks = unique(summary_CBECS_sum_Cooling_Energy$SUB_REGION), 
                   labels = addline_format(summary_CBECS_sum_Cooling_Energy$SUB_REGION)) +
  guides(fill = guide_legend(title = "")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = font_size - 2),
                     axis.text.y = element_text(size = font_size),
                     axis.title = element_text(size = font_size),
                     plot.title = element_text(size = font_size),
                     legend.position = "bottom", legend.direction = "horizontal",
                     legend.text = element_text(size = font_size))

# grid.arrange(p1, p2, nrow = 2)
imagefname <- paste0("2_output/Energyuse_Script/NNR_RAW/Paper_Images/CoolingEnergy_and_Diff_CenDiv500m.jpg")
jpeg(filename = imagefname, width = 10, height = 8, units = "in", bg = "white", res = 320)
grid.arrange(p1, p2, nrow = 2)
dev.off()

summary_CBECS_sum_Cooling_Energy %>% mutate(PercentDiff = (Difference / UHI) * 100) %>% 
  write.csv(paste0("2_output/Energyuse_Script/NNR_RAW/Paper_Images/CoolingEnergy_CBECS_and_Diff_CenDiv500m.csv"))
summary_RECS_sum_Cooling_Energy %>% mutate(PercentDiff = (Difference / UHI) * 100) %>% 
  write.csv(paste0("2_output/Energyuse_Script/NNR_RAW/Paper_Images/CoolingEnergy_RECS_and_Diff_CenDiv500m.csv"))

## Percent differences

# Total Cooling Difference
p1 <- summary_CBECS_sum_Cooling_Energy %>% mutate(PercentDiff = (Difference / UHI) * 100) %>% 
  select(SUB_REGION, PercentDiff) %>% 
  ggplot(aes(x = SUB_REGION, y = PercentDiff)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "", y = "Difference in \nCooling Energy Consumption [%]", title = "Commercial Buildings") +
  scale_x_discrete(breaks = unique(summary_CBECS_sum_Cooling_Energy$SUB_REGION), 
                   labels = addline_format(summary_CBECS_sum_Cooling_Energy$SUB_REGION)) +
  guides(fill = guide_legend(title = "")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = font_size - 2),
                     axis.text.y = element_text(size = font_size),
                     axis.title = element_text(size = font_size - 2),
                     plot.title = element_text(size = font_size),
                     legend.position = "bottom", legend.direction = "horizontal",
                     legend.text = element_text(size = font_size))

p2 <- summary_RECS_sum_Cooling_Energy %>% mutate(PercentDiff = (Difference / UHI) * 100) %>% 
  select(SUB_REGION, PercentDiff) %>% 
  ggplot(aes(x = SUB_REGION, y = PercentDiff)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "", y = "Difference in \nCooling Energy Consumption [%]", title = "Residential Buildings") +
  scale_x_discrete(breaks = unique(summary_RECS_sum_Cooling_Energy$SUB_REGION), 
                   labels = addline_format(summary_RECS_sum_Cooling_Energy$SUB_REGION)) +
  guides(fill = guide_legend(title = "")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = font_size - 2),
                     axis.text.y = element_text(size = font_size),
                     axis.title = element_text(size = font_size - 2),
                     plot.title = element_text(size = font_size),
                     legend.position = "bottom", legend.direction = "horizontal",
                     legend.text = element_text(size = font_size))

# grid.arrange(p1, p2, nrow = 2)
imagefname <- paste0("2_output/Energyuse_Script/NNR_RAW/Paper_Images/CoolingEnergy_Percent_Diff_CenDiv500m.jpg")
jpeg(filename = imagefname, width = 10, height = 8, units = "in", bg = "white", res = 320)
grid.arrange(p1, p2, nrow = 2)
dev.off()

save.image(file = "EnergyuseTotal_Cooling_CBECS_RECS_w_urb_mask.RData")

#### DO NOT NEED TO RUN BEYOND THIS POINT #######

df.cbecs_total_elec <- merge(merge_CBECS_mean_Total_Energy_UHI, merge_CBECS_mean_Total_Energy_nonUHI, 
                               c("ID", "CenDiv", "SUB_REGION", "Join_Count", "TARGET_FID", "GRIDCODE")) #
df.cbecs_total_elec <- df.cbecs_total_elec %>% select(ID, DN, Area_km2, Area_m2, CenDiv, SUB_REGION, MeanElecUse.x, MeanElecUse.y)
colnames(df.cbecs_total_elec)[7:8] <- c("UHI", "Non-UHI")

merge_RECS_mean_Total_Energy_UHI <- merge(RECS_mean_Total_Energy_UHI, df.urban_all, by = "ID")
colnames(merge_RECS_mean_Total_Energy_UHI)[2] <- "MeanElecUse"
merge_RECS_mean_Total_Energy_nonUHI <- merge(RECS_mean_Total_Energy_nonUHI, df.urban_all, by = "ID")
colnames(merge_RECS_mean_Total_Energy_nonUHI)[2] <- "MeanElecUse"

df.recs_total_elec <- merge(merge_RECS_mean_Total_Energy_UHI, merge_RECS_mean_Total_Energy_nonUHI, 
                             c("ID", "DN", "Area_km2", "Area_m2", "CenDiv", "SUB_REGION")) #
df.recs_total_elec <- df.recs_total_elec %>% select(ID, DN, Area_km2, Area_m2, CenDiv, SUB_REGION, MeanElecUse.x, MeanElecUse.y)
colnames(df.recs_total_elec)[7:8] <- c("UHI", "Non-UHI")

###*** Filter urban areas lower than 500m

merge_elev_mean_all_500m <- merge_elev_mean_all %>% filter(MeanElev <= 500)
merge_elev_mean_all_500m <- merge_elev_mean_all_500m %>% select(ID, DN, Area_km2, Area_m2, CenDiv, SUB_REGION, MeanElev)

df.cbecs_total_elec_500m <- merge(df.cbecs_total_elec, merge_elev_mean_all_500m)
df.recs_total_elec_500m <- merge(df.recs_total_elec, merge_elev_mean_all_500m)

summary.cbecs_total_elec_500m <- df.cbecs_total_elec_500m %>% group_by(SUB_REGION) %>% summarize(UHI = sum(UHI, na.rm = T), NonUHI = sum(`Non-UHI`, na.rm = T))
summary.recs_total_elec_500m <- df.recs_total_elec_500m %>% group_by(SUB_REGION) %>% summarize(UHI = sum(UHI, na.rm = T), NonUHI = sum(`Non-UHI`, na.rm = T))

sumdata.cbecs_total_elec_500m <- data.frame(value = apply(summary.cbecs_total_elec_500m[2:3], 2, sum))
sumdata.recs_total_elec_500m <- data.frame(value = apply(summary.recs_total_elec_500m[2:3], 2, sum))

sumdata.cbecs_total_elec_500m$type <- rownames(sumdata.cbecs_total_elec_500m)
sumdata.recs_total_elec_500m$type <- rownames(sumdata.recs_total_elec_500m)

ggplot(data = sumdata.cbecs_total_elec_500m, aes(x = type, y = value)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Total Energy Use \n[thou. BTU]", title = "Commercial Buildings") +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = 15),
                     axis.text.y = element_text(size = 15),
                     axis.title = element_text(size = 15))

ggplot(data = sumdata.recs_total_elec_500m, aes(x = type, y = value)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Total Energy Use \n[thou. BTU]", title = "Residential Buildings") +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = 15),
                     axis.text.y = element_text(size = 15),
                     axis.title = element_text(size = 15))

# Median total consumption by Census Divisions
summary.cbecs_total_elec_500m %>% reshape2::melt(id = "SUB_REGION") %>% 
  ggplot(aes(x = SUB_REGION, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "", y = "Total Energy Use \n[thou. BTU]", title = "Commercial Buildings") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, size = 15),
                     axis.text.y = element_text(size = 15),
                     axis.title = element_text(size = 15))

summary.recs_total_elec_500m %>% reshape2::melt(id = "SUB_REGION") %>% 
  ggplot(aes(x = SUB_REGION, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "", y = "Total Energy Use \n[thou. BTU]", title = "Residential Buildings") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, size = 15),
                     axis.text.y = element_text(size = 15),
                     axis.title = element_text(size = 15))

# Percent median total energy use

# gsub(" ", "\n", levels(birds$effect))
p <- summary.cbecs_total_elec_500m %>% mutate(Diff = ((UHI - NonUHI)/UHI) * 100) %>% 
  select(SUB_REGION, Diff) %>% 
  ggplot(aes(x = SUB_REGION, y = Diff)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "", y = "Difference in \nTotal Energy Consumption [%]", title = "Commercial Buildings") +
  scale_x_discrete(breaks = unique(summary.cbecs_total_elec_500m$SUB_REGION), 
                   labels = addline_format(summary.cbecs_total_elec_500m$SUB_REGION)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = font_size),
                     axis.text.y = element_text(size = font_size),
                     axis.title = element_text(size = font_size),
                     plot.title = element_text(size = font_size))
imagefname <- paste0("2_output/Energyuse_Script/NNR_RAW/Paper_Images/CBECS_TotalPercentDiff_TotalEnergy_CenDiv500m.jpg")
jpeg(filename = imagefname, width = 10, height = 6, units = "in", bg = "white", res = 320)
print(p)
dev.off()


p <- summary.recs_total_elec_500m %>% mutate(Diff = ((UHI - NonUHI)/UHI) * 100) %>% 
  select(SUB_REGION, Diff) %>% 
  ggplot(aes(x = SUB_REGION, y = Diff)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "", y = "Difference in \nTotal Energy Consumption [%]", title = "Residential Buildings") +
  scale_x_discrete(breaks = unique(summary.cbecs_total_elec_500m$SUB_REGION), 
                   labels = addline_format(summary.cbecs_total_elec_500m$SUB_REGION)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, size = font_size),
                     axis.text.y = element_text(size = font_size),
                     axis.title = element_text(size = font_size),
                     plot.title = element_text(size = font_size))
imagefname <- paste0("2_output/Energyuse_Script/NNR_RAW/Paper_Images/RECS_TotalPercentDiff_TotalEnergy_CenDiv500m.jpg")
jpeg(filename = imagefname, width = 10, height = 6, units = "in", bg = "white", res = 320)
print(p)
dev.off()
