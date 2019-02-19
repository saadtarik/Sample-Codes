# Sort station data for multiple years (get the yearly number of records, quality index, CDDs, and so on)
#
# Author: Saad Tarik
# Created: Sep 21, 2017

# V2.0 - includes the files from polynomial fitting estimations

## NOTE: Run this script after running the "Process_GHCN_Data_For_Temp_Timeseries.R" script that produces yearly csv files with CDDs

setwd("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/Polynomial_Reg_CDD")

rm(list = ls())

processing_year <- seq(from = 2008, to = 2017, by = 1)
library(xlsx)
library(dplyr)
library(ggplot2)
library(lubridate)
df_yearData <- c()
nfiles <- list.files(pattern = "CDD_stations_w_polyreg_") # Find the .csv files in the directory
nfiles <- nfiles[26:35] ############## CHANGE based on processing years #####################
fullData <- lapply(nfiles, read.csv) # Concatenate into a list

# Concatenate into one data frame
allData <- do.call("rbind", fullData) # Convert the list items to dataframe by rows
head(allData)

## Find the data for all stations
unique_station_ID <- as.character(unique(allData$St_name))

# ind_bad_stations <- match(c("USC00303033", "USC00048655", "USR0000MKIL", "USR0000COAM", "USR0000INFK", "USW00013750", "USC00138808", "USC00093028", "USR0000CBLR"), unique_station_ID)
# unique_station_ID <- unique_station_ID[-ind_bad_stations]

# Initialize data
St_ID <- rep(NA, length(unique_station_ID))
lat_deg <- rep(NA, length(unique_station_ID))
lon_deg <- rep(NA, length(unique_station_ID))
elev_m <- rep(NA, length(unique_station_ID))
mean_records <- rep(NA, length(unique_station_ID))
mean_qual_ind <- rep(NA, length(unique_station_ID))
CDD_C_mean <- rep(NA, length(unique_station_ID))
CDD_F_mean <- rep(NA, length(unique_station_ID))
years_record <- rep(NA, length(unique_station_ID))
CDD_F_median <- rep(NA, length(unique_station_ID))
CDD_F_highlow <- rep(NA, length(unique_station_ID))

CDD_threshold <- 1.5
CDD_upper_thres <- rep(NA, length(unique_station_ID))
CDD_lower_thres <- rep(NA, length(unique_station_ID))

xind <- seq(1, length(processing_year), by = 1)

df_QC_station <- data.frame()

stations_na <- 0 # Initialize the # of stations with NA values of CDD in a year
ind_na <- c() # Initialize the indices of the stations with NA values of CDD in a year

for (i in 1:length(unique_station_ID)) { # length(unique_station_ID)
  cat("\nProcessing station #", i, "of", length(unique_station_ID))
  selected_station_data <- allData %>% filter(St_name == unique_station_ID[i])
  
  if(any(is.na(selected_station_data$CDD_st_F))) {
    stations_na <- stations_na + 1
    ind_na[stations_na] <- i
    next
  }
  
  ## Debug step ## USC00265400 USC00042319 USC00049035 USR0000WDIR USC00303033
  # selected_station_data <- allData %>% filter(St_name == "USR0000COAM")
  
  ## Debug step ##
  # selected_station_data$ind <- xind
  # selected_station_data$yyyy <- format(strptime(processing_year, "%Y"), "%Y")
  
  # ## Debug Plot ##
  # ggplot(data = selected_station_data, mapping = aes(x = yyyy)) + geom_point(aes(y = CDD_st_C), size = 2) +
  #   # geom_point(aes(y = qual_index), size = 2) + scale_y_continuous(sec.axis = sec_axis(~./150, name = "Quality Index"))
  #   ggtitle(unique_station_ID[i]) + labs(x = "Year", y = "CDD [deg C]") + theme_gray(base_size = 14) +
  #   geom_line(stat = "hline", yintercept = "mean")
  # ## Debug Plot ##
  
  St_ID[i] <- unique_station_ID[i]
  lat_deg[i] <- selected_station_data$lat_st[1]
  lon_deg[i] <- selected_station_data$lon_st[1]
  elev_m[i] <- selected_station_data$Elev_st[1]
  years_record[i] <- nrow(selected_station_data)
  mean_records[i] <- mean(selected_station_data$nrecords)
  mean_qual_ind[i] <- mean(selected_station_data$qual_index)
  CDD_C_mean[i] <- mean(selected_station_data$CDD_st_C, na.rm = T)
  CDD_F_mean[i] <- mean(selected_station_data$CDD_st_F, na.rm = T)
  CDD_F_median[i] <- median(selected_station_data$CDD_st_F, na.rm = T)
  CDD_upper_thres[i] <- CDD_F_median[i] + CDD_threshold * CDD_F_median[i]
  CDD_lower_thres[i] <- CDD_F_median[i] - CDD_threshold * CDD_F_median[i]

  selected_station_data$yrecords <- nrow(selected_station_data)
  # selected_station_data$CDD_F_lower <- CDD_lower_thres[i]
  # selected_station_data$CDD_F_median <- CDD_F_median[i]
  # selected_station_data$CDD_F_upper <- CDD_upper_thres[i]
  
  # if (any(selected_station_data$CDD_st_F > selected_station_data$CDD_F_upper | selected_station_data$CDD_st_F < selected_station_data$CDD_F_lower)) {
  #   CDD_F_highlow[i] <- "Yes"
  # } else {
  #   CDD_F_highlow[i] <- "No"
  # }
  
  # selected_station_data$CDD_outside_threshold <- CDD_F_highlow[i]
  df_QC_station <- rbind(df_QC_station, selected_station_data)
}

df_mean_CDD <- data.frame(St_ID, lat_deg, lon_deg, elev_m, years_record, mean_records, mean_qual_ind, CDD_C_mean, CDD_F_mean)
head(df_mean_CDD)

qual_ind_thres <- 0.70
stations_na <- unique_station_ID[ind_na] # Stations with NA values

# test_df <- df_mean_CDD %>% filter(years_record == length(processing_year) & mean_qual_ind > qual_ind_thres)
# max(test_df$CDD_F_mean, na.rm = T)
# # Filter data for stations with records for all year and every day, i.e., mean_qual_ind == 1
# filtered_df_mean_CDD <- df_mean_CDD %>% filter(years_record == length(processing_year) & mean_qual_ind > qual_ind_thres)
# write.xlsx(filtered_df_mean_CDD, file = paste0("GHCN_Long-Term_Mean_CDD_polyreg_qualind", qual_ind_thres * 100, "_", processing_year[1],"to",processing_year[length(processing_year)], ".xlsx"))

## Additional filter
# Filter stations with records for full processing year
df_QC_station_original <- df_QC_station # Keep the original data frame
df_QC_station <- df_QC_station %>% filter(yrecords == length(processing_year))

# Find the stations names with qual index less than the threshold
station_low_qualind <- df_QC_station %>% filter(qual_index < qual_ind_thres)

unique_stations_low_qual <- as.character(unique(station_low_qualind$St_name))

df_QC_station_threshold <- df_QC_station[!df_QC_station$St_name %in% unique_stations_low_qual, ] # Stations passing the 10 year and the 70% criteria

# df_QC_station_thresholdpolyfit <- df_QC_station_threshold %>% filter(CDD_outside_threshold == "No") # Stations passing the polyfit criteria
# 
# df_QC_station_thresholdpolyfit <- df_QC_station_thresholdpolyfit[!df_QC_station_thresholdpolyfit$St_name %in% stations_na, ]

# filtered_data_all %>% group_by(Dates) %>% summarise(sd = sd(glucose), meanbg = mean(glucose))

## Summarize by station with means
# mean_QC_station_threshold <- df_QC_station_thresholdpolyfit %>% group_by(St_name) %>% 
#   summarise_each(funs(mean))

mean_QC_station_threshold <- df_QC_station_threshold %>% group_by(St_name) %>% 
  summarise_each(funs(mean))

mean_QC_station_threshold$X <- NULL
# out_fname <- paste0("GHCN_Long-Term_Mean_CDD_polyreg_qualind_with_medRemove_", qual_ind_thres * 100, "_", processing_year[1],"to",processing_year[length(processing_year)], ".xlsx")
# write.xlsx(mean_QC_station_threshold, file = out_fname)

out_fname <- paste0("QC_GHCN_stations_Passing_Qual_", qual_ind_thres * 100, "_", processing_year[1],"to",processing_year[length(processing_year)], ".xlsx")
write.xlsx(mean_QC_station_threshold, file = out_fname)
