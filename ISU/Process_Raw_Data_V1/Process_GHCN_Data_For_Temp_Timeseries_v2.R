# Script for processing GHCN data from multiple years
# 
# Author: Saad Tarik
# Created: Sep 19, 2017

# V2.0 - includes missing dates and values estimation using polynomial fitting

setwd("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009")

rm(list = ls())

processing_year <- seq(from = 2013, to = 2017, by = 1)

library(dplyr)
library(xlsx)
library(R.matlab)
library(ggplot2)
library(lubridate)

tic <- proc.time()
for (i in 1:length(processing_year)) {
  cat("Processing: ", processing_year[i])
  fname <- paste0("RAW_DATA/airUS_", processing_year[i], ".csv")
  alldata <- read.csv(file = fname, header = TRUE, sep = ",")
  unique_station_ID <- unique(as.character(alldata$name))
  
  newdf <- data.frame(alldata$name, alldata$time, alldata$type, alldata$value)
  colnames(newdf) <- c("Station_All", "Date_vec", "TM", "TVAL")
  
  # Get Station info
  fname <- paste0("RAW_DATA/USstation", processing_year[i], ".csv")
  station_info <- read.csv(file = fname, header = TRUE, sep = ",")
  
  # write.xlsx(newdf, file = paste0("RAW_DATA/US_", processing_year[i], "Temp_data.xlsx"))
  # writeMat(paste0("RAW_DATA/US_", processing_year[i], "Temp_data.mat"), All_station_data = newdf)
  qual_index <- rep(NA, length(unique_station_ID))
  nrecords <- rep(NA, length(unique_station_ID))
  CDD_st_C <- rep(NA, length(unique_station_ID))
  CDD_st_F <- rep(NA, length(unique_station_ID))
  St_name <- rep(NA, length(unique_station_ID))
  lon_st <- rep(NA, length(unique_station_ID))
  lat_st <- rep(NA, length(unique_station_ID))
  Elev_st <- rep(NA, length(unique_station_ID))
  
  for (j in 1:length(unique_station_ID)) {
    
    cat("\n Processing Station:", j, "of", length(unique_station_ID), "in Year", processing_year[i])
    
    current_station_id <- unique_station_ID[j]
    filtered_station_data <- alldata %>% filter(name == unique_station_ID[j])
    filtered_station_info <- station_info %>% filter(name == unique_station_ID[j])
    nrecords[j] <- nrow(filtered_station_data)
    
    St_name[j] <- as.character(filtered_station_info$name)
    lon_st[j] <- filtered_station_info$x
    lat_st[j] <- filtered_station_info$y
    Elev_st[j] <- filtered_station_info$z
    
    # Check if leap year
    if((processing_year[i] %% 4) == 0) {
      if((processing_year[i] %% 100) == 0) {
        if((processing_year[i] %% 400) == 0) {
          leap_year <- "TRUE"
        } else {
          leap_year <- "FALSE"
        }
      } else {
        leap_year <- "TRUE"
      }
    } else {
      leap_year <- "FALSE"
    }
    # Calculate quality index for each station
    if (leap_year == "TRUE") {
      qual_index[j] <- nrow(filtered_station_data) / 732
    } else if (leap_year == "FALSE") {
      qual_index[j] <- nrow(filtered_station_data) / 730
    }
    
    filtered_station_data$timechar <- as.character(filtered_station_data$time)
    # ########### DEBUG Figure ###########
    # ggplot(data = filtered_station_data, mapping = aes(x = as.Date(as.character(time), '%Y%m%d'), y = value, color = type)) +
    #   geom_line(size = 1.5) + labs(x = "Date", y = "Temp. [C]") + ggtitle(paste0(unique_station_ID[j], " (", processing_year[i], ")"))
    # ########### DEBUG Figure ###########
    
    # Calculate CDDs
    station_TMAX <- filtered_station_data %>% filter(type == "TMAX")
    station_TMIN <- filtered_station_data %>% filter(type == "TMIN")
    if (nrow(station_TMAX) != 0 & nrow(station_TMIN) != 0) {  
      ## Estimate missing values
      # Find missing dates and merge with the dataset
      station_TMAX$timestampasdate <- (as.Date(as.character(station_TMAX$time), format = "%Y%m%d"))
      station_TMIN$timestampasdate <- (as.Date(as.character(station_TMIN$time), format = "%Y%m%d"))
      start_date <- as.Date(paste0(processing_year[i], "0101"), format = "%Y%m%d")
      finish_date <- as.Date(paste0(processing_year[i], "1231"), format = "%Y%m%d")
      
      allDates <- data.frame(timestampasdate = seq.Date(start_date, finish_date, by = "day"))
      
      # Merge with original dataframe
      all_station_TMAX <- merge(station_TMAX, allDates, by = "timestampasdate", all = TRUE)
      all_station_TMIN <- merge(station_TMIN, allDates, by = "timestampasdate", all = TRUE)
      
      ## Fit polynomial to estimate missing value
      regobj_TMAX <- lm(all_station_TMAX$value ~ poly(all_station_TMAX$timestampasdate, 5))
      regobj_TMIN <- lm(all_station_TMIN$value ~ poly(all_station_TMIN$timestampasdate, 5))
      all_station_TMAX$polypredtemp = predict(regobj_TMAX, newdata = all_station_TMAX$timestampasdate)
      all_station_TMIN$polypredtemp = predict(regobj_TMIN, newdata = all_station_TMIN$timestampasdate)
      all_station_TMAX$Temp_all = with(all_station_TMAX, ifelse(is.na(value) == T, polypredtemp, value))
      all_station_TMIN$Temp_all = with(all_station_TMIN, ifelse(is.na(value) == T, polypredtemp, value))
      
      all_station_TMAX$type = with(all_station_TMAX, ifelse(is.na(type) == T, "Pred", "TMAX"))
      all_station_TMIN$type = with(all_station_TMIN, ifelse(is.na(type) == T, "Pred", "TMIN"))
      # ########### DEBUG Figure ###########
      # ggplot(data = all_station_TMAX, mapping = aes(x = timestampasdate, y = Temp_all, color = type)) +
      #   geom_line(size = 1.5) + labs(x = "Date", y = "Temp. [C]") + ggtitle(paste0(unique_station_ID[j], " (", processing_year[i], ")"))
      # ggplot(data = all_station_TMIN, mapping = aes(x = timestampasdate, y = Temp_all, color = type)) +
      #   geom_line(size = 1.5) + labs(x = "Date", y = "Temp. [C]") + ggtitle(paste0(unique_station_ID[j], " (", processing_year[i], ")"))
      # ########### DEBUG Figure ###########
      
      df_full_Year_data <- data.frame(all_station_TMIN$timestampasdate, all_station_TMAX$Temp_all, all_station_TMIN$Temp_all)
      colnames(df_full_Year_data) <- c("Dates", "TMAX", "TMIN")
      
      # df_station_TMEAN <- filtered_station_data %>% select(timechar, value)
      # station_TMEAN <- df_station_TMEAN %>% group_by(timechar) %>% summarise(meantemp = mean(value))
      # station_TMEAN$CDD_daily <- station_TMEAN$meantemp - 18
      # station_TMEAN$CDD_daily[station_TMEAN$CDD_daily < 0] = NA
      df_full_Year_data$meantemp <- rowMeans(subset(df_full_Year_data, select = c(TMAX, TMIN)))
      df_full_Year_data$CDD_daily <- df_full_Year_data$meantemp - 18
      df_full_Year_data$CDD_daily[df_full_Year_data$CDD_daily < 0] = NA
      
      CDD_st_C[j] <- sum(df_full_Year_data$CDD_daily, na.rm = T)
      CDD_st_F[j] <- CDD_st_C[j] * 1.8
    }
    
  }
  
  df_CDD_station <- data.frame(St_name, lat_st, lon_st, Elev_st, nrecords, qual_index, CDD_st_C, CDD_st_F)
  
  # Output to spreadsheet
  outfname <- paste0("Processed_Data/Polynomial_Reg_CDD/CDD_stations_w_polyreg_", processing_year[i], ".csv")
  # write.csv(df_CDD_station, outfname)
}

## Display computational run-time
toc <- proc.time()
elapsed_time <- toc - tic
cat("\n Total elapsed time: ", elapsed_time[3]," seconds/", elapsed_time[3]/60," minutes/", elapsed_time[3]/3600, "hours")