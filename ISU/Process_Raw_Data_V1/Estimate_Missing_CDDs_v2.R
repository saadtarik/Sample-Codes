# Script to produce yearly timeseries of GHCN data and estimate missing temperatures
# 
# Author: Saad Tarik
# Created: Oct 25, 2017

# V2.0: uses the ratio of long-term CDD and available CDD to estimate the missing CDDs for ALL stations

setwd("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/RAW_DATA")

rm(list = ls())

processing_year <- seq(from = 2008, to = 2017, by = 1)

library(dplyr)
library(lubridate)
library(raster)
library(ggplot2)

## Specify user-defined variables
latlim <- c(23, 52)
lonlim <- c(-130, -61)

nfiles <- list.files(pattern = "airUS_") # Find the corresponding file(s) that mathces the pattern in the directory
nfiles <- nfiles[26:35] #### CHANGE based on the processing_year
fullData <- lapply(nfiles, read.csv) # Concatenate into a list

## Concatenate into one dataframe
allData <- do.call("rbind", fullData)  # Convert the list items to dataframe by rows
unique_station_ID <- unique(as.character(allData$name))

new_allData <- data.frame(Station_ID = allData$name, Time = allData$time, TM = allData$type, TVAL = allData$value)

## Read the quality controlled station data and match with new_allData
fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/Polynomial_Reg_CDD/QC_GHCN_stations_Passing_Qual_70_", processing_year[1],"to",processing_year[length(processing_year)], ".xlsx")
library(xlsx)
QC_Station_Data <- read.xlsx(fname, sheetIndex = 1)
QC_Station_names <- as.character(QC_Station_Data$St_name)

## Filter new_allData based on the QC_Station_names
filtered_new_allData <- new_allData[new_allData$Station_ID %in% QC_Station_names, ]

## Date objects
start_date_full <- as.Date(paste0(processing_year[1], "0101"), format = "%Y%m%d")
finish_date_full <- as.Date(paste0(processing_year[length(processing_year)], "1231"), format = "%Y%m%d")
all_dates_full <- data.frame(Date = seq.Date(start_date_full, finish_date_full, by = "day"))

## Loop through the stations
df_QC_stations <- data.frame() # Initialize data frame object
filtered_new_allData$Station_ID <- as.character(filtered_new_allData$Station_ID)
library(reshape2)
station_data <- list()

longTerm_Mean_CDD_station <- rep(NA, length(QC_Station_names))
YearCDD_Ratio <- data.frame(matrix(ncol = length(processing_year) + 1, nrow = length(QC_Station_names)))

St_name_vec <- rep(NA, length(QC_Station_names))
lat_vec <- rep(NA, length(QC_Station_names))
lon_vec <- rep(NA, length(QC_Station_names))

for (j in 1:length(QC_Station_names)) { #length(QC_Station_names)
  # station_of_interest <- "USW00003954" # USR0000WDIR USC00111497 USC00417028 USR0000WDIR USC00303033 USC00138808
  station_of_interest <- QC_Station_names[j]
  # station_data <- lapply(filtered_new_allData, subset, Station_ID == station_of_interest)
  station_data <- filtered_new_allData %>% filter(Station_ID == station_of_interest)
  df_station_data <- station_data # do.call("rbind", station_data)
  
  df_station_data$Date <- as.Date(as.character(df_station_data$Time), '%Y%m%d')
  df_station_data$year <- as.numeric(format(df_station_data$Date,'%Y'))
  
  # Sort TMAX and TMIN dataset
  df_station_data_TMAX <- df_station_data %>% filter(TM == "TMAX")
  df_station_data_TMIN <- df_station_data %>% filter(TM == "TMIN")
  
  df_station_data_combined <- merge(df_station_data_TMAX, df_station_data_TMIN, by = "Date", all = TRUE) # Combines the data.frame with TMAX and TMIN
  
  df_station_data_FULL <- merge(df_station_data_combined, all_dates_full, by = "Date", all = TRUE) # Combines the data.frame with dates and NA TMAX values
  
  df_station_data_FULL <- df_station_data_FULL[, -c(4, 7:9, 11)]
  colnames(df_station_data_FULL)[2:6] <- c("Station_ID", "Time", "TMAX", "Year", "TMIN")
  
  df_station_data_FULL$TAVG <- (df_station_data_FULL$TMAX + df_station_data_FULL$TMIN) / 2
  df_station_data_FULL$dailyCDD_raw <- df_station_data_FULL$TAVG - 18 # Includes negative CDDs, i.e, HDDs
  df_station_data_FULL$dailyCDD_actual <- ifelse(!is.na(df_station_data_FULL$dailyCDD_raw) & df_station_data_FULL$dailyCDD_raw < 0, 0, df_station_data_FULL$dailyCDD_raw)
  
  df_station_data_FULL$day <- day(df_station_data_FULL$Date)
  df_station_data_FULL$month <- month(df_station_data_FULL$Date)
  df_station_data_FULL$year <- year(df_station_data_FULL$Date)
  
  df_station_CDD <- df_station_data_FULL[, c(1, 7:12)]
  
  df_dailyCDD_FULL <- dcast(df_station_CDD, month + day ~ year, value.var = "dailyCDD_actual") 
  df_dailyCDD_FULL$DayMonth <- paste(df_dailyCDD_FULL$day, month.abb[df_dailyCDD_FULL$month], sep = "-")
  
  # m.df_dailyCDD_FULL <- df_dailyCDD_FULL[, -c(1,2)]
  # m.df_dailyCDD_FULL <- melt(m.df_dailyCDD_FULL)
  
  df_dailyCDD_FULL$dailyavgCDD <- rowMeans(x = df_dailyCDD_FULL[, 3:12], na.rm = T) # Long-term avg CDD
  
  df_monthly_mean_CDD <- data.frame(df_dailyCDD_FULL$month, df_dailyCDD_FULL$dailyavgCDD)
  colnames(df_monthly_mean_CDD) <- c("month", "LongTermAvgCDD")
  
  df_daily_mean_CDD <- data.frame(df_dailyCDD_FULL$month, df_dailyCDD_FULL$day, df_dailyCDD_FULL$dailyavgCDD)
  colnames(df_daily_mean_CDD) <- c("month", "day", "LongTermAvgCDD")
  
  sum_daily_mean_CDD <- sum(df_daily_mean_CDD$LongTermAvgCDD)
  
  df_longterm_meanCDD <- aggregate(.~month, data = df_monthly_mean_CDD, sum)
  
  yearly_CDD_vec <- rep(NA, length(processing_year))
  
  MonthlyCDD_Ratio <- data.frame(matrix(ncol = length(processing_year) + 1, nrow = 12))
  # Now loop through each year and fill the missing values
  for (i in 1:length(processing_year)) {
    cat("\nProcessing Station #", j, "of", length(QC_Station_names), "in Year", processing_year[i])
    df_yearly_CDD <- df_station_CDD %>% filter(year == processing_year[i])
    
    # Find the locations of NA in df_yearly_CDD
    df_merged <- merge(df_daily_mean_CDD, df_yearly_CDD, by = c("day", "month"))
    
    df_merged_new <- data.frame(day = df_merged$day, month = df_merged$month, year = df_merged$year, Date = df_merged$Date,
                                LongTerm_AvgCDD = df_merged$LongTermAvgCDD, year_CDD = df_merged$dailyCDD_actual)
    
    df_yearCDD_no_NA <- df_merged_new[!is.na(df_merged_new$year_CDD), ]
    
    sum_LongTerm_AvgCDD <- sum(df_yearCDD_no_NA$LongTerm_AvgCDD)
    sum_year_CDD <- sum(df_yearCDD_no_NA$year_CDD)
    
    ratioCDD <- sum_LongTerm_AvgCDD / sum_year_CDD
    
    estimatedCDD <- df_merged_new$LongTerm_AvgCDD / ratioCDD
    
    ## Fill in the missing CDDs
    df_merged_new$full_CDD <- with(df_merged_new, ifelse(is.na(year_CDD) == T, estimatedCDD, year_CDD))
    
    yearly_CDD_vec[i] <- sum(df_merged_new$full_CDD)
    
    monthly_totCDD <- df_merged_new %>% group_by(month) %>% summarise(sum(full_CDD))
    colnames(monthly_totCDD)[2] <- "CDD"
    MonthlyCDD_Ratio[, i+1] <- monthly_totCDD$CDD
    
    YearCDD_Ratio[j, i+1] <- yearly_CDD_vec[i]
    # dailyCDD_year <- df_dailyCDD_FULL[, i+2]
    # mm <- df_dailyCDD_FULL$month
    # df_yearCDD <- data.frame(mm, dailyCDD_year)
    # sum_df_yearCDD <- df_yearCDD %>% group_by(mm) %>% summarise(sum(dailyCDD_year))
    # sum_df_yearCDD <- data.frame(sum_df_yearCDD)
    # 
    # year_df <- data.frame(df_longterm_meanCDD$month, sum_df_yearCDD$sum.dailyCDD_year., df_longterm_meanCDD$LongTermAvgCDD)
    # colnames(year_df) <- c("month", "YearCDD", "LongTermCDD")
    # 
    # sum_longTermCDD <- sum(year_df$LongTermCDD[is.na(year_df$YearCDD) == FALSE])
    # sum_yearCDD <- sum(year_df$YearCDD, na.rm = T)
    # 
    # ratio <- sum_longTermCDD / sum_yearCDD
    # 
    # year_df$estimatedCDD <- ratio * year_df$LongTermCDD # all_station_TMAX$Temp_all <- with(all_station_TMAX, ifelse(is.na(value) == T, polypredtemp, value))
    # 
    # year_df$allnewCDD <- with(year_df, ifelse(is.na(YearCDD) == T, estimatedCDD, YearCDD))
    # 
    # longTerm_CDD_Ratio <- sum_df_yearCDD$sum.dailyCDD_year. / df_longterm_meanCDD$LongTermAvgCDD
    # longTerm_CDD_Ratio[is.na(longTerm_CDD_Ratio)] <- -9999 # Replace NAs with artificially inflated value of -9999
    # 
    # new_monthly_CDD <- longTerm_CDD_Ratio * df_longterm_meanCDD$LongTermAvgCDD
  }
  MonthlyCDD_Ratio$X1 <- monthly_totCDD$month
  colnames(MonthlyCDD_Ratio) <- c("Month", paste0("Y",processing_year))
  write.csv(MonthlyCDD_Ratio, file = paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/MonthlyCDD_all_Stations/GHCN_", station_of_interest, "_monthlyCDDs_est_allStations_",
                                          processing_year[1],"to",processing_year[length(processing_year)], ".csv"))
  
  longTerm_Mean_CDD_station[j] <- mean(yearly_CDD_vec, na.rm = T)
  YearCDD_Ratio[j, 1] <- QC_Station_names[j]
  # df_station_CDD_FULL <- rbind(QC_Station_names[j], yearly_CDD_vec)
  
}

colnames(YearCDD_Ratio)[1] <- "St_name"

df_station_longTermAvgCDD <- data.frame(QC_Station_names, longTerm_Mean_CDD_station)
df_station_longTermAvgCDD$CDD_LT_stF <- df_station_longTermAvgCDD$longTerm_Mean_CDD_station * 1.8 # Convert to F
colnames(df_station_longTermAvgCDD)[1:2] <- c("St_name", "CDD_LT_stC")

df_new_station_longTermAvgCDD <- merge(QC_Station_Data, df_station_longTermAvgCDD, by = "St_name")

out_fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/GHCN_LongTerm_MeanCDD_",
                    processing_year[1],"to",processing_year[length(processing_year)], "w_MeanRatioCDD.xlsx")

write.xlsx(df_new_station_longTermAvgCDD, file = out_fname)

out_fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/GHCN_LongTerm_MeanCDD_",
                    processing_year[1],"to",processing_year[length(processing_year)], "w_MeanRatioCDD.Rda")

save(df_new_station_longTermAvgCDD, file = out_fname)

YearCDD_Ratio_new <- merge(YearCDD_Ratio, df_new_station_longTermAvgCDD, by = "St_name")

## Find the stations within CONUS
df_CONUS_new_station_longTermAvgCDD <- df_new_station_longTermAvgCDD %>% filter(lat_st > latlim[1] & lat_st < latlim[2] & lon_st > lonlim[1] & lon_st < lonlim[2])
df_CONUS_new_station_longTermAvgCDD_noNA <- df_CONUS_new_station_longTermAvgCDD[!is.na(df_CONUS_new_station_longTermAvgCDD$CDD_LT_stC), ]
ggplot(df_CONUS_new_station_longTermAvgCDD_noNA, aes(x = lon_st, y = lat_st)) + geom_point(aes(color = CDD_LT_stC))

out_fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/GHCN_CONUS_LongTerm_MeanCDD_",
                    processing_year[1],"to",processing_year[length(processing_year)], "w_MeanRatioCDD.xlsx")

write.xlsx(df_CONUS_new_station_longTermAvgCDD_noNA, file = out_fname)

out_fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/GHCN_CONUS_LongTerm_MeanCDD_",
                    processing_year[1],"to",processing_year[length(processing_year)], "w_MeanRatioCDD.csv")

write.csv(df_CONUS_new_station_longTermAvgCDD_noNA, file = out_fname)

out_fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/GHCN_CONUS_LongTerm_MeanCDD_",
                    processing_year[1],"to",processing_year[length(processing_year)], "w_MeanRatioCDD.Rda")

save(df_CONUS_new_station_longTermAvgCDD_noNA, file = out_fname)

YearCDD_Ratio_new_CONUS <- YearCDD_Ratio_new %>% filter(lat_st > latlim[1] & lat_st < latlim[2] & lon_st > lonlim[1] & lon_st < lonlim[2])

YearCDD_Ratio_new_CONUS_noNA <- YearCDD_Ratio_new_CONUS[!is.na(YearCDD_Ratio_new_CONUS$CDD_LT_stC), ]

out_fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/GHCN_CONUS_YearCDD_w_LT_Mean_",
                    processing_year[1],"to",processing_year[length(processing_year)], "w_MeanRatioCDD.csv")

write.csv(YearCDD_Ratio_new_CONUS_noNA, file = out_fname)

out_fname <- paste0("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009/Processed_Data/GHCN_CONUS_YearCDD_w_LT_Mean_",
                    processing_year[1],"to",processing_year[length(processing_year)], "w_MeanRatioCDD.Rda")

save(YearCDD_Ratio_new_CONUS_noNA, file = out_fname)

## Find Yearly Anomaly
YearCDD_Ratio_new_CONUS_noNA_ANOMALY <- YearCDD_Ratio_new_CONUS_noNA
for (k in 1:length(processing_year)) {
  YearCDD_Ratio_new_CONUS_noNA_ANOMALY[, k + 1] <- YearCDD_Ratio_new_CONUS_noNA_ANOMALY[, k + 1] - YearCDD_Ratio_new_CONUS_noNA_ANOMALY$CDD_LT_stC
}

colnames(YearCDD_Ratio_new_CONUS_noNA_ANOMALY)[2:11] <- paste0("Y", processing_year[1]:processing_year[length(processing_year)])

m.YearCDD_Ratio_new_CONUS_noNA_ANOMALY <- YearCDD_Ratio_new_CONUS_noNA_ANOMALY %>% 
  dplyr::select(St_name, lat_st, lon_st, Y2008, Y2009, Y2010, Y2011, Y2012, Y2013, Y2014, Y2015, Y2016, Y2017) %>% 
  melt(id = c("St_name", "lat_st", "lon_st"))

library(ggplot2)
library(ggmap)
library(ggrepel)
basemap <- get_map(location = c(lon = -97.51, lat = 39.35), zoom = 4, scale = "auto", maptype = "toner", source = "stamen")
# Humid Continental c(lon = -86.7, lat = 39.6), Semiarid c(lon = -108.0, lat = 39.3),
plotmap <- ggmap(basemap)

plotmap +
  geom_point(aes(x = lon_st, y = lat_st, color = value), data = m.YearCDD_Ratio_new_CONUS_noNA_ANOMALY, size = 1) +
  scale_color_gradient2(low = "darkblue", midpoint = 0, high = "darkred") +
  facet_wrap(~variable, ncol = 5) + ggtitle("") +
  theme_grey(base_size = 15)
