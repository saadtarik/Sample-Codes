setwd("C:/Users/starik/Documents/Data_Repository/GHCND/GHCND_from_ST/Data/Data_2009")

library(xlsx)
station_of_interest <- "USW00023293"
st_dat <- read.xlsx2(file = paste0(station_of_interest, "_station_temp_allData_w_NCEP.xlsx"), sheetIndex = 1)
# Convert factors to date
st_dat$date = as.Date(as.character(st_dat$time), '%Y%m%d')
st_dat$tavg = as.numeric(as.character(st_dat$tavg))
st_dat$doe = as.numeric(as.character(st_dat$doe))
st_dat$NNR_995 = as.numeric(as.character(st_dat$NNR_995))
st_dat$NNR_2m = as.numeric(as.character(st_dat$NNR_2m))

library(ggplot2)
ggplot(data = st_dat, mapping = aes(x = date)) + 
  geom_line(mapping = aes(y = tavg), color = 'purple', size = 1.5) +
  geom_line(mapping = aes(y = doe), color = 'orange', size = 1.5) +
  # geom_line(mapping = aes(y = NNR_995), color = 'green', size = 1.5) +
  # geom_line(mapping = aes(y = NNR_2m), color = 'red', size = 1.5) +
  labs(x = "Date", y = "Temperature [C]") + ggtitle("Daily Mean Temperature") + theme_bw(base_size = 14)
