load("Raw_Data/landsat.rda")
R <- 5
G <- 4
B <- 3

library(raster)
plotRGB(all_2019_st_crop, r = R, g = G, b = B, stretch = "lin", axes = TRUE,
        main = paste("Landsat composite w/ Bands", R,",", G,", and", B, "in 2019"))
box(col = "white")
training_2019 <- shapefile("Train_2019.shp")
training_1989 <- shapefile("Train_1989.shp")
# 1989
classes_1989 <- rasterize(training_1989, all_1989_st_crop, field = "Class_ID")
names(classes_1989) <- "Classes"
mask_train_1989 <- mask(all_1989_st_crop, classes_1989)
# Combine the classes raster with the masked raster
mask_train_1989_class <- addLayer(mask_train_1989, classes_1989)

# 2019
classes_2019 <- rasterize(training_2019, all_2019_st_crop, field = "Class_ID")
names(classes_2019) <- "Classes"
mask_train_2019 <- mask(all_2019_st_crop, classes_2019)
# Combine the classes raster with the masked raster
mask_train_2019_class <- addLayer(mask_train_2019, classes_2019)
# 1989
mask_train_1989_values <- getValues(mask_train_1989_class) # Get the values 
mask_train_1989_values <- na.omit(mask_train_1989_values) # Remove NAs
df.mask_train_1989_values <- data.frame(mask_train_1989_values) # Convert to data frame
df.mask_train_1989_values$Classes <- as.factor(df.mask_train_1989_values$Classes) # Convert classes to factor type

# 2019
mask_train_2019_values <- getValues(mask_train_2019_class) # Get the values 
mask_train_2019_values <- na.omit(mask_train_2019_values) # Remove NAs
df.mask_train_2019_values <- data.frame(mask_train_2019_values) # Convert to data frame
df.mask_train_2019_values$Classes <- as.factor(df.mask_train_2019_values$Classes) # Convert classes to factor type
library(randomForest)
model_rf_1989 <- randomForest(x = df.mask_train_1989_values[, c(1:6)],
                              y =df.mask_train_1989_values$Classes,
                              importance = TRUE)

colnames(model_rf_1989$confusion) <- c("No River", "River", "Class Error")
rownames(model_rf_1989$confusion) <- c("No River", "River")

model_rf_2019 <- randomForest(x = df.mask_train_2019_values[, c(1:7)],
                              y =df.mask_train_2019_values$Classes,
                              importance = TRUE)

colnames(model_rf_2019$confusion) <- c("No River", "River", "Class Error")
rownames(model_rf_2019$confusion) <- c("No River", "River")

# Prediction
predict_rf_1989 <- predict(all_1989_st_crop, model = model_rf_1989, na.rm = TRUE)
predict_rf_2019 <- predict(all_2019_st_crop, model = model_rf_2019, na.rm = TRUE)

# Visualilze
plot(predict_rf_2019, col = c("red", "blue"))
legend("bottomleft", legend = c("No River", "River"), fill = c("red", "blue"), bg = "white")
# Set No River as NA export as polygon
predict_rf_1989[predict_rf_1989 == 1] <- NA
predict_rf_2019[predict_rf_2019 == 1] <- NA

# Combined Plot
library(quickPlot)
clearPlot()
Plot(predict_rf_1989, col = "red", legend = FALSE, new = TRUE, title = "")
Plot(predict_rf_2019, col = "blue", legend = FALSE, addTo = "predict_rf_1989", title = FALSE)
legend("bottomleft", legend = c("1989", "2019"), fill = c("red", "blue"), bg = "white")
