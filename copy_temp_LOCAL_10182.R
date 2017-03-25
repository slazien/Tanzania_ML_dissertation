set.seed(1)

library(ggplot2)
library(ggmap)
library(data.table)
library(RCurl)
library(RJSONIO)
library(plyr)
library(glmnet)
library(randomForest)
library(dplyr)
library(gbm)
library(xgboost)
library(zoo)
library(Ckmeans.1d.dp) 
library(colorRamps)
library(caret)
library(RColorBrewer)
library(scales) #needed?

# CLEAN WORKSPACE?
clean_workspace = T
if(clean_workspace) {
      rm(list = ls())
      gc()
}



# Options
OS = "Windows"
user = "ThinkPad"
load_workspace = F



if(OS == "Linux") {
      setwd("~/OneDrive/DISSERTATION")
} else {
      setwd(paste0("C:/Users/", user, "/OneDrive/DISSERTATION"))
}

if(load_workspace) load("exploratory_workspace.RData")

# ========== LOAD DATA ========== #

if(OS == "Linux") {
      wpd <- fread(paste0("~/OneDrive/DISSERTATION/data/wpdx/Water_Point_Data_Exchange_Complete_Dataset.csv"), 
                   header = TRUE, sep = ",", strip.white = TRUE)
} else {
      wpd <- fread(paste0("C:/Users/", user, "/OneDrive/DISSERTATION/data/wpdx/Water_Point_Data_Exchange_Complete_Dataset.csv"), 
                   header = TRUE, sep = ",", strip.white = TRUE)
}

colnames(wpd) <- gsub("#", "", colnames(wpd))
tz <- wpd[country_id == "TZ"]

# Remove some unecessary columns
rownames(tz) <- seq(1:nrow(tz))
tz <- tz[, !c(2, 6, 15, 17, 18, 19, 20, 21, 22, 23, 26, 27), with = FALSE]
tz$groundwater_depth <- as.factor(readRDS("closest_feature.Rds"))
tz$groundwater_productivity <- as.factor(readRDS("groundwater_productivity.Rds"))
tz$groundwater_storage <- as.factor(readRDS("groundwater_storage.Rds"))

# Fix water_source factor levels
tz$water_source <- gsub("Machine Drilled Borehole", "Machine-drilled borehole", tz$water_source)

# Subset to include only relevant point water source types
relevant_types <- c("Shallow well", "Machine-drilled borehole", "Spring", "Hand-drilled tube well")
tz <- tz[tz$water_source %in% relevant_types]
tz <- as.data.frame(tz)
tz <- na.omit(tz)

# SUBSET BY WATER_SOURCE TYPE, one of:
# 1. Shallow well
# 2. Machine-drilled borehole
# 3. Spring
# 4. Hand-drilled tube well
WATER_SOURCE_SEL = relevant_types[c(1,2,3,4)]

if(WATER_SOURCE_SEL == "") {
      # Do nothing
} else if (WATER_SOURCE_SEL %in% relevant_types) {
      tz <- tz[tz$water_source == WATER_SOURCE_SEL,]
} else {
      stop("Water source type not valid")
}

water <- read.csv("water_data.csv", header = T, skip = 3)

functional_by_year <- readRDS("functional_by_year.Rds")
functional_by_year <- functional_by_year[-c(1:10),]
functional_by_year$diff <- functional_by_year$functional - functional_by_year$not_functional
functional_by_year_melt <- melt(functional_by_year)

# ========== DEFINE FUNCTIONS ========== #

# Utility
get_ellipsis_args <- function(...)
{
      dots <- list(...)
      return(dots)
}

# Distance
dist_euclid <- function(point_lat, point_lon, city_coords) {
      dist <- numeric(length(point_lat))
      dist <- 111 * sqrt((point_lat - city_coords[1])^2 + (point_lon - city_coords[2])^2)
      return(dist)
}

to_radians <- function(x) {
      return(x * pi/180)
}

dist_haversine <- function(point_lat, point_lon, city_coords) {
      phi1 <- to_radians(point_lat)
      phi2 <- to_radians(city_coords[1])
      delta_lat <- phi2 - phi1
      lambda1 <- to_radians(point_lon)
      lambda2 <- to_radians(city_coords[2])
      delta_lon <- lambda2 - lambda1
      distance_hav <- 2 * 6371 * asin(sqrt(sin(delta_lat/2)^2 + cos(phi1) * cos(phi2) * sin(delta_lon/2)^2))
      return(distance_hav)
}

dist_euclid2 <- function(errors_lat, errors_lon) {
      return(sqrt((errors_lat^2 + errors_lon^2)*111))
}

# Lat/lon google api
url <- function(address, return.call = "json", sensor = "false") {
      root <- "http://maps.google.com/maps/api/geocode/"
      u <- paste(root, return.call, "?address=", address, "&sensor=", sensor, sep = "")
      return(URLencode(u))
}

geoCode <- function(address, verbose = FALSE) {
      if (verbose) cat(address, "\n")
      doc <- getURL(url(address))
      x <- fromJSON(doc, simplify = FALSE)
      if (x$status == "OK") {
            lat <- x$results[[1]]$geometry$location$lat
            lng <- x$results[[1]]$geometry$location$lng
            location_type <- x$results[[1]]$geometry$location_type
            formatted_address <- x$results[[1]]$formatted_address
            return(c(lat, lng, location_type, formatted_address))
            Sys.sleep(0.5)
      } else {
            return(c(NA, NA, NA, NA))
      }
}

# Data transformation
standardise <- function(x) {
      return((x - mean(x))/sd(x))
}

add_jitter <- function(actual_lat, actual_lon, predicted_lat, predicted_lon, N = seq(0, 1, by = 0.0001)) {
      mse_lat <- numeric()
      mse_lon <- numeric()
      j = 1
      for(i in N) {
            mse_lat[j] <- mean((actual_lat - jitter(predicted_lat, amount = i))^2)
            mse_lon[j] <- mean((actual_lon - jitter(predicted_lon, amount = i))^2)
            j = j + 1
      }
      return(data.frame(mse_lat = mse_lat, mse_lon = mse_lon))
}

errors_df <- function(predicted_lat, predicted_lon, actual_lat, actual_lon) {
      df <- data.frame(lat = (predicted_lat - actual_lat), lon = (predicted_lon - actual_lon))
      return(df)
}

# Model metrics
mae <- function(preds, test) {
      return(mean(abs(preds - test)))
}

mse <- function(preds, test) {
      return(mean((preds - test)^2))
}

rmse <- function(preds, test) {
      return(sqrt(mean((preds-test)^2)))
}

return_metrics <- function(pred.lat, pred.lon, actual.lat, actual.lon) {
      loss_df <- data.frame(MAE_lat = mae(pred.lat, actual.lat),
                            MAE_lon = mae(pred.lon, actual.lon),
                            MSE_lat = mse(pred.lat, actual.lat),
                            MSE_lon = mse(pred.lon, actual.lon),
                            RMSE_lat = rmse(pred.lat, actual.lat),
                            RMSE_lon = rmse(pred.lon, actual.lon))
      return(loss_df)
}

error_cdf <- function(errors_col, x) {
      return(sum(abs(errors_col) <= x)/length(errors_col))
}

# Plot
wssplot <- function(data, nc = 15, seed = 1234) {
      wss <- numeric(nc)
      for (i in 1:nc) {
            set.seed(seed)
            wss[i] <- sum(kmeans(data, centers = i)$withinss)
      }
      tmp_df <- data.frame(number_clust = 1:nc, wss = wss)
      ggplot(tmp_df, aes(x = number_clust, y = wss)) + geom_point() + geom_line() + 
            xlab("Number of clusters") + ylab("Within group sum of squares") + ggtitle("Within group sum of squares as a function of k") + 
            theme_minimal() + theme(plot.title = element_text(face = "bold", size = 20), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + 
            geom_vline(xintercept = 20, col = "Red", lwd = 1, lty = 3)
}

plot_with_jitter <- function(actual_lat, actual_lon, predicted_lat, predicted_lon, jitter_amount = 0.1, model = " NULL MODEL", alpha = 0.5) {
      ggplot() + theme_minimal() +
            geom_point(aes(x = actual_lat, y = actual_lon), color = "darkgreen") + 
            geom_point(aes(x = jitter(predicted_lat, amount = jitter_amount), 
                           y = jitter(predicted_lon, amount = jitter_amount)), color = "darkred", alpha = alpha) +
            xlab("Lat") + ylab("Lon") + ggtitle(paste0("Actual vs predicted lat/lon for jitter (", jitter_amount, ") ", model, " predictions" ))
}

compare_error_dist <- function(..., breaks = 200) {
      args_e <- get_ellipsis_args(...)
      df_n <- numeric(1)
      for(i in 1:length(args_e)) {
            if(is.data.frame(args_e[[i]])) df_n <- df_n + 1
      }
      par(mfrow = c(df_n, 2))
      for(i in 1:df_n) {
            hist(abs(args_e[[i]]$lat), breaks, 
                 xlab = "Absolute error in lat", ylab = "Frequency", main = paste0(args_e[[length(args_e) - df_n + i]], " - lat"))
            hist(abs(args_e[[i]]$lon), breaks, 
                 xlab = "Absolute error in lon", ylab = "Frequency", main = paste0(args_e[[length(args_e) - df_n + i]], " - lon"))
      }
}

error_cdf_plot <- function(errors_lat, errors_lon, model = "model", save_cdf = FALSE, AUC_max_km = 15) {
      dists <- dist_euclid2(errors_lat, errors_lon)
      km_thresh <- seq(0, max(dists), by = 0.001)
      prob <- numeric(length(km_thresh))
      j = 1
      for(i in km_thresh) {
            prob[j] <- error_cdf(dists, i)
            j = j + 1
      }
      km_cdf <- data.frame(km = km_thresh, prob = prob)
      if(AUC_max_km > max(km_cdf$km)) {
            message(paste0("AUC_max_km is greater than maximum error (", max(km_cdf$km) ,") changing to maximum"))
            AUC_max_km <- max(km_cdf$km)
      }
      if(save_cdf == TRUE) {
            tmp_cdf_df <<- km_cdf
            message("CDF data frame saved to 'tmp_cdf_df' object")
      }
      km_cdf <- km_cdf[1:which(km_cdf$km == AUC_max_km),]
      AUC <- sum(diff(km_cdf$km[order(km_cdf$km)]) * rollmean(km_cdf$prob[order(km_cdf$km)], 2))
      return(ggplot(km_cdf, aes(x = km, y = prob)) + 
                   geom_line(lwd = 1, col = "steelblue") + 
                   theme_minimal() + 
                   scale_x_continuous(breaks = seq(0, max(km_thresh), by = 1)) + 
                   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
                   xlab("Error (km)") + ylab("Probability") + ggtitle(paste0("CDF of straight-line error probability for ", model)) +
                   annotate("text", x = AUC_max_km-3, y = 0.85, label = paste0("AUC = ", round(AUC, 2)), size = 8))
}

# Machine learning
grid_search <- function(param1_seq, param2_seq, 
                        param1_name, param2_name, 
                        params_list = list(eta = 0.05, max_depth = 17, 
                                           colsample_bytree = 0.85, min_child_weight = 1, 
                                           subsample = 0.8, colsample_bylevel = 1, lambda = 0.2), 
                        nrounds_cv = 300, nrounds_grid = 300, save_grids = FALSE) {
      
      len1 <- length(param1_seq)
      len2 <- length(param2_seq)
      rmse_grid_lat <- matrix(nrow = len1, ncol = len2)
      rmse_grid_lon <- rmse_grid_lat
      rownames(rmse_grid_lat) <- param1_seq
      rownames(rmse_grid_lon) <- param1_seq
      colnames(rmse_grid_lat) <- param2_seq
      colnames(rmse_grid_lon) <- param2_seq
      
      color_ramp <- blue2red(len1 * len2)
      par(mfrow = c(1,2))
      
      pb <- txtProgressBar(min = 0, max = len1 * len2, style = 3)
      iter = 1
      
      for(i in 1:len1) {
            for(j in 1:len2) {
                  
                  params = params_list
                  
                  params[[which(names(params) == param1_name)]] <- param1_seq[i]
                  params[[which(names(params) == param2_name)]] <- param2_seq[j]
                  
                  # MODEL
                  
                  xgb.lat <- xgboost(data = model.mat[train,], label = target.lat[train], 
                                     nrounds = nrounds_grid, params = params, metrics = "rmse", verbose = 0, objective = "reg:linear")
                  
                  xgb.lon <- xgboost(data = model.mat[train,], label = target.lon[train], 
                                     nrounds = nrounds_grid, params = params, metrics = "rmse", verbose = 0, objective = "reg:linear")
                  
                  # PREDICT
                  
                  xgb.lat.pred <- predict(xgb.lat, model.mat[-train,])
                  xgb.lon.pred <- predict(xgb.lon, model.mat[-train,])
                  
                  rmse_grid_lat[i, j] <- return_metrics(xgb.lat.pred, xgb.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])$RMSE_lat
                  rmse_grid_lon[i, j] <- return_metrics(xgb.lat.pred, xgb.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])$RMSE_lon
                  
                  # VIS
                  
                  image(rmse_grid_lat, main = "Lat", xlab = param1_name, ylab = param2_name, axes = F, col = color_ramp)
                  axis(1, at = seq(0, 1, length.out = len1), labels = param1_seq)
                  axis(2, at = seq(0, 1, length.out = len2), labels = param2_seq)
                  image(rmse_grid_lon, main = "Lon", xlab = param1_name, ylab = param2_name, axes = F, col = color_ramp)
                  axis(1, at = seq(0, 1, length.out = len1), labels = param1_seq)
                  axis(2, at = seq(0, 1, length.out = len2), labels = param2_seq)
                  
                  setTxtProgressBar(pb, iter)
                  iter = iter + 1
            }
      }
      
      grid_sum <- rmse_grid_lat + rmse_grid_lon
      
      min_idx <- arrayInd(which.min(grid_sum), dim(grid_sum))
      
      lowest_sum_rmse <- grid_sum[min_idx[1], min_idx[2]]
      
      print(paste("Lowest summed RMSE: ", lowest_sum_rmse))
      
      best_param1 <- param1_seq[min_idx[,1]]
      best_param2 <- param2_seq[min_idx[,2]]
      
      print(paste("Best", param1_name, ":", best_param1))
      print(paste("Best", param2_name, ":", best_param2))
      
      best_param_list <- params_list
      best_param_list[[param1_name]] <- best_param1
      best_param_list[[param2_name]] <- best_param2
      
      xgb.lat.cv <- xgb.cv(data = model.mat[train,], label = target.lat[train], 
                           nrounds = nrounds_cv, params = best_param_list, metrics = "rmse", verbose = 0, objective = "reg:linear", nfold = 10)
      
      xgb.lon.cv <- xgb.cv(data = model.mat[train,], label = target.lon[train], 
                           nrounds = nrounds_cv, params = best_param_list, metrics = "rmse", verbose = 0, objective = "reg:linear", nfold = 10)
      
      min_rmse_lat <- min(xgb.lat.cv$test.rmse.mean)
      min_rmse_lon <- min(xgb.lon.cv$test.rmse.mean)
      
      min_rmse_ntrees_lat <- which(xgb.lat.cv$test.rmse.mean == min_rmse_lat)
      min_rmse_ntrees_lon <- which(xgb.lon.cv$test.rmse.mean == min_rmse_lon)
      
      plot(1:nrow(xgb.lat.cv), xgb.lat.cv$test.rmse.mean, xlab = "Ntrees", ylab = "RMSE", main = "CV - lat", type = "l")
      plot(1:nrow(xgb.lon.cv), xgb.lon.cv$test.rmse.mean, xlab = "Ntrees", ylab = "RMSE", main = "CV - lon", type = "l")
      
      print(paste0("Ntrees for lat: ", min_rmse_ntrees_lat))
      print(paste0("Ntrees for lon: ", min_rmse_ntrees_lon))
      
      if(save_grids) {
            tmp_grid_lat <<- rmse_grid_lat
            tmp_grid_lon <<- rmse_grid_lon
            message("Grids saved to 'tmp_grid_[lat/lon]' objects")
      }
      
      par(mfrow = c(1, 1))
}

# ========== PLOT DATA ========== #

# Plot the points
map <- get_map("tanzania", maptype = "terrain", zoom = 6)
ggmap(map) + geom_point(data = wpd[wpd$country_name == "Tanzania, United Republic of", ], aes(x = lon_deg, y = lat_deg), alpha = 0.1, col = "blue") + 
      xlab(expression(`Longitude `(degree))) + ylab(expression(`Latitude `(degree))) + 
      ggtitle("Locations of mapped water points") + theme(axis.title.x = element_text(size = 16), 
                                                          axis.title.y = element_text(size = 16), 
                                                          plot.title = element_text(face = "bold", size = 20), 
                                                          legend.title = element_text(face = "bold"))

# Water access world map
map.world <- map_data(map = "world")
map.world$percentage <- numeric(nrow(map.world))
country_names <- unique(water$Country.Name)

for (i in 1:length(country_names)) {
      map.world$percentage[map.world$region == water$Country.Name[i]] <- water$X2015[i]
}

# Manually add several data values
map.world$percentage[map.world$percentage == 0] <- NA
map.world$percentage[map.world$region == "USA"] <- 99.2
map.world$percentage[map.world$region == "Russia"] <- 96.9
map.world$percentage[map.world$region == "UK"] <- 100
map.world$percentage[map.world$region == "Slovakia"] <- 100
map.world$percentage[map.world$region == "Bulgaria"] <- 99
map.world$percentage[map.world$region == "Kosovo"] <- 80
map.world$percentage[map.world$region == "Serbia"] <- 99
map.world$percentage[map.world$region == "Macedonia"] <- 99
map.world$percentage[map.world$region == "Saudi Arabia"] <- 97
map.world$percentage[map.world$region == "Iran"] <- 99

# Plot the map
gg <- ggplot()
gg <- gg + geom_map(data = map.world, map = map.world, aes(x = long, y = lat, map_id = region), fill = "#ffffff", 
                    color = NA)
gg <- gg + geom_map(data = map.world, map = map.world, color = "white", size = 0.15, aes(fill = percentage, 
                                                                                         map_id = region))
gg <- gg + theme_minimal() + xlab(expression(`Longitude `(degree))) + ylab(expression(`Latitude `(degree))) + 
      ggtitle("Access to an improved water source in 2015") + 
      theme(axis.title = element_blank(), plot.title = element_text(size = 24, face = "bold"), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), axis.text = element_blank(), 
            legend.title = element_text(face = "bold"), legend.position = c("bottom")) + 
      scale_fill_continuous(name = "% of population with access to\nan improved water source")
print(gg)

# ========== DATA QUALITY ASSESSMENT ========== #

# a) contextual - check how many rows are incomplete (at least 1 NA)
na_mat <- sapply(tz, is.na)

# ========== EXPLORATORY ANALYSIS ========== #

# Get locations of 10 biggest cities (by population)
address <- c("Dar es Salaam", "Mwanza", "Zanzibar City", "Arusha", "Mbeya", "Morogoro", "Tanga", "Dodoma", 
             "Kigoma", "Moshi, Kilimanjaro")
locations <- ldply(address, function(x) geoCode(x))
names(locations) <- c("lat", "lon", "location_type", "formatted")

# Create a distance matrix for top 10 cities
dist_mat_hav <- matrix(nrow = nrow(locations), ncol = nrow(tz))
dist_mat_euclid <- matrix(nrow = nrow(locations), ncol = nrow(tz))

for (i in 1:nrow(locations)) {
      dist_mat_hav[i, ] <- dist_haversine(point_lat = tz$lat_deg, point_lon = tz$lon_deg, city_coords = as.numeric(locations[i, c(1:2)]))
      rownames(dist_mat_hav) <- unlist(strsplit(locations$formatted, split = ","))[seq(1, (nrow(locations)) * 2, by = 2)]
}

for (i in 1:nrow(locations)) {
      dist_mat_euclid[i, ] <- dist_euclid(point_lat = tz$lat_deg, point_lon = tz$lon_deg, city_coords = as.numeric(locations[i, c(1:2)]))
      rownames(dist_mat_euclid) <- unlist(strsplit(locations$formatted, split = ","))[seq(1, (nrow(locations)) * 2, by = 2)]
}

dist_mat_euclid <- t(dist_mat_euclid)
dist_mat_hav <- t(dist_mat_hav)

# Difference between hav and euclid dist
errors_distance_df <- as.data.frame(abs(dist_mat_euclid - dist_mat_hav))
errors_distance_df_melt <- melt(errors_distance_df)

ggplot(errors_distance_df_melt, aes(x = variable, y = value)) + geom_boxplot() + theme_minimal() + xlab("City") + ylab("Difference(km)") + 
      ggtitle("Absolute difference between euclidean and haversine distance for top 10 cities") + 
      theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), plot.title = element_text(face = "bold", size = 20)) + 
      scale_y_continuous(breaks = seq(0, max(errors_distance_df_melt$value), by = 1))

# How does install_year influence status_functional?
functional_by_year$color[functional_by_year$diff < 0] <- "red"
functional_by_year$color[functional_by_year$diff >= 0] <- "green"

ggplot(functional_by_year_melt[functional_by_year_melt$variable != "diff",], aes(x = year, y = value)) + 
      geom_bar(aes(fill = variable), stat = "identity", position = "dodge") + xlab("Year") + ylab("Number of point water sources") + 
      ggtitle("Functional and non-functional wells by year of installation") + coord_flip() + 
      theme_minimal() +theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.title = element_text(size = 24, face = "bold")) +
      scale_fill_manual("Status", labels = c("Functional", "Non-functional"), values = c("green3", "red3")) 
#doesn't tell us much, let's see the difference

ggplot(functional_by_year_melt[functional_by_year_melt$variable == "diff",], aes(x = year, y = value)) + 
      geom_bar(aes(fill = functional_by_year$color), stat = "identity", position = "dodge") + xlab("Year") + ylab("Difference") + 
      ggtitle("Difference between functional and non-functional point water sources by year of installation") + 
      theme_minimal() + theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.title = element_text(size = 24, face = "bold")) + 
      scale_fill_manual("Functional - not functional", labels = c("More functional", "More non-functional"), values = c("green3", "red3")) + 
      geom_vline(xintercept = 39, lwd = 1, lty = 2, col = "steelblue3") + coord_flip()
# now that's better - in past years non-functional dominate, while in recent years more functional, so it is seems like the more recently-built wells are indeed more likely to be functional

# Let's now see distribution of each water source type
ggplot(tz, aes(x = lon_deg, y = lat_deg)) + 
      geom_point(aes(color = water_source), alpha = 0.3) + theme_minimal() + 
      xlab("Lat") + ylab("Lon") + ggtitle("Water point data by source type") + 
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold")) + 
      guides(col = guide_legend(title = "Water source type"), face = "bold")


# ========== TRANSFORM DATA ========== #

# Nearest city and distance cols
nearest_city <- numeric(nrow(tz))
nearest_dist <- numeric(nrow(tz))
for(i in 1:nrow(tz)) {
      nearest_city[i] <- names(which.min(dist_mat_hav[i,]))
      nearest_dist[i] <- min(dist_mat_hav[i,])
}

tz$nearest_city <- as.factor(nearest_city)
tz$nearest_dist <- nearest_dist

# Get well status
tz$status_functional <- numeric(nrow(tz))
tz$status_functional[grep("status:functional", tz$status, ignore.case = T)] <- 1
tz$status_functional[grep("status:not functional", tz$status, ignore.case = T)] <- 0

# Get well quality
tz$status_quality <- character(nrow(tz))
status_vec <- unique(unlist(strsplit(tz$status, split = "\\|")))
quality_vec <- status_vec[grep("^quality", status_vec, ignore.case = T)]
for(i in 1:length(quality_vec)) {
      tz$status_quality[grep(quality_vec[i], tz$status, ignore.case = T)] <- quality_vec[i]
}
tz$status_quality[tz$status_quality == ""] <- "Quality not indicated"

tz$status_quality <- as.factor(tz$status_quality)

# Get well quantity
tz$status_quantity <- character(nrow(tz))
quantity_vec <- status_vec[grep("^quantity", status_vec, ignore.case = T)]
for(i in 1:length(quantity_vec)) {
      tz$status_quantity[grep(quantity_vec[i], tz$status, ignore.case = T)] <- quantity_vec[i]
}
tz$status_quantity[tz$status_quantity == ""] <- "Quantity not indicated"
tz$status_quantity <- as.factor(tz$status_quantity)

# Is depth < 50m
idx <- tz$groundwater_depth == "VS" | tz$groundwater_depth == "S" | tz$groundwater_depth == "SM"
tz$is50less <- integer(nrow(tz))
tz$is50less[idx] <- 1

# Normalise - mean 0 and st.dev 1
tz$lat_deg_stand <- standardise(tz$lat_deg)
tz$lon_deg_stand <- standardise(tz$lon_deg)

# See how many clusters for k-means
wssplot(tz[, c("lon_deg_stand", "lat_deg_stand")], 50)

k <- 20
set.seed(1)
clust <- kmeans(tz[, c("lat_deg_stand", "lon_deg_stand")], k, iter.max = 1000, nstart = 100)
ggplot(tz, aes(x = lon_deg, y = lat_deg)) + geom_point(aes(color = as.factor(clust$cluster))) + theme_minimal() +
      xlab("Lat") + ylab("Lon") + ggtitle(paste0("Clusters for k = ", k)) + guides(col = guide_legend(ncol = 2, title = "Cluster ID")) +
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold"))

# Add cluster id + size features
tz$cluster_id <- as.factor(clust$cluster)
tz$cluster_size <- NULL
tz$cluster_center_lat <- NULL
tz$cluster_centre_lon <- NULL
for (i in 1:length(clust$size)) {
      tz$cluster_size[tz$cluster_id == i] <- clust$size[i]
      tz$cluster_center_lat[tz$cluster_id == i] <- clust$centers[i,1]
      tz$cluster_center_lon[tz$cluster_id == i] <- clust$centers[i,2]
}
tz$cluster_size <- as.factor(tz$cluster_size)

tz$is_1992_later <- ifelse(tz$install_year >= 1992, 1, 0)

# Df of type by cluster
type_cluster_df <- melt(as.data.frame(table(tz$water_source, tz$cluster_id)))
colnames(type_cluster_df) <- c("water_source", "cluster_id", "variable", "value")

ggplot(type_cluster_df, aes(x = cluster_id, y = value, fill = reorder(water_source, -value))) + geom_bar(stat = "identity") + theme_minimal() +
      xlab("Cluster ID") + ylab("Number of data points") + ggtitle("Number of water sources by types for each cluster") +
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold")) +
      scale_fill_manual(name = "Water source type", values = brewer.pal(4, "Set2"))

# Some data type conversion
tz$groundwater_depth <- as.character(tz$groundwater_depth)
tz$groundwater_productivity <- as.character(tz$groundwater_productivity)
tz$groundwater_storage <- as.character(tz$groundwater_storage)
tz$nearest_city <- as.character(tz$nearest_city)

tz$groundwater_depth[tz$groundwater_depth == "VS"] <- "1"
tz$groundwater_depth[tz$groundwater_depth == "S"] <- "2"
tz$groundwater_depth[tz$groundwater_depth == "SM"] <- "3"
tz$groundwater_depth[tz$groundwater_depth == "M"] <- "4"
tz$groundwater_depth[tz$groundwater_depth == "D"] <- "5"
tz$groundwater_depth[tz$groundwater_depth == "VD"] <- "6"
tz$groundwater_depth <- as.numeric(tz$groundwater_depth)

tz$groundwater_productivity[tz$groundwater_productivity == "VL"] <- "1"
tz$groundwater_productivity[tz$groundwater_productivity == "L"] <- "2"
tz$groundwater_productivity[tz$groundwater_productivity == "LM"] <- "3"
tz$groundwater_productivity[tz$groundwater_productivity == "M"] <- "4"
tz$groundwater_productivity[tz$groundwater_productivity == "H"] <- "5"
tz$groundwater_productivity[tz$groundwater_productivity == "VH"] <- "6"
tz$groundwater_productivity <- as.numeric(tz$groundwater_productivity)

tz$groundwater_storage[tz$groundwater_storage == "0"] <- "0"
tz$groundwater_storage[tz$groundwater_storage == "L"] <- "1"
tz$groundwater_storage[tz$groundwater_storage == "LM"] <- "2"
tz$groundwater_storage[tz$groundwater_storage == "M"] <- "3"
tz$groundwater_storage[tz$groundwater_storage == "H"] <- "4"
tz$groundwater_storage[tz$groundwater_storage == "VH"] <- "5"
tz$groundwater_storage <- as.numeric(tz$groundwater_storage)

tz$cluster_id <- as.numeric(tz$cluster_id)
tz$cluster_size <- as.numeric(tz$cluster_size)

cities <- unique(tz$nearest_city)

for(i in 1:length(cities)) {
      tz$nearest_city[tz$nearest_city == cities[i]] <- as.character(i)
}

tz$nearest_city <- as.numeric(tz$nearest_city)

# ========== ML ========== #
set.seed(1)
# Partition into training and test sets
train_id <- sample(nrow(tz), 0.8*nrow(tz)) # 80% train + validation sets
tz_test <- tz[-train_id,] # 20% test set
tz <- tz[train_id,] # Train + validation set
set.seed(1)
train <- sample(nrow(tz), 0.8*nrow(tz)) # Train idx (80% of 80% = 64% of all data is train and 16% is validation)

# Summary of splits (% of original data):
# Train: 64% 
# Validation: 16%
#----------------
# Train + validation: 80%
# Test: 20%

# ---------- Random forest ---------- #

mtry = 5
ntree = 2000
nodesize = 5

# Train
rf.lat <- randomForest(lat_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                             nearest_city + nearest_dist + is_1992_later, tz[train,], importance = TRUE, replace = FALSE, 
                       mtry = mtry, ntree = ntree, nodesize = nodesize)
rf.lon <- randomForest(lon_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                             nearest_city + nearest_dist + is_1992_later, tz[train,], importance = TRUE, replace = FALSE, 
                       mtry = 5, ntree = ntree, nodesize = nodesize)
rf.lat.pred <- predict(rf.lat, tz[-train,])
rf.lon.pred <- predict(rf.lon, tz[-train,])

# Metrics
return_metrics(rf.lat.pred, rf.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])

preds_rf_df <- data.frame(actual_lat = tz$lat_deg[-train], actual_lon = tz$lon_deg[-train], predicted_lat = rf.lat.pred, predicted_lon = rf.lon.pred)
errors_rf <- errors_df(rf.lat.pred, rf.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])

# Plot predictions against actual locations
ggplot() + 
      geom_point(data = preds_rf_df, aes(x = actual_lon, y = actual_lat, color = "darkgreen")) + 
      geom_point(data = preds_rf_df, aes(x = predicted_lon, y = predicted_lat, color = "darkred"), alpha = 0.5) +
      xlab("Lon") + ylab("Lat") + ggtitle("Actual and predicted locations for Random Forest") + theme_minimal() +
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold")) + 
      scale_color_manual(name = "Location type", values = c("darkgreen", "darkred"), labels = c("Actual location", "Predicted location"))

# Perhaps add some jitter?
# rf_mse_jitter <- add_jitter(actual_lat = preds_rf_df$actual_lat, actual_lon = preds_rf_df$actual_lon, 
#                             predicted_lat = preds_rf_df$predicted_lat, predicted_lon = preds_rf_df$predicted_lon)
# 
# plot_with_jitter(actual_lat = preds_rf_df$actual_lat, actual_lon = preds_rf_df$actual_lon, 
#                  predicted_lat = preds_rf_df$predicted_lat, predicted_lon = preds_rf_df$predicted_lon, 
#                  jitter_amount = 0.1, model = "random forest", alpha = 0.3) +
#       theme(axis.title.x = element_text(size = 16), 
#             axis.title.y = element_text(size = 16), 
#             plot.title = element_text(face = "bold", size = 20), 
#             legend.title = element_text(face = "bold"))

# Plot error as a function of 'm' (number of variables considered at each random split)?

# N <- seq(1, 7, by = 1)
# errors_big_df <- data.frame(MAE_lat = numeric(length(N)),
#                             MAE_lon = numeric(length(N)),
#                             MSE_lat = numeric(length(N)),
#                             MSE_lon = numeric(length(N)),
#                             RMSE_lat = numeric(length(N)),
#                             RMSE_lon = numeric(length(N)))
# iter = 1
# pb <- txtProgressBar(min = 1, max = length(N), style = 3)
# for(i in N){
#       # Random forest
#       rf.lat <- randomForest(lat_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size +
#                                    nearest_city + nearest_dist, tz[train,], mtry = i)
#       rf.lon <- randomForest(lon_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size +
#                                    nearest_city + nearest_dist, tz[train,], mtry = i)
#       rf.lat.pred <- predict(rf.lat, tz[-train,])
#       rf.lon.pred <- predict(rf.lon, tz[-train,])
#       # MAE, MSE, RMSE
#       loss_rf <- data.frame(MAE_lat = mae(rf.lat.pred, tz$lat_deg[-train]),
#                             MAE_lon = mae(rf.lon.pred, tz$lon_deg[-train]),
#                             MSE_lat = mse(rf.lat.pred, tz$lat_deg[-train]),
#                             MSE_lon = mse(rf.lon.pred, tz$lon_deg[-train]),
#                             RMSE_lat = rmse(rf.lat.pred, tz$lat_deg[-train]),
#                             RMSE_lon = rmse(rf.lon.pred, tz$lon_deg[-train]))
#       errors_big_df[iter,] <- loss_rf
#       iter = iter + 1
#       setTxtProgressBar(pb, i)
# }
# close(pb)

errors_m <- readRDS("errors_rf_m.Rds")
errors_melt <- melt(errors_m)
errors_melt$m <- rep(c(1:7), 3)

# Plot the errors as a function of 'm'
ggplot(errors_melt, aes(x = m, y = value, color = variable)) + 
      geom_line(lwd = 1) + theme_minimal() + xlab("M") + ylab("Value of metric") + ggtitle("Loss function for varying values of M") +
      scale_x_continuous(breaks = seq(1, 7, by = 1)) + 
      scale_y_continuous(breaks = seq(0, 0.3, by = 0.02)) + 
      scale_color_manual(name = "Variable", 
                         labels = c("MAE lat", "MAE lon", "MSE lat", "MSE lon", "RMSE lat", "RMSE lon"), 
                         values = c("gold", "gold3", "steelblue1", "steelblue4", "green3", "darkgreen")) +
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold"))

error_cdf_plot(errors_rf$lat, errors_rf$lon, model = "Random Forest", AUC_max_km = 12)



# ------- GBM ------- #

# quick
n.trees = 1000
int.depth = 10
n.minobsinnode = 5
shrinkage = 0.25

# slightly better but takes ages
n.trees = 10000
int.depth = 30
n.minobsinnode = 20
shrinkage = 0.01

# Train
gbm.lat <- gbm(lat_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                     nearest_city + nearest_dist, data = tz[train,], 
               distribution = "gaussian", n.trees = n.trees, interaction.depth = int.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage)
gbm.lon <- gbm(lon_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                     nearest_city + nearest_dist, data = tz[train,], 
               distribution = "gaussian", n.trees = n.trees, interaction.depth = int.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage)

gbm.lat.pred <- predict(gbm.lat, tz[-train,], n.trees = n.trees)
gbm.lon.pred <- predict(gbm.lon, tz[-train,], n.trees = n.trees)

# Metrics
return_metrics(gbm.lat.pred, gbm.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])
preds_gbm_df <- data.frame(actual_lat = tz$lat_deg[-train], actual_lon = tz$lon_deg[-train], 
                           predicted_lat = gbm.lat.pred, predicted_lon = gbm.lon.pred)
errors_gbm <- errors_df(gbm.lat.pred, gbm.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])

# Plot predictions against actual locations
ggplot() + 
      geom_point(data = preds_gbm_df, aes(x = actual_lon, y = actual_lat, color = "darkgreen")) + 
      geom_point(data = preds_gbm_df, aes(x = predicted_lon, y = predicted_lat, color = "darkred")) +
      xlab("Lon") + ylab("Lat") + ggtitle("Actual and predicted locations for GBM") + theme_minimal() +
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold")) + 
      scale_color_manual(name = "Location type", values = c("darkgreen", "darkred"), labels = c("Actual location", "Predicted location"))

# Perhaps add some jitter?
# gbm_mse_jitter <- add_jitter(actual_lat = preds_gbm_df$actual_lat, actual_lon = preds_gbm_df$actual_lon, 
#                              predicted_lat = preds_gbm_df$predicted_lat, predicted_lon = preds_gbm_df$predicted_lon)
# 
# plot_with_jitter(actual_lat = preds_gbm_df$actual_lat, actual_lon = preds_gbm_df$actual_lon, 
#                  predicted_lat = preds_gbm_df$predicted_lat, predicted_lon = preds_gbm_df$predicted_lon, 
#                  jitter_amount = 0.1, model = "GBM", alpha = 0.3) +
#       theme(axis.title.x = element_text(size = 16), 
#             axis.title.y = element_text(size = 16), 
#             plot.title = element_text(face = "bold", size = 20), 
#             legend.title = element_text(face = "bold"))

error_cdf_plot(errors_gbm$lon, errors_gbm$lat, model = "GBM", AUC_max_km = 12)



# ---------- Lasso-regularised regression ---------- #

# # Train
# x.lat <- model.matrix(lat_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id +
#                             cluster_size + nearest_city + nearest_dist, tz)
# x.lon <- model.matrix(lon_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id +
#                             cluster_size + nearest_city + nearest_dist, tz)
# y.lat <- tz$lat_deg
# y.lon <- tz$lon_deg
# lasso.lat <- glmnet(x.lat[train,], y.lat[train], alpha = 1)
# lasso.lon <- glmnet(x.lon[train,], y.lon[train], alpha = 1)
# cv.lat <- cv.glmnet(x.lat[train,], y.lat[train], alpha = 1)
# cv.lon <- cv.glmnet(x.lon[train,], y.lon[train], alpha = 1)
# bestlam.lat <- cv.lat$lambda.min
# bestlam.lon <- cv.lon$lambda.min
# lasso.lat.pred <- predict(lasso.lat, s = bestlam.lat, newx = x.lat[-train, ])
# lasso.lon.pred <- predict(lasso.lon, s = bestlam.lon, newx = x.lon[-train, ])
# 
# # Metrics
# return_metrics(lasso.lat.pred, lasso.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])
# errors_las <- errors_df(lasso.lat.pred, lasso.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])
# names(errors_las) <- c("lat", "lon")
# preds_las_df <- data.frame(actual_lat = tz$lat_deg[-train], actual_lon = tz$lon_deg[-train],
#                            predicted_lat = lasso.lat.pred, predicted_lon = lasso.lon.pred)
# names(preds_las_df)[3:4] <- c("predicted_lat", "predicted_lon")
# 
# # Plot predictions against actual locations
# ggplot() +
#       geom_point(data = preds_las_df, aes(x = actual_lat, y = actual_lon), color = "darkgreen") +
#       geom_point(data = preds_las_df, aes(x = predicted_lat, y = predicted_lon), color = "darkred") +
#       xlab("Lat") + ylab("Lon") + ggtitle("Actual vs predicted lat/lon for no-jitter lasso predictions") + theme_minimal() +
#       theme(axis.title.x = element_text(size = 16),
#             axis.title.y = element_text(size = 16),
#             plot.title = element_text(face = "bold", size = 20),
#             legend.title = element_text(face = "bold"))
# 
# # Perhaps add some jitter?
# las_mse_jitter <- add_jitter(actual_lat = preds_las_df$actual_lat, actual_lon = preds_las_df$actual_lon,
#                              predicted_lat = preds_las_df$predicted_lat, predicted_lon = preds_las_df$predicted_lon)
# 
# plot_with_jitter(actual_lat = preds_las_df$actual_lat, actual_lon = preds_las_df$actual_lon,
#                  predicted_lat = preds_las_df$predicted_lat, predicted_lon = preds_las_df$predicted_lon,
#                  jitter_amount = 0.1, model = "lasso", alpha = 0.2) +
#       theme(axis.title.x = element_text(size = 16),
#             axis.title.y = element_text(size = 16),
#             plot.title = element_text(face = "bold", size = 20),
#             legend.title = element_text(face = "bold"))
# 
# error_cdf_plot(errors_las$lat, errors_las$lon, model = "Lasso")



# ------- XGBOOST ------- #

model.mat <- as.matrix(tz[,c("groundwater_depth", "groundwater_productivity", "groundwater_storage", 
                             "cluster_id", "cluster_size", "nearest_city", "nearest_dist", 
                             "is50less", "is_1992_later")])

params = list(eta = 0.05,
              max_depth = 17,
              colsample_bytree = 1,
              min_child_weight = 4,
              subsample = 1,
              colsample_bylevel = 0.5, 
              lambda = 0.7)

target.lat <- as.matrix(tz[,"lat_deg"])
target.lon <- as.matrix(tz[,"lon_deg"])

# # CV - nrounds =~ 
# xgb.lat.cv <- xgb.cv(data = model.mat[train,], label = target.lat[train],
#                    nrounds = 1000, params = params, metrics = "rmse", nfold = 10)
# 
# xgb.lon.cv <- xgb.cv(data = model.mat[train,], label = target.lon[train],
#                    nrounds = 1000, params = params, metrics = "rmse", nfold = 10)

# Train
xgb.lat <- xgboost(data = model.mat[train,], label = target.lat[train], 
                   nrounds = 226, params = params, metrics = "rmse", verbose = 0, objective = "reg:linear")

xgb.lon <- xgboost(data = model.mat[train,], label = target.lon[train], 
                   nrounds = 356, params = params, metrics = "rmse", verbose = 0, objective = "reg:linear")

xgb.lat.pred <- predict(xgb.lat, model.mat[-train,])
xgb.lon.pred <- predict(xgb.lon, model.mat[-train,])

# Metrics
return_metrics(xgb.lat.pred, xgb.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])
preds_xgb_df <- data.frame(actual_lat = tz$lat_deg[-train], actual_lon = tz$lon_deg[-train], predicted_lat = xgb.lat.pred, predicted_lon = xgb.lon.pred)
errors_xgb <- errors_df(xgb.lat.pred, xgb.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])

# Variable importance
xgb_importance_lat <- xgb.importance(colnames(model.mat), "xgb.lat.dump", xgb.lat)
xgb_importance_lon <- xgb.importance(colnames(model.mat), "xgb.lon.dump", xgb.lon)
xgb.plot.importance(xgb_importance_lat)
xgb.plot.importance(xgb_importance_lon)

# Plot predictions against actual locations
ggplot() + 
      geom_point(data = preds_xgb_df, aes(x = actual_lon, y = actual_lat, color = "darkgreen")) + 
      geom_point(data = preds_xgb_df, aes(x = predicted_lon, y = predicted_lat, color = "darkred")) +
      xlab("Lon") + ylab("Lat") + ggtitle("Actual and predicted locations for XGB") + theme_minimal() +
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold")) + 
      scale_color_manual(name = "Location type", values = c("darkgreen", "darkred"), labels = c("Actual location", "Predicted location"))

error_cdf_plot(errors_xgb$lon, errors_xgb$lat, model = "xgboost", AUC_max_km = 12)

# # ANN
# x.lat <- model.matrix(~ lat_deg + groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id +
#                             cluster_size + nearest_city + nearest_dist, tz)
# x.lon <- model.matrix(~ lon_deg + groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id +
#                             cluster_size + nearest_city + nearest_dist, tz)
# 
# colnames(x.lat)[54] <- "nearest_cityDaresSalaam"
# n <- colnames(x.lat)[-c(1:2)]
# f <- as.formula(paste("lat_deg ~", paste(n[!n %in% "lat_deg"], collapse = " + ")))
# nn <- neuralnet(f, x.lat[1:1000,], hidden = c(2,2), learningrate = NULL)
# pr.nn <- compute(nn, x.lat[1000:2000,-c(1:2)])
# mse(pr.nn$net.result, x.lat[1000:2000, 2])


# WHAT IF WE TAKE A WEIGHTED AVERAGE OF PREDICTIONS? (hint: it's even better!)

weight_rf <- seq(0, 1, by = 0.01)
weight_xgb <- seq(0, 1, by = 0.01)
len_rf <- length(weight_rf)
len_xgb <- length(weight_xgb)
lat_grid <- matrix(nrow = len_rf, ncol = len_xgb)
lon_grid <- lat_grid

pb <- txtProgressBar(min = 0, max = len_rf * len_xgb, style = 3)

color_ramp <- blue2red(len_rf * len_xgb)

iter = 1

for(i in 1:len_rf) {
      for(j in 1:len_xgb) {
            mean.lat.pred <- (weight_rf[i] * rf.lat.pred + weight_xgb[j] * xgb.lat.pred)/(weight_rf[i] + weight_xgb[j])
            mean.lon.pred <- (weight_rf[i] * rf.lon.pred + weight_xgb[j] * xgb.lon.pred)/(weight_rf[i] + weight_xgb[j])
            lat_grid[i, j] <- return_metrics(mean.lat.pred, mean.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])$RMSE_lat
            lon_grid[i, j] <- return_metrics(mean.lat.pred, mean.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])$RMSE_lon
            iter = iter + 1
            setTxtProgressBar(pb, iter)
      }
}

image(lat_grid, xlab = "RF weight", ylab = "XGBoost weight", main = "RMSE for Latitude", col = color_ramp, axes = F)
axis(1, at = seq(0, 1, length.out = len_rf), labels = weight_rf)
axis(2, at = seq(0, 1, length.out = len_xgb), labels = weight_xgb)
image(lon_grid, xlab = "RF weight", ylab = "XGBoost weight", main = "RMSE for Longitude", col = color_ramp, axes = F)
axis(1, at = seq(0, 1, length.out = len_rf), labels = weight_rf)
axis(2, at = seq(0, 1, length.out = len_xgb), labels = weight_xgb)

lat_best <- arrayInd(which.min(lat_grid), dim(lat_grid))
lon_best <- arrayInd(which.min(lon_grid), dim(lon_grid))

lat_rf_best <- weight_rf[lat_best[1]]
lat_xgboost_best <- weight_rf[lat_best[2]]

lon_rf_best <- weight_rf[lon_best[1]]
lon_xgboost_best <- weight_rf[lon_best[2]]

best_df <- data.frame(rf = c(lat_rf_best, lon_rf_best), xgboost = c(lat_xgboost_best, lon_xgboost_best))
rownames(best_df) <- c("lat", "lon")
best_df

# Metrics
avg.lat.pred <- (best_df[1, 1] * rf.lat.pred + best_df[1, 2] * xgb.lat.pred)/sum(best_df[1,])
avg.lon.pred <- (best_df[2, 1] * rf.lon.pred + best_df[2, 2] * xgb.lon.pred)/sum(best_df[2,])
return_metrics(pred.lat = avg.lat.pred, pred.lon = avg.lon.pred, actual.lat = tz$lat_deg[-train], actual.lon = tz$lon_deg[-train])
preds_avg_df <- data.frame(actual_lat = tz$lat_deg[-train], actual_lon = tz$lon_deg[-train], predicted_lat = avg.lat.pred, predicted_lon = avg.lon.pred)
errors_avg <- errors_df(predicted_lat = avg.lat.pred, predicted_lon = avg.lon.pred, actual_lat = tz$lat_deg[-train], actual_lon = tz$lon_deg[-train])

# Plot predictions against actual locations
ggplot() + 
      geom_point(data = preds_avg_df, aes(x = actual_lon, y = actual_lat, color = "darkgreen")) + 
      geom_point(data = preds_avg_df, aes(x = predicted_lon, y = predicted_lat, color = "darkred")) +
      xlab("Lon") + ylab("Lat") + ggtitle("Actual and predicted locations for weighted average\nof Random Forest and XGB models") + theme_minimal() +
      theme(axis.title.x = element_text(size = 16), 
            axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), 
            legend.title = element_text(face = "bold")) + 
      scale_color_manual(name = "Location type", values = c("darkgreen", "darkred"), labels = c("Actual location", "Predicted location"))

error_cdf_plot(errors_avg$lat, errors_avg$lat, model = "Weighted average", AUC_max_km = 12)

# Compare all model's distribution of errors
compare_error_dist(errors_rf, errors_xgb, errors_gbm, errors_avg, 
                   model1 = "Random Forest", model2 = "XGB", model3 = "GBM", model4 = "Averaged RF & XGB")

# ========== RESULTS ANALYSIS ========== #

# Make a data frame of all metrics for all models
all_metrics_df <- data.frame(MAE_lat = numeric(4), MAE_lon = numeric (4), MSE_lat = numeric(4), 
                             MSE_lon = numeric(4), RMSE_lat = numeric(4), RMSE_lon = numeric(4))
all_metrics_df[1,] <- unname(return_metrics(rf.lat.pred, rf.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train]))
all_metrics_df[2,] <- unname(return_metrics(gbm.lat.pred, gbm.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train]))
all_metrics_df[3,] <- unname(return_metrics(xgb.lat.pred, xgb.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train]))
all_metrics_df[4,] <- unname(return_metrics(avg.lat.pred, avg.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train]))
rownames(all_metrics_df) <- c("random_forest", "GBM", "XGB", "weighted_average")

all_metrics_df_melt <- melt(all_metrics_df)
all_metrics_df_melt$model <- c("random_forest", "GBM", "XGB", "weighted_average")
all_metrics_df_melt$variable <- gsub("_", " ", all_metrics_df_melt$variable)

ggplot(all_metrics_df_melt, aes(x = variable, y = value, fill = model)) + geom_bar(stat = "identity", position = "dodge") + 
      xlab("Metric") + ylab("Value") + ggtitle("Summary of performance of all models") + theme_minimal() + 
      theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
            plot.title = element_text(face = "bold", size = 20), legend.title = element_text(face = "bold")) +
      scale_fill_manual(name = "Model", values = c("steelblue1", "darkorange", "gray", "darkgreen"), labels = c("GBM", "RF", "Weighted average", "XGB"))