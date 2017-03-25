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