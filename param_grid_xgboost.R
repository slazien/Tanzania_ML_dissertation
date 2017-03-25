seq_subsample <- seq(0.2, 1, by = 0.1)
seq_colsample_bytree <- seq(0.2, 1, by = 0.1)
len <- length(seq_subsample)

param_mat <- t(matrix(c(seq_subsample, seq_colsample_bytree), ncol = 2))
rownames(param_mat) <- c("subsample", "colsample_bytree")
param_mat

loss_lat <- matrix(nrow = len, ncol = len)
loss_lon <- matrix(nrow = len, ncol = len)

target.lat <- as.matrix(tz[,"lat_deg"])
target.lon <- as.matrix(tz[,"lon_deg"])

pb <- txtProgressBar(min = 0, max = len^2, style = 3)

iter = 1

for(i in 1:len) {
      for(j in 1:len) {
            params = list(eta = 0.05,
                          max_depth = 15,
                          colsample_bytree = param_mat[2, j],
                          min_child_weight = 10,
                          subsample = param_mat[1, i],
                          colsample_bylevel = 1,
                          lambda = 2)
            
            xgb.lat <- xgboost(data = model.mat[train,], label = target.lat[train], 
                               nrounds = 500, params = params, metrics = "rmse", verbose = 0)
            
            xgb.lon <- xgboost(data = model.mat[train,], label = target.lon[train], 
                               nrounds = 500, params = params, metrics = "rmse", verbose = 0)
            
            xgb.lat.pred <- predict(xgb.lat, model.mat[-train,])
            xgb.lon.pred <- predict(xgb.lon, model.mat[-train,])
            loss_lat[i, j] <- unname(return_metrics(xgb.lat.pred, xgb.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])$RMSE_lat)
            loss_lon[i, j] <- unname(return_metrics(xgb.lat.pred, xgb.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])$RMSE_lon)
            iter = iter + 1
            setTxtProgressBar(pb, iter)
      }
}