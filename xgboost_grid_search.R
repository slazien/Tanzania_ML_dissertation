grid_search <- function(param1_seq, param2_seq, 
                        param1_name, param2_name, 
                        params_list = list(eta = 0.05, max_depth = 17, 
                                           colsample_bytree = 0.85, min_child_weight = 1, 
                                           subsample = 0.8, colsample_bylevel = 1, alpha = 0.2), 
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
                  axis(1, at = seq(0, max(param1_seq), length.out = len1), labels = param1_seq)
                  axis(2, at = seq(0, max(param2_seq), length.out = len2), labels = param2_seq)
                  image(rmse_grid_lon, main = "Lon", xlab = param1_name, ylab = param2_name, axes = F, col = color_ramp)
                  axis(1, at = seq(0, max(param1_seq), length.out = len1), labels = param1_seq)
                  axis(2, at = seq(0, max(param2_seq), length.out = len2), labels = param2_seq)
                  
                  setTxtProgressBar(pb, iter)
                  iter = iter + 1
            }
      }
      
      grid_sum <- rmse_grid_lat + rmse_grid_lon
      
      min_idx <- arrayInd(which.min(grid_sum), dim(grid_sum))
      
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
      
      gc()
}