<<<<<<< HEAD
tz_mat <- model.matrix(~ -1 + lat_deg + lon_deg + groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                                        nearest_city + nearest_dist + is_1992_later + is50less, tz)

maxs <- apply(tz_mat, 2, max)
mins <- apply(tz_mat, 2, min)

tz_mat_scaled <- scale(tz_mat, center = mins, scale = maxs - mins)

# train_ <- tz_mat_scaled[train, ]
# test_ <- tz_mat_scaled[-train, ]

nn.lat <- neuralnet(lat_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                          nearest_city + nearest_dist + is_1992_later + is50less, data = tz_mat_scaled[train,], 
                    hidden = c(10, 10), linear.output = T, stepmax = 1e6, lifesign = "full", lifesign.step = 500, threshold = 0.1)


nn.lon <- neuralnet(lon_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                          nearest_city + nearest_dist + is_1992_later + is50less, data = tz_mat_scaled[train,], 
                    hidden = c(10, 10), linear.output = T, stepmax = 1e6, lifesign = "full", lifesign.step = 500, threshold = 0.1)

nn.lat.pred <- compute(nn.lat, tz_mat_scaled[-train, -c(1, 2)])$net.result
nn.lat.pred <- nn.lat.pred*(max(tz$lat_deg)-min(tz$lat_deg))+min(tz$lat_deg)

nn.lon.pred <- compute(nn.lon, tz_mat_scaled[-train, -c(1, 2)])$net.result
nn.lon.pred <- nn.lon.pred*(max(tz$lon_deg)-min(tz$lon_deg))+min(tz$lon_deg)

=======
tz_mat <- model.matrix(~ -1 + lat_deg + lon_deg + groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                                        nearest_city + nearest_dist + is_1992_later + is50less, tz)

maxs <- apply(tz_mat, 2, max)
mins <- apply(tz_mat, 2, min)

tz_mat_scaled <- scale(tz_mat, center = mins, scale = maxs - mins)

# train_ <- tz_mat_scaled[train, ]
# test_ <- tz_mat_scaled[-train, ]

nn.lat <- neuralnet(lat_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                          nearest_city + nearest_dist + is_1992_later + is50less, data = tz_mat_scaled[train,], 
                    hidden = c(10, 10), linear.output = T, stepmax = 1e6, lifesign = "full", lifesign.step = 500, threshold = 0.1)


nn.lon <- neuralnet(lon_deg ~ groundwater_depth + groundwater_productivity + groundwater_storage + cluster_id + cluster_size + 
                          nearest_city + nearest_dist + is_1992_later + is50less, data = tz_mat_scaled[train,], 
                    hidden = c(10, 10), linear.output = T, stepmax = 1e6, lifesign = "full", lifesign.step = 500, threshold = 0.1)

nn.lat.pred <- compute(nn.lat, tz_mat_scaled[-train, -c(1, 2)])$net.result
nn.lat.pred <- nn.lat.pred*(max(tz$lat_deg)-min(tz$lat_deg))+min(tz$lat_deg)

nn.lon.pred <- compute(nn.lon, tz_mat_scaled[-train, -c(1, 2)])$net.result
nn.lon.pred <- nn.lon.pred*(max(tz$lon_deg)-min(tz$lon_deg))+min(tz$lon_deg)

>>>>>>> 400d31d6835eca3597657e8c2efa993c2510df1f
return_metrics(nn.lat.pred, nn.lon.pred, tz$lat_deg[-train], tz$lon_deg[-train])