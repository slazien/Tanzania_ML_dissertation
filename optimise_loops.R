find_closest_feature <- function(feature_lat, feature_lon, feature_cat) {
      
}



for(i in 1:nrow(groundwater_depth)) {
      closest_feature[i] <- dist_haversine(tz$lat_deg, tz$lon_deg, as.vector(groundwater_depth[i,c(1:2)]))
}


closest_feature2 <- character(nrow(tz))

cat_vec <- groundwater_storage$cat
groundwater_lat <- groundwater_storage$lat
groundwater_lon <- groundwater_storage$lon
groundwater_feature <- groundwater_storage$cat
tz_latlon <- tz[, .SD, .SDcols = c("lat_deg", "lon_deg")]
tz_latlon <- data.frame(tz_latlon)

pb <- txtProgressBar(0, nrow(tz), style = 3)

for(i in 1:nrow(tz)) {
      closest_feature2[i] <- as.character(groundwater_feature[which.min(dist_haversine(groundwater_lat, groundwater_lon, as.numeric(tz_latlon[i,])))])
      setTxtProgressBar(pb, i)
}