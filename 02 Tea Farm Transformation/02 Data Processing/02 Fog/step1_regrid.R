# Install necessary packages if you haven't already
install.packages("ncdf4")
install.packages("terra")
install.packages("akima") # Optional, for specific interpolation methods

library(ncdf4)
library(terra)
library(akima) # For interp() function, if needed

# --- 1. Read the NetCDF data ---
nc_file <- "C:/Users/user/Downloads/wc_v0_8_fog_ratio_year_2018.nc" # Replace with your NetCDF file path
# nc_file <- "C:/Users/user/Downloads/pgw_v0_8_fog_ratio_year_2018.nc" # Replace with your NetCDF file path
nc <- nc_open(nc_file)

# Get variable information
var_name <- "yearly_fog_ratio" # Replace with your data variable name
data_var <- ncvar_get(nc, var_name)

# Get irregular lat and lon
# Assuming lat and lon are 2D arrays corresponding to your data
lat_var <- ncvar_get(nc, "lat") # Replace with your latitude variable name
lon_var <- ncvar_get(nc, "lon") # Replace with your longitude variable name

# Close the NetCDF file
nc_close(nc)

# --- 2. Prepare the irregular data for interpolation ---
# Flatten the 2D lat/lon and data arrays into vectors of points
# Ensure the dimensions match. If your data_var is (lon_dim, lat_dim) and
# lat_var, lon_var are also (lon_dim, lat_dim), they should align naturally
# when flattened. You might need to adjust 'as.vector' order depending on
# how your NetCDF stores the data (e.g., column-major vs row-major).
# Often, NetCDF data is stored in (lon, lat) order.
# If your data is structured differently, you might need to transpose
# data_var, lat_var, and lon_var using t()
# For example:
# lat_var <- t(lat_var)
# lon_var <- t(lon_var)
# data_var <- t(data_var)

points_df <- data.frame(
  lon = as.vector(lon_var),
  lat = as.vector(lat_var),
  value = as.vector(data_var)
)

# Remove any NA values if present, as interpolation functions might struggle
points_df <- na.omit(points_df)

# --- 3. Define the target regular grid ---
# Determine the extent of your data to define the new grid
lon_min <- min(points_df$lon)
lon_max <- max(points_df$lon)
lat_min <- min(points_df$lat)
lat_max <- max(points_df$lat)

# Define the resolution for your regular grid (e.g., 0.25 degrees)
resolution_lon <- 0.00100
resolution_lat <- 0.00100
# resolution_lon <- 0.018018 # Approximately 2km in longitude
# resolution_lat <- 0.018018 # Approximately 2km in latitude

# Create a new, empty SpatRaster for the target grid
# Use the extent from your irregular data, or define a custom extent
target_rast <- rast(
  xmin = lon_min, xmax = lon_max,
  ymin = lat_min, ymax = lat_max,
  res = c(resolution_lon, resolution_lat),
  crs = "EPSG:4326" # WGS84 geographic coordinate system
)

# --- 4. Interpolate the data onto the regular grid ---

# Method 1: Using terra::interpIDW (Inverse Distance Weighted) or terra::interpNear (Nearest Neighbor)
# These are suitable for point data
# Be aware that for interpolation, you might need to handle large datasets
# more carefully, potentially by chunking or using parallel processing.
# For simplicity, we'll do it directly here.

# Convert points_df to a SpatVector object
points_spatvec <- vect(points_df, geom = c("lon", "lat"), crs = "EPSG:4326")

# Perform interpolation
# interpIDW is good for general spatial data. You might need to tune 'p' and 'max_dist'
# If you have time series data, you'll need to loop through each time step
# and interpolate separately.
# interpolated_rast <- interpIDW(points_spatvec, field = "value", rast = target_rast, p = 2) # p=2 is common for IDW

# Method 2: Using akima::interp for more sophisticated interpolation (e.g., bilinear, bicubic)
# This method works well if your irregular data can be conceptualized as scattered points
# that you want to interpolate onto a grid.
# Note: akima::interp requires ordered x and y, and unique values.
# If your irregular grid is truly a mesh, this might be more complex.
# For truly scattered points (e.g., from weather stations), this is very useful.

# If your data is truly an irregular grid (like a curvilinear grid) and you want
# to transform it, you might need a different approach, potentially involving
# defining polygons for each grid cell and then rasterizing/averaging.
# However, if it's more like scattered observations that you want to grid:

 interp_result <- interp(x = points_df$lon,
                         y = points_df$lat,
                         z = points_df$value,
                         xo = unique(xFromCol(target_rast)), # target longitudes
                         yo = unique(yFromRow(target_rast)), # target latitudes
                         linear=TRUE) # or "spline"

# # Convert the akima output to a SpatRaster
 interpolated_rast_akima <- rast(t(interp_result$z[,]),
                                 extent = ext(lon_min, lon_max, lat_min, lat_max),
                                 crs = "EPSG:4326")

#fog_2018_pgw <- interpolated_rast_akima
#fog_2018_wc <- interpolated_rast_akima


#fog_2018_dif <- fog_2018_pgw - fog_2018_wc
# --- 5. (Optional) Write the new NetCDF file ---
output_file <- "D:/TEST/fog_2018_100m_wc.nc"
# output_file <- "D:/TEST/fog_2018_100m_pgw.nc"
writeCDF(interpolated_rast_akima, output_file, varname = var_name, overwrite = TRUE)

message(paste("Successfully re-gridded data to:", output_file))

# You can also plot the result to visualize
plot(interpolated_rast_akima, main = paste("Interpolated Yearly Fog Ratio for 2018 WC"))




