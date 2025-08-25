# regrid_rain_year.R
# This script processes all monthly rain data for 2018, calculates the yearly total,
# regrids it to a 2km resolution, and saves it as an intermediate NetCDF file.

# Load necessary packages
library(ncdf4)
library(terra)
library(akima)

# --- 1. Define Input/Output Files and Paths ---
# Directory containing all monthly NetCDF files for 2018
input_nc_dir <- "D:/20250818_中研院環變周報附件/00 附件/lab resource/2018 2km 氣象資料 ( 年降雨 )/" # <--- 請修改為您的資料夾路徑
output_file_regrid <- "D:/TEST_RAIN/rain_2018_yearly_2km.nc" # <--- 中繼 NetCDF 檔案

# --- 2. Process Monthly Rain Data and Calculate Yearly Total ---
cat("--- Starting Yearly Rain Data Regridding ---\n")
nc_files <- list.files(path = input_nc_dir, pattern = "\\.nc$", full.names = TRUE)

if (length(nc_files) == 0) {
  stop(paste("Error: No NetCDF files found in the directory:", input_nc_dir))
}

# Initialize a raster to accumulate yearly total rain
yearly_rain_total <- NULL
yearly_rain_extent <- NULL
yearly_rain_crs <- NULL

for (file in nc_files) {
  cat(paste("Processing file:", basename(file), "\n"))
  
  # Open the NetCDF file
  nc <- nc_open(file)
  
  # Extract variables
  rain_var_name <- "RAIN"
  lat_var_name <- "LAT"
  lon_var_name <- "LON"
  
  data_var_raw <- ncvar_get(nc, rain_var_name)
  lat_var <- ncvar_get(nc, lat_var_name)
  lon_var <- ncvar_get(nc, lon_var_name)
  
  # Close the NetCDF file
  nc_close(nc)
  
  # Calculate Monthly Total Rainfall
  if (length(dim(data_var_raw)) == 3) {
    data_var_2d <- apply(data_var_raw, c(1, 2), sum, na.rm = TRUE)
  } else if (length(dim(data_var_raw)) == 2) {
    data_var_2d <- data_var_raw
  } else {
    warning(paste("Skipping file", basename(file), "due to unexpected dimensions."))
    next # Skip to next file
  }
  
  # Prepare Data for Interpolation (same logic as before)
  lon_grid <- rep(lon_var, each = length(lat_var))
  lat_grid <- rep(lat_var, times = length(lon_var))
  
  points_df <- data.frame(
    lon = lon_grid,
    lat = lat_grid,
    value = as.vector(data_var_2d)
  )
  points_df <- na.omit(points_df)
  
  # --- 3. Define the Target Grid for Interpolation (2km) ---
  # We define the extent from the first file and reuse it for all subsequent files
  if (is.null(yearly_rain_extent)) {
    yearly_rain_extent <- ext(min(points_df$lon), max(points_df$lon),
                              min(points_df$lat), max(points_df$lat))
    yearly_rain_crs <- "EPSG:4326"
    
    resolution_lon <- 0.0200
    resolution_lat <- 0.0200
    
    target_rast_2km <- rast(
      extent = yearly_rain_extent,
      res = c(resolution_lon, resolution_lat),
      crs = yearly_rain_crs
    )
  }
  
  # --- 4. Interpolate the monthly data onto the 2km grid ---
  interp_result <- interp(x = points_df$lon, y = points_df$lat, z = points_df$value,
                          xo = seq(ext(target_rast_2km)[1], ext(target_rast_2km)[2], by = resolution_lon),
                          yo = seq(ext(target_rast_2km)[3], ext(target_rast_2km)[4], by = resolution_lat),
                          linear = TRUE)
  
  interpolated_rast <- rast(t(interp_result$z),
                            extent = yearly_rain_extent,
                            crs = yearly_rain_crs)
  
  # Accumulate the yearly total
  if (is.null(yearly_rain_total)) {
    yearly_rain_total <- interpolated_rast
  } else {
    yearly_rain_total <- yearly_rain_total + interpolated_rast
  }
}
cat("\nAll monthly files processed. Calculating yearly total...\n")
names(yearly_rain_total) <- "yearly_rain_mm"

# --- 5. Write the final yearly NetCDF file ---
cat(paste("Exporting regridded yearly total rain data to:", output_file_regrid, "\n"))
writeCDF(yearly_rain_total, output_file_regrid, varname = "yearly_rain_mm", overwrite = TRUE)

cat("\n--- Regridding yearly rain step complete! ---\n")
plot(yearly_rain_total, main = "Regridded Yearly Rain (2018) on 2km Grid")