# subset_rain_year.R
# This script takes the regridded yearly NetCDF file, projects, resamples, and subsets it
# to match a target landuse map, then exports the final GeoTIFF.

# Load necessary packages
library(terra)

# --- 1. Define Input/Output Files ---
# Input regridded yearly rain data (from previous step)
regridded_rain_file <- "D:/20250818_中研院環變周報附件/02 找 2018t2024 解/V 2018 降雨資料/rain_2018_yearly_2km.nc" # <--- 請確認路徑

# Landuse map to define extent and target resolution
landuse_map_path <- "C:/Users/user/Downloads/1718_229_198_tree20.tif" # <--- 請確認路徑

# Output path for the final GeoTIFF
output_tif_path <- "D:/20250818_中研院環變周報附件/02 找 2018t2024 解/V 2018 降雨資料/rain_2018_yearly_subset.tif" # <--- 最終輸出的路徑

# --- 2. Read Landuse Map and get its properties ---
cat("--- Starting Yearly Rain Data Subsetting ---\n")
if (!file.exists(landuse_map_path)) {
  stop(paste("Error: Landuse map file not found at:", landuse_map_path))
}
lu.map <- rast(landuse_map_path)

target_crs_obj <- crs(lu.map)
target_extent <- ext(lu.map)
target_crs_string <- crs(lu.map, proj=TRUE)
target_res <- res(lu.map)

cat(paste("Target CRS:", target_crs_string, "\n"))
cat(paste("Target Extent:", paste(as.vector(target_extent), collapse=", "), "\n"))
cat(paste("Target Resolution:", target_res[1], "m\n"))


# --- 3. Import regridded yearly rain data ---
if (!file.exists(regridded_rain_file)) {
  stop(paste("Error: Regridded rain file not found at:", regridded_rain_file))
}
rain_rast <- rast(regridded_rain_file)
cat(paste("Loaded yearly rain data from:", regridded_rain_file, "\n"))

# --- 4. Project, Resample, and Subset the projected rain data ---
cat("Projecting, resampling, and subsetting rain data...\n")
# Create an empty raster with the exact properties of the landuse map
target_grid_6m <- rast(
  extent = target_extent,
  res = target_res,
  crs = target_crs_string
)

# Use project() with the target grid to perform all steps in one go
subset_rain <- project(rain_rast, target_grid_6m, method='bilinear')
cat("Processing complete.\n")


# --- 5. Export as GeoTIFF ---
cat(paste("Exporting final GeoTIFF to:", output_tif_path, "\n"))
writeRaster(subset_rain, filename = output_tif_path, overwrite = TRUE)

cat("\n--- Subsetting step complete! ---\n")

# Plot the final result
plot(subset_rain, main = "Yearly Rain in 2018")