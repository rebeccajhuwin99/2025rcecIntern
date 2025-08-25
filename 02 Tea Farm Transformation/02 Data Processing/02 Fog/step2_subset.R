# `terra` is the best tool for handling spatial data (like rasters)
library(terra)

# --- 1. Define Input/Output File Paths ---
# Input WC NetCDF file
wc_100m_file <- "D:/TEST/fog_2018_100m_wc.nc" 

# Define a reference map to set the final output's CRS, extent, and resolution
# This file should be consistent with your other simulation inputs (e.g., land use map)
landuse_map_path <- "C:/Users/user/Downloads/1718_229_198_tree20.tif"

# Define the output path for the final TIFF file
# The filename will be "fog_2018_100m_wc_subset.tif"
output_tif_file <- "D:/TEST/fog_2018_100m_wc_subset.tif"

# Define the output path for the PNG plot file
output_png_file <- "D:/TEST/fog_2018_100m_wc_subset.png"


# --- 2. Read the reference map and get its properties ---
message(paste("正在讀取參考地圖:", landuse_map_path))
if (!file.exists(landuse_map_path)) {
  stop(paste("錯誤: 參考地圖檔案未找到於:", landuse_map_path))
}
lu.map <- rast(landuse_map_path)

# Extract the CRS, extent, and resolution we need from the reference map
target_crs_obj <- crs(lu.map)
target_extent <- ext(lu.map)
target_res <- res(lu.map)

message(paste("目標 CRS:", crs(lu.map, proj=TRUE)))
message(paste("目標範圍:", paste(as.vector(target_extent), collapse=", ")))
message(paste("目標解析度:", target_res[1], "m"))

# --- 3. Read the NetCDF file into a SpatRaster object ---
message(paste("正在讀取 WC 檔案:", wc_100m_file))
if (!file.exists(wc_100m_file)) {
  stop(paste("錯誤: WC NetCDF 檔案未找到於:", wc_100m_file))
}
fog_100m_wc <- rast(wc_100m_file)

# --- 4. Project, Resample, and Crop ---
# We will use the `project()` function to perform all steps at once
# This will transform fog_100m_wc to match the CRS, extent, and resolution of lu.map
message("正在進行投影、重取樣與裁剪...")
fog_subset_resampled <- project(fog_100m_wc, lu.map, method='bilinear')

# --- 5. Write the SpatRaster object to a TIFF file ---
message(paste("正在將資料寫入 TIFF 檔案:", output_tif_file))
writeRaster(fog_subset_resampled, filename = output_tif_file, overwrite = TRUE)

message(paste("成功生成 WC 霧氣 TIFF 檔案:", output_tif_file))


# --- 6. (Optional) Plot the result for visual verification and save as PNG ---
# This section directly writes the plot to a PNG file to avoid the "figure margins too large" error
message(paste("正在將圖表輸出為 PNG 檔案:", output_png_file))

# Open a PNG graphics device
png(filename = output_png_file, width = 1200, height = 900, res = 150)

# Create the plot with an English title
plot(fog_subset_resampled, main="WC Fog Frequency in 2018")

# Close the graphics device to save the file
dev.off()

message(paste("成功生成 WC 霧氣 PNG 檔案:", output_png_file))
