# 載入必要的套件
library("terra")

# --- 設定您的絕對路徑基本資料夾 ---
# 請將這裡的路徑替換為您電腦中存放所有資料的實際路徑
base_dir <- "D:/TEST"

# 組合檔案的絕對路徑
lu_map_path <- "C:/Users/user/Downloads/2024_229_198_tree20.tif"
fog_dif_path <- "D:/TEST/fog_2018_100m_pgw_wc_dif.nc"

# --- 載入資料 ---
# 使用 terra::rast() 載入資料，這是處理 raster 資料更現代的方法
lu.map <- rast(lu_map_path)
fog <- rast(fog_dif_path)

cat("lu.map 的原始坐標系統:", crs(lu.map), "\n")
cat("fog 的原始坐標系統:", crs(fog), "\n")

# --- 修正步驟：直接進行重新投影，確保坐標系統一致 ---
# 不再做條件判斷，而是直接將 fog 圖層重新投影到 lu.map 的坐標系統
cat("\n--> 正在將 fog 圖層重新投影到 lu.map 的坐標系統...\n")
fog <- project(fog, lu.map)

# 使用 crop() 函數，以 lu.map 的範圍來裁切 fog 圖層
# 這個步驟會先將 fog 圖層的範圍限制在 lu.map 的邊界內
cat("--> 正在以 lu.map 的範圍來裁切 fog 圖層...\n")
sub.fog <- crop(fog, ext(lu.map))

# 使用 resample() 函數，將 sub.fog 圖層的解析度和網格與 lu.map 對齊
# 'bilinear' 方法是一種內插方法，用於處理數值型資料
cat("--> 正在將 sub.fog 圖層的解析度與 lu.map 對齊...\n")
sub.fog <- resample(sub.fog, lu.map, method = 'bilinear')

# --- 輸出 GeoTIFF 檔案 ---
# 建立輸出的檔案名稱
# 修正 fname 未定義的問題，從 fog_dif_path 提取一個有意義的名稱
fname <- basename(fog_dif_path)
output_filename <- paste0("new_", substr(fname, start = 1, stop = nchar(fname) - 3), "_subset.tif")

# 組合輸出的完整路徑
output_path <- file.path(base_dir, output_filename)

# 寫入 GeoTIFF 檔案
writeRaster(sub.fog, filename = output_path, overwrite = TRUE)

# --- 繪圖檢視並同時輸出 PNG 檔案 ---
# 建立 PNG 輸出的完整路徑
output_png_path <- file.path(base_dir, paste0("new_", substr(fname, start = 1, stop = nchar(fname) - 3), "_subset.png"))

# 繪製並儲存 PNG 檔案
cat(paste0("\n--> 正在繪製並儲存視覺化圖檔至：", output_png_path, "\n"))
png(filename = output_png_path, width = 800, height = 600)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(sub.fog, main = paste0(substr(fname, start = 1, stop = nchar(fname) - 3), "_change."))
dev.off()

# 顯示成功訊息
cat(paste0("\n處理完成！新檔案已儲存至：", output_path, "\n"))
cat(paste0("視覺化圖檔已儲存至：", output_png_path, "\n"))
