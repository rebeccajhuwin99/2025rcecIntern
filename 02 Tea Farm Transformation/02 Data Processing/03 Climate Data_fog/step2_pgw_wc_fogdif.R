# 載入必要的套件
library(ncdf4)
library(terra) # 主要用於 raster 處理
# akima 套件在此步驟不需要，因為我們直接處理 SpatRaster 對象

# --- 1. 定義 100m 解析度的 WC 和 PGW 檔案路徑 ---
# 請務必將這些路徑替換為您實際的檔案位置
# 這些檔案應該是您之前運行第一個程式碼段（將原始 NetCDF 重新網格化到 100m）後生成的。
wc_100m_file <- "D:/TEST/fog_2018_100m_wc.nc" # 假設這是您的 2019 WC 檔案
pgw_100m_file <- "D:/TEST/fog_2018_100m_pgw.nc" # 假設這是您的 2019 PGW 檔案

# --- 2. 讀取 100m 的 WC 和 PGW SpatRaster ---
message(paste("正在讀取 WC 檔案:", wc_100m_file))
fog_100m_wc <- rast(wc_100m_file)
message(paste("正在讀取 PGW 檔案:", pgw_100m_file))
fog_100m_pgw <- rast(pgw_100m_file)

# --- 3. 檢查檔案的一致性 (範圍和 CRS) ---
# 為了確保計算正確，PGW 和 WC 資料必須具有相同的地理範圍和坐標參考系統。
# 如果您是使用相同的原始程式碼生成這兩個 100m 檔案，它們應該會一致。
if (ext(fog_100m_wc) != ext(fog_100m_pgw) || crs(fog_100m_wc) != crs(fog_100m_pgw)) {
  stop("錯誤：WC 和 PGW 檔案的範圍或坐標參考系統 (CRS) 不一致。\n請檢查它們是如何生成的，確保它們是完全對齊的。")
} else {
  message("WC 和 PGW 檔案的範圍和 CRS 一致，可以進行差異計算。")
}

# --- 4. 計算差異 (PGW - WC) ---
# 這將生成一個與 fog_100m_wc 具有相同 100m 解析度和範圍的 SpatRaster。
message("正在計算 PGW 和 WC 的差異...")
fog_dif_100m <- fog_100m_pgw - fog_100m_wc

# --- 5. 輸出結果為新的 NetCDF 檔案 ---
# 輸出檔案名稱確認為 fog_2019_100m_pgw_wc_dif.nc
output_dif_file <- "D:/TEST/fog_2018_100m_pgw_wc_dif.nc"

message(paste("正在將差異資料寫入:", output_dif_file))
writeCDF(fog_dif_100m, output_dif_file, varname = "yearly_fog_ratio_diff", overwrite = TRUE)

message(paste("成功生成 100m 解析度的 2019 年霧頻率差異資料:", output_dif_file))

# --- 6. (可選) 繪製結果 ---
plot(fog_dif_100m, main="Change in Annual Fog Frequency in 2018")