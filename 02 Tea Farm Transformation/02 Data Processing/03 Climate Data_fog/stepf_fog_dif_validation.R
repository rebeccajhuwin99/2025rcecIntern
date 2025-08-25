# 載入必要的套件
library("terra")
library("ggplot2")

# --- 設定您的絕對路徑基本資料夾 ---
# 請將這裡的路徑替換為您電腦中存放所有資料的實際路徑
base_dir <- "D:/20250818_中研院環變周報附件/03 氣候變遷 vs 無/cc 霧氣資料"

# --- 設定所有檔案的絕對路徑 ---
# 載入作為裁切與對齊基準的 lu.map
lu_map_path <- "C:/Users/user/Downloads/2024_229_198_tree20.tif"

# 1. 原始、全台灣範圍的 2018 年霧氣資料 (這是你驗證的起點)
fog_2018_path <- file.path(base_dir, "fog_2018_100m_wc.nc")
# 2. 原始、全台灣範圍的 PGW 霧氣資料 (這是你驗證的目標)
pgw_fog_path <- file.path(base_dir, "fog_2018_100m_pgw.nc")
# 3. 原始、全台灣範圍的差異資料 (dif) (這是你想要校驗的檔案)
fog_dif_path <- file.path(base_dir, "fog_2018_100m_pgw_wc_dif.nc")

# --- 載入資料 ---
# 使用 terra::rast() 載入資料
lu.map <- rast(lu_map_path)
fog_2018_full <- rast(fog_2018_path)
pgw_fog_full <- rast(pgw_fog_path)
fog_dif_full <- rast(fog_dif_path)

# --- **核心驗證步驟：在全範圍進行運算** ---
cat("--> Starting the validation process: calculating PGW based on your dif file...\n")
# 步驟一：根據 2018 資料與你的 dif 檔案，計算出一個「模擬的 PGW」
calculated_pgw_full <- fog_2018_full + fog_dif_full
cat("--> Successfully calculated the PGW simulation.\n")

# 步驟二：將「模擬的 PGW」與原始的 PGW 進行比較，計算兩者的差異
difference_full <- calculated_pgw_full - pgw_fog_full
cat("--> Calculated the difference between the simulated and original PGW.\n")

# --- 裁切與對齊最終的差異圖層 ---
# 定義一個函數來簡化裁切和對齊的過程
align_and_crop <- function(raster_full, target_raster) {
  # 重新投影到 lu.map 的坐標系統
  projected_raster <- project(raster_full, target_raster)
  # 重新取樣以對齊解析度和網格
  resampled_raster <- resample(projected_raster, target_raster, method = 'bilinear')
  # 使用 lu.map 的範圍進行裁切
  cropped_raster <- crop(resampled_raster, ext(target_raster))
  return(cropped_raster)
}

# 使用函數處理最終的差異圖層
cat("--> Aligning and cropping the final difference layer to the lu.map extent and resolution...\n")
difference <- align_and_crop(difference_full, lu.map)

# --- 輸出統計摘要 ---
# 取得差異圖層的統計值
diff_values <- values(difference)
diff_summary <- summary(diff_values)

cat("\n--- Summary Statistics of the Difference Layer ---\n")
print(diff_summary)

# --- 繪圖檢視 ---
# 將差異圖層轉換為 data.frame 以便用 ggplot2 繪圖
difference_df <- as.data.frame(difference, xy=TRUE)
names(difference_df)[3] <- "Difference"

# 使用 ggplot2 繪製差異圖
cat("\n--> Plotting the difference map. Values closer to 0 indicate higher accuracy...\n")
p <- ggplot(data = difference_df, aes(x = x, y = y, fill = Difference)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") + # 使用 viridis 調色盤，讓視覺效果更好
  labs(
    title = "Difference between Calculated PGW and Original PGW",
    x = NULL,
    y = NULL,
    fill = "Difference Value"
  ) +
  coord_equal() +
  theme_minimal()

print(p)

# 顯示成功訊息
cat(paste0("\nValidation complete! Please check the plot and summary statistics.\n"))
