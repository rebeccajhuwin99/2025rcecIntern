# =======================================================
# 氣候變遷影響分析 (針對選定組合)
# 根據手動挑選的最佳參數組合，計算氣候變遷情境下的 Grass 數量
# =======================================================
library(terra)

# 設定隨機種子，確保每次模擬結果可重現
set.seed(42)

# --- Define the fun_ca_process function for simulation ---
# (此函式已從先前的腳本複製，用於後續的模擬)
fun_ca_process <- function(in_arr, it_step, tg_value, final_lu_map, rain_priority, fog_priority, rain_weight, priority_influence, threshold_method, threshold_value) {
  dim_xy <- dim(in_arr)
  nx <- dim_xy[1]
  ny <- dim_xy[2]
  
  mat <- in_arr
  mat[in_arr == tg_value] <- 1
  mat[in_arr != tg_value] <- 0
  temp_mat <- mat
  
  for (k in 1:it_step) {
    for (i in 1:nx) {
      for (j in 1:ny) {
        E <- i + 1; W <- i - 1; N <- j - 1; S <- j + 1
        if (E > nx) { E = 1 }
        if (W < 1) { W = nx }
        if (N < 1) { N = ny }
        if (S > ny) { S = 1 }
        
        numb_alive = mat[W, N] + mat[i, N] + mat[E, N] +
          mat[W, j] + mat[E, j] +
          mat[W, S] + mat[i, S] + mat[E, S]
        
        if (mat[i, j] == 1 && numb_alive < 2) { temp_mat[i, j] = 0 }
        if (mat[i, j] == 1 && numb_alive > 3) { temp_mat[i, j] = 0 }
        if (mat[i, j] == 1 && (numb_alive == 2 || numb_alive == 3)) { temp_mat[i, j] = 1 }
        if (mat[i, j] == 0 && numb_alive == 3) { temp_mat[i, j] = 1 }
      }
    }
    mat = temp_mat
  }
  
  ca_candidates_arr <- mat
  ca_candidates_arr[in_arr == tg_value] <- 0
  
  num_ca_candidates_raw <- sum(ca_candidates_arr == 1, na.rm = TRUE)
  if (num_ca_candidates_raw == 0) {
    return(in_arr)
  }
  
  normalize_raster <- function(r) {
    r_values <- values(r, na.rm = TRUE)
    if (length(r_values) == 0 || min(r_values) == max(r_values)) {
      return(r * 0)
    }
    min_val <- min(r_values)
    max_val <- max(r_values)
    return((r - min_val) / (max_val - min_val))
  }
  
  normalized_rain <- normalize_raster(rain_priority)
  normalized_fog <- normalize_raster(fog_priority)
  
  random_raster <- rast(matrix(runif(nx * ny), nrow = nx, ncol = ny),
                        extent = ext(final_lu_map), crs = crs(final_lu_map))
  
  priority_factors_combined <- (rain_weight * normalized_rain) + ((1 - rain_weight) * normalized_fog)
  
  combined_suitability <- (priority_influence * priority_factors_combined) + ((1 - priority_influence) * random_raster)
  
  ca_mask_rast <- rast(ca_candidates_arr, extent = ext(final_lu_map), crs = crs(final_lu_map))
  combined_suitability_masked <- mask(combined_suitability, ca_mask_rast, maskvalues = 0, inverse = FALSE, updatevalue = NA)
  
  suitability_values <- values(combined_suitability_masked, na.rm = TRUE)
  
  if (length(suitability_values) == 0) {
    return(in_arr)
  }
  
  current_threshold <- 0
  if (threshold_method == "fixed") {
    current_threshold <- threshold_value
  } else if (threshold_method == "quantile") {
    current_threshold <- quantile(suitability_values, probs = threshold_value, na.rm = TRUE)
  }
  
  cells_to_convert_rast <- ifel(combined_suitability_masked >= current_threshold, 1, NA)
  
  out_pixel_arr <- in_arr
  cells_to_convert_arr <- as.array(cells_to_convert_rast)[,,1]
  
  out_pixel_arr[ !is.na(cells_to_convert_arr) ] <- tg_value
  
  return(out_pixel_arr)
}

# --- 設定輸入檔案路徑 ---
lu_1718_path <- "C:/Users/user/Downloads/1718_229_198_tree20.tif"
rain_2018_path <- "C:/Users/user/Downloads/rain_2018_yearly_subset.tif"
fog_2018_path <- "C:/Users/user/Downloads/final_6mfog_2018wc_229_198.tif"
fog_dif_2018_path <- "C:/Users/user/Downloads/final_6mfogdif_2018wc_229_198.tif"

# 儲存結果的基準路徑
output_base_dir <- "D:/"

# 檢查資料夾是否存在，不存在則建立
if (!dir.exists(output_base_dir)) {
  dir.create(output_base_dir, recursive = TRUE)
  cat(paste("已建立輸出資料夾:", output_base_dir, "\n"))
}

# --- 載入資料 ---
lu_1718_map <- rast(lu_1718_path)
rain_pri <- rast(rain_2018_path)
fog_pri <- rast(fog_2018_path)
fog_dif_pri <- rast(fog_dif_2018_path)

# --- 確保所有 raster 空間對齊 ---
rain_pri <- resample(rain_pri, lu_1718_map, method='bilinear')
fog_pri <- resample(fog_pri, lu_1718_map, method='bilinear')
fog_dif_pri <- resample(fog_dif_pri, lu_1718_map, method='bilinear')

# --- 氣候變遷情境下的霧氣圖層 ---
fog_with_change <- fog_pri + fog_dif_pri

# --- **請在此處填入你挑選的 5 組參數組合** ---
selected_combinations <- data.frame(
  it_step = c(2, 1, 2),  # 替換成你挑選的 it_step 值
  rain_weight = c(0.5, 0.75, 0.25), # 替換成你挑選的 rain_weight 值
  priority_influence = c(0.75, 0.75, 0.75), # 替換成你挑選的 priority_influence 值
  threshold_method = c("fixed", "fixed", "fixed"), # 通常保持為 "fixed"
  threshold_value = c(0.9, 0.9, 0.9) # 替換成你挑選的 threshold_value 值
)

# --- 準備儲存結果的資料框 ---
results_df_climate <- data.frame(
  it_step = integer(),
  rain_weight = numeric(),
  priority_influence = numeric(),
  threshold_value = numeric(),
  grass_count_sim_climate = numeric()
)

# --- 執行氣候變遷情境模擬 ---
cat("\n--- 開始對選定的 5 組組合進行氣候變遷情境模擬 ---\n")
target_value <- 6 # 假設最終轉換成草地

for (i in 1:nrow(selected_combinations)) {
  # 從資料框中取出當前組合的參數
  combo <- selected_combinations[i, ]
  
  cat(paste0("  模擬進度：", round(i / nrow(selected_combinations) * 100), "% - 正在處理第 ", i, " 組組合...\n"))
  
  # 執行「有氣候變遷」情境模擬
  simulated_lu_arr_climate <- tryCatch({
    fun_ca_process(
      in_arr = as.array(lu_1718_map)[,,1],
      it_step = combo$it_step,
      tg_value = target_value,
      final_lu_map = lu_1718_map,
      rain_priority = rain_pri,
      fog_priority = fog_with_change,
      rain_weight = combo$rain_weight,
      priority_influence = combo$priority_influence,
      threshold_method = combo$threshold_method,
      threshold_value = combo$threshold_value
    )
  }, error = function(e) {
    cat(paste("-> 氣候變遷模擬錯誤發生:", e$message, "\n"))
    return(as.array(lu_1718_map)[,,1])
  })
  
  # 計算模擬結果中的草地數量
  grass_count_climate <- sum(values(rast(simulated_lu_arr_climate)) == target_value, na.rm = TRUE)
  
  # 將結果儲存回資料框
  new_row <- data.frame(
    it_step = combo$it_step,
    rain_weight = combo$rain_weight,
    priority_influence = combo$priority_influence,
    threshold_value = combo$threshold_value,
    grass_count_sim_climate = grass_count_climate
  )
  results_df_climate <- rbind(results_df_climate, new_row)
}

cat("\n--- 氣候變遷情境模擬完成 ---\n")

# --- 顯示與儲存結果 ---
cat("\n--- 選定組合的氣候變遷模擬結果 ---\n")
print(results_df_climate)

# 將結果儲存到 CSV 檔案
output_path_csv <- file.path(output_base_dir, "selected_climate_impact_results.csv")
write.csv(results_df_climate, file = output_path_csv, row.names = FALSE)
cat(paste0("\n結果已儲存至:", output_path_csv, "\n"))