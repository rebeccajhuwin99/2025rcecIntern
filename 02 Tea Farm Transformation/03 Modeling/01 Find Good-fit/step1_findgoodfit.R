# =======================================================
# CA 參數校準
# 執行所有參數組合的校準，並輸出排序後的結果
# 輸出：一個包含所有校準結果的 CSV 檔案
# =======================================================
library(terra)

# 設定隨機種子，確保每次模擬結果可重現
set.seed(42)

# --- 定義計算 Kappa 係數的輔助函式 ---
# 這個函式會將兩個 raster 轉為矩陣，然後計算它們之間的 Kappa 係數。
# Kappa 係數是用來衡量兩個分類結果一致性的指標，數值越接近 1 代表一致性越高。
calculate_kappa <- function(simulated_raster, reference_raster) {
  # 確保兩個 raster 空間對齊
  if (!compareGeom(simulated_raster, reference_raster)) {
    reference_raster <- resample(reference_raster, simulated_raster, method = 'near')
  }
  
  # 將 raster 轉換為向量以計算混淆矩陣
  sim_vals <- values(simulated_raster, na.rm=TRUE)
  ref_vals <- values(reference_raster, na.rm=TRUE)
  
  # 移除 NA 值以確保向量長度相同
  common_indices <- which(!is.na(sim_vals) & !is.na(ref_vals))
  sim_vals <- sim_vals[common_indices]
  ref_vals <- ref_vals[common_indices]
  
  # 建立混淆矩陣
  # table() 會統計兩個向量中各類別的交叉數量
  conf_matrix <- table(sim_vals, ref_vals)
  
  # 計算總樣本數
  total_samples <- sum(conf_matrix)
  
  # 計算觀察到的一致性 (Po)
  # Po 是混淆矩陣對角線元素（模擬與參考一致）的總和，除以總樣本數
  po <- sum(diag(conf_matrix)) / total_samples
  
  # 計算期望的一致性 (Pe)
  # Pe 是透過每個類別的邊際總和計算的
  row_sums <- rowSums(conf_matrix)
  col_sums <- colSums(conf_matrix)
  pe <- sum(row_sums * col_sums) / (total_samples^2)
  
  # 計算 Kappa 係數
  if (pe == 1) { # 避免分母為 0 的情況
    kappa_value <- 1
  } else {
    kappa_value <- (po - pe) / (1 - pe)
  }
  
  return(kappa_value)
}

# --- Define the fun_ca_process function for calibration ---
# 這個版本接受更多的參數，並移除了強制從林地(1)轉換的限制
fun_ca_process <- function(in_arr, it_step, tg_value, final_lu_map, rain_priority, fog_priority, rain_weight, priority_influence, threshold_method, threshold_value) {
  dim_xy <- dim(in_arr)
  nx <- dim_xy[1]
  ny <- dim_xy[2]
  
  # --- Step 1: Run CA to get a mask of potentially activated cells ---
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
  ca_candidates_arr[in_arr == tg_value] <- 0 # Candidates cannot be the target type already
  
  num_ca_candidates_raw <- sum(ca_candidates_arr == 1, na.rm = TRUE)
  if (num_ca_candidates_raw == 0) {
    return(in_arr)
  }
  
  # --- Step 2: Incorporate Rainfall and Fog as Priorities ---
  normalize_raster <- function(r) {
    r_values <- values(r, na.rm = TRUE)
    if (length(r_values) == 0 || min(r_values) == max(r_values)) {
      return(r * 0)
    }
    min_val <- min(r_values)
    max_val <- max(r_values)
    return((r - min_val) / (max_val - min_val))
  }
  
  # 正規化降雨和霧氣的圖層
  normalized_rain <- normalize_raster(rain_priority)
  normalized_fog <- normalize_raster(fog_priority)
  
  # 建立一個隨機雜訊圖層，範圍為 0 到 1
  random_raster <- rast(matrix(runif(nx * ny), nrow = nx, ncol = ny),
                        extent = ext(final_lu_map), crs = crs(final_lu_map))
  
  # 根據 priority_influence 參數來結合優先因素和隨機性
  priority_factors_combined <- (rain_weight * normalized_rain) + ((1 - rain_weight) * normalized_fog)
  
  combined_suitability <- (priority_influence * priority_factors_combined) + ((1 - priority_influence) * random_raster)
  
  # --- Step 3: Apply the CA mask to the combined suitability raster ---
  ca_mask_rast <- rast(ca_candidates_arr, extent = ext(final_lu_map), crs = crs(final_lu_map))
  combined_suitability_masked <- mask(combined_suitability, ca_mask_rast, maskvalues = 0, inverse = FALSE, updatevalue = NA)
  
  suitability_values <- values(combined_suitability_masked, na.rm = TRUE)
  
  if (length(suitability_values) == 0) {
    return(in_arr)
  }
  
  # --- Step 4: Calculate threshold and select cells to convert based on method ---
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
lu_2024_path <- "C:/Users/user/Downloads/2021_229_198_tree20.tif"
rain_2018_path <- "C:/Users/user/Downloads/rain_2018_yearly_subset.tif"
fog_2018_path <- "C:/Users/user/Downloads/final_6mfog_2018wc_229_198.tif"
# 新增 fog_dif 檔案路徑
# fog_dif_2018_path <- "C:/Users/user/Downloads/final_6mfogdif_2018wc_229_198.tif"

# 儲存所有模擬結果的基準路徑
output_base_dir <- "D:/"

# 檢查資料夾是否存在，不存在則建立
if (!dir.exists(output_base_dir)) {
  dir.create(output_base_dir, recursive = TRUE)
  cat(paste("已建立輸出資料夾:", output_base_dir, "\n"))
}

# --- 載入資料 ---
lu_1718_map <- rast(lu_1718_path)
lu_2024_map <- rast(lu_2024_path)
rain_pri <- rast(rain_2018_path)
fog_pri <- rast(fog_2018_path)
# fog_dif_pri <- rast(fog_dif_2018_path)

# --- 確保所有 raster 空間對齊 ---
rain_pri <- resample(rain_pri, lu_1718_map, method='bilinear')
fog_pri <- resample(fog_pri, lu_1718_map, method='bilinear')
# fog_dif_pri <- resample(fog_dif_pri, lu_1718_map, method='bilinear')
lu_2024_map <- resample(lu_2024_map, lu_1718_map, method='near')


# --- **參數網格設定** ---
it_steps <- 1:5 # 測試 1 到 5 個步數
rain_weights <- seq(0, 1, by = 0.25) # 測試 0 到 1 的所有權重
priority_influences <- seq(0, 1, by = 0.25) # 測試 0 到 1 的所有影響力
fixed_thresholds <- seq(0, 0.9, by = 0.1) # 測試 0 到 0.9 的所有固定閾值
threshold_methods <- c("fixed")

# --- 準備儲存結果的資料框 ---
results_df <- data.frame(
  it_step = integer(),
  rain_weight = numeric(),
  priority_influence = numeric(),
  threshold_method = character(),
  threshold_value = numeric(),
  kappa_value = numeric(),
  grass_count_sim_normal = numeric(),
  grass_count_reference_1718 = numeric(),
  grass_count_reference_2024 = numeric()
)

# --- 模擬參數 ---
target_value <- 6 # 假設最終轉換成草地
ref_grass_count_1718 <- sum(values(lu_1718_map) == target_value, na.rm = TRUE)
ref_grass_count_2024 <- sum(values(lu_2024_map) == target_value, na.rm = TRUE)


# --- 開始進行 CA 參數校準實驗 ---
cat("\n--- 開始進行 CA 參數校準實驗（只執行無氣候變遷情境）---\n")
total_combinations <- length(it_steps) * length(rain_weights) * length(priority_influences) * length(fixed_thresholds)
current_combo <- 0

# --- 進行多層嵌套迴圈以測試所有參數組合 ---
for (it_step in it_steps) {
  for (weight in rain_weights) {
    for (influence in priority_influences) {
      for (method in threshold_methods) {
        thresholds_to_test <- fixed_thresholds
        
        for (thresh_val in thresholds_to_test) {
          current_combo <- current_combo + 1
          
          # --- 執行「無氣候變遷」情境模擬 ---
          simulated_lu_arr_normal <- tryCatch({
            fun_ca_process(
              in_arr = as.array(lu_1718_map)[,,1],
              it_step = it_step,
              tg_value = target_value,
              final_lu_map = lu_1718_map,
              rain_priority = rain_pri,
              fog_priority = fog_pri,
              rain_weight = weight,
              priority_influence = influence,
              threshold_method = method,
              threshold_value = thresh_val
            )
          }, error = function(e) {
            cat(paste(" -> 錯誤發生:", e$message, "\n"))
            return(as.array(lu_1718_map)[,,1])
          })
          
          simulated_lu_raster_normal <- rast(simulated_lu_arr_normal, extent = ext(lu_1718_map), crs = crs(lu_1718_map))
          kappa_result <- calculate_kappa(simulated_lu_raster_normal, lu_2024_map)
          grass_count_normal <- sum(values(simulated_lu_raster_normal) == target_value, na.rm = TRUE)
          
          # 將結果儲存到資料框
          new_row <- data.frame(
            it_step = it_step,
            rain_weight = weight,
            priority_influence = influence,
            threshold_method = method,
            threshold_value = thresh_val,
            kappa_value = kappa_result,
            grass_count_sim_normal = grass_count_normal,
            grass_count_reference_1718 = ref_grass_count_1718,
            grass_count_reference_2024 = ref_grass_count_2024
          )
          results_df <- rbind(results_df, new_row)
          
          # 新增：顯示即時進度與 Kappa 值
          cat(paste0(" 進度：", round(current_combo / total_combinations * 100), "% - 組合: Step=", it_step, ", Weight=", weight, ", Influence=", influence, ", Value=", thresh_val, " -> Kappa: ", round(kappa_result, 4), "\n"))
        }
      }
    }
  }
}
cat("\n--- 參數校準實驗結束 ---\n")

# --- 排序並顯示結果 ---
results_df_sorted <- results_df[order(-results_df$kappa_value), ]
cat("\n--- 依據 Kappa 係數排序後的前 10 名組合 ---\n")
print(head(results_df_sorted, 10))

# --- 將結果儲存為 CSV 檔案以供後續分析 ---
output_path_csv <- file.path(output_base_dir, "ca_calibration_results_sorted.csv")
write.csv(results_df_sorted, file = output_path_csv, row.names = FALSE)
cat(paste0("\n所有結果已儲存至:", output_path_csv, "\n"))
