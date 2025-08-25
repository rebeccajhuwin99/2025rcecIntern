# ============================================
# CA + 環境適宜度（雨/霧）氣候變遷比較版
# 目的：對最佳參數組合進行有無氣候變遷的 LU 模擬，並記錄 Grass 像元數
# 輸出：總共 8 個 tif、8 個 png，以及一個包含 Grass 像元數的 CSV
# 新增功能：將模擬後新生成的 Grass 像元標記為 7-pink
# 更新功能：計算 7-pink 的像元數量並記錄到 CSV
# ============================================
library(terra)

# 僅在腳本開始時設定一次隨機種子，確保後續模擬的隨機性不同
set.seed(42)

# ---------- 輸入檔案路徑 ----------
# 請確保以下路徑正確無誤
lu_1718_path      <- "C:/Users/user/Downloads/1718_229_198_tree20.tif"
lu_2024_path      <- "C:/Users/user/Downloads/2021_229_198_tree20.tif"
rain_2018_path    <- "C:/Users/user/Downloads/rain_2018_yearly_subset.tif"
fog_2018_path     <- "C:/Users/user/Downloads/final_6mfog_2018wc_229_198.tif"
fog_dif_2018_path <- "C:/Users/user/Downloads/final_6mfogdif_2018wc_229_198.tif"


# ---------- 輸出資料夾 ----------
out_dir_path    <- "D:/"
if (!dir.exists(out_dir_path)) dir.create(out_dir_path, recursive = TRUE)

# ---------- 顏色（1~6），新增 7-pink ----------
cols_1to6 <- c(
  "darkgreen",    # 1-Forest
  "brown",        # 2-Built-up
  "lightblue",    # 3-Water
  "yellow",       # 4-Agri
  "gray",         # 5-Unkn
  "yellowgreen"   # 6-Grass
)
names(cols_1to6) <- 1:6

# 針對新增的 7 號類別，準備新的顏色向量
cols_1to7 <- c(cols_1to6, "pink")
names(cols_1to7) <- 1:7

# 為 plot() 設定完整的顏色和圖例標籤
plot_levels_1to6 <- data.frame(
  ID = 1:6,
  class = c("1-Forest", "2-Built-up", "3-Water", "4-Agri", "5-Unkn", "6-Grass")
)
# 新增 7 號類別的圖例標籤
plot_levels_1to7 <- data.frame(
  ID = 1:7,
  class = c("1-Forest", "2-Built-up", "3-Water", "4-Agri", "5-Unkn", "6-Grass", "7-New Grass")
)


# ============================================
# 輔助函式 (Function)
# ============================================
# 對齊到 LU
align_to_lu <- function(x, template) {
  if (!compareGeom(x, template, stopOnError = FALSE)) {
    if (!same.crs(x, template)) {
      x <- project(x, template, method = "bilinear")
    }
    if (!compareGeom(x, template, stopOnError = FALSE)) {
      x <- resample(x, template, method = "bilinear")
    }
  }
  return(x)
}

# 正規化到 [0,1]
norm01 <- function(r) {
  rmin <- global(r, "min", na.rm = TRUE)[1,1]
  rmax <- global(r, "max", na.rm = TRUE)[1,1]
  if (is.na(rmin) || is.na(rmax) || rmax == rmin) {
    return(r*0)
  } else {
    (r - rmin) / (rmax - rmin)
  }
}

# 核心 CA 運算函式
# 將 fun_ca_process 整合到這裡，確保單一檔案的完整性
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

# ============================================
# 主運算與輸出流程
# ============================================

# ---------- 您提供的 goodfit 參數組合 ----------
best_combinations <- data.frame(
  run_id = 1:2,
  it_step = c(2, 1),
  rain_weight = c(0.5, 0.75),
  priority_influence = c(0.75, 0.75),
  threshold_value = c(0.9, 0.9)
)
best_combinations$threshold_method <- "fixed" # 根據您的輸入，門檻法都是 fixed


# 讀取所有必要資料
lu_1718_map  <- rast(lu_1718_path)
lu_2024_map  <- rast(lu_2024_path)
rain_pri     <- rast(rain_2018_path)
fog_pri      <- rast(fog_2018_path)
fog_dif_pri  <- rast(fog_dif_2018_path)

# 確保所有 raster 空間對齊
rain_pri    <- align_to_lu(rain_pri, lu_1718_map)
fog_pri     <- align_to_lu(fog_pri, lu_1718_map)
fog_dif_pri <- align_to_lu(fog_dif_pri, lu_1718_map)
lu_2024_map <- align_to_lu(lu_2024_map, lu_1718_map)

# 將 lu_1718_map 轉為 array，方便後續的比較運算
lu_1718_arr <- as.array(lu_1718_map)[,,1]

# 準備 Grass 像元數的記錄資料框
# 新增一個欄位來記錄新增的 Grass 像元數
grass_counts_df <- data.frame(
  run_id = integer(),
  description = character(),
  initial_grass_count = integer(),
  simulated_grass_count = integer(),
  new_grass_pixels = integer()
)

# 記錄參考資料的 Grass 像元數
initial_grass_count_ref <- sum(lu_1718_arr == 6, na.rm = TRUE)
final_grass_count_ref <- sum(values(lu_2024_map) == 6, na.rm = TRUE)

# 將參考數據加入資料框，新欄位設為 NA
grass_counts_df <- rbind(grass_counts_df, data.frame(
  run_id = NA,
  description = "Reference_2018",
  initial_grass_count = initial_grass_count_ref,
  simulated_grass_count = NA,
  new_grass_pixels = NA
))
grass_counts_df <- rbind(grass_counts_df, data.frame(
  run_id = NA,
  description = "Reference_2021",
  initial_grass_count = initial_grass_count_ref,
  simulated_grass_count = final_grass_count_ref,
  new_grass_pixels = NA
))


# 針對每一組最佳組合進行模擬
for (i in 1:nrow(best_combinations)) {
  combo <- best_combinations[i, ]
  
  # ---------- 模擬 1：無氣候變遷 (No Climate Change) ----------
  cat(paste0("\n=== 開始模擬：Run ", combo$run_id, " - 無氣候變遷 ===\n"))
  
  sim_arr_no_cc <- fun_ca_process(
    in_arr = lu_1718_arr,
    it_step = combo$it_step,
    tg_value = 6,
    final_lu_map = lu_1718_map,
    rain_priority = rain_pri,
    fog_priority = fog_pri,
    rain_weight = combo$rain_weight,
    priority_influence = combo$priority_influence,
    threshold_method = combo$threshold_method,
    threshold_value = combo$threshold_value
  )
  sim_rast_no_cc <- rast(sim_arr_no_cc, extent = ext(lu_1718_map), crs = crs(lu_1718_map))
  
  # 存檔
  base_tag_no_cc <- paste0("no_cc_run", combo$run_id,
                           "_it", combo$it_step,
                           "_rw", sprintf("%.2f", combo$rain_weight),
                           "_pi", sprintf("%.2f", combo$priority_influence),
                           "_thr", sprintf("%.2f", combo$threshold_value))
  tif_path_no_cc <- file.path(out_dir_path, paste0(base_tag_no_cc, ".tif"))
  png_path_no_cc <- file.path(out_dir_path, paste0(base_tag_no_cc, ".png"))
  
  writeRaster(sim_rast_no_cc, tif_path_no_cc, overwrite = TRUE)
  png(png_path_no_cc, width = 1600, height = 1400, res = 200)
  plot(sim_rast_no_cc, col = cols_1to6, main = paste0("LU Simulation: No Climate Change (Run ", combo$run_id, ")"))
  dev.off()
  
  # 記錄 Grass 像元數
  grass_count_no_cc <- sum(values(sim_rast_no_cc) == 6, na.rm = TRUE)
  
  # 新增功能：找出新長出的 Grass，並將其標記為 7-pink
  # 找出初始不為 Grass 且模擬後變成 Grass 的像元
  new_grass_pixels_no_cc <- (lu_1718_arr != 6) & (sim_arr_no_cc == 6)
  
  # 記錄新長出的 Grass 像元數
  new_grass_count_no_cc <- sum(new_grass_pixels_no_cc, na.rm = TRUE)
  
  # 複製模擬結果，並將新長出的 Grass 標記為 7
  sim_arr_no_cc_pink <- sim_arr_no_cc
  sim_arr_no_cc_pink[new_grass_pixels_no_cc] <- 7
  
  # 轉回 raster
  sim_rast_no_cc_pink <- rast(sim_arr_no_cc_pink, extent = ext(lu_1718_map), crs = crs(lu_1718_map))
  
  # 存檔，檔名加上 _pink
  tif_path_no_cc_pink <- file.path(out_dir_path, paste0(base_tag_no_cc, "_pink.tif"))
  png_path_no_cc_pink <- file.path(out_dir_path, paste0(base_tag_no_cc, "_pink.png"))
  
  writeRaster(sim_rast_no_cc_pink, tif_path_no_cc_pink, overwrite = TRUE)
  png(png_path_no_cc_pink, width = 1600, height = 1400, res = 200)
  plot(sim_rast_no_cc_pink, col = cols_1to7, main = paste0("LU Simulation: No Climate Change (Run ", combo$run_id, ") - New Grass Highlighted"))
  dev.off()
  
  # 更新資料框，加入新長出的 Grass 數量
  grass_counts_df <- rbind(grass_counts_df, data.frame(
    run_id = combo$run_id,
    description = "No_Climate_Change",
    initial_grass_count = initial_grass_count_ref,
    simulated_grass_count = grass_count_no_cc,
    new_grass_pixels = new_grass_count_no_cc
  ))
  
  
  # ---------- 模擬 2：有氣候變遷 (With Climate Change) ----------
  cat(paste0("\n=== 開始模擬：Run ", combo$run_id, " - 有氣候變遷 ===\n"))
  
  # 計算有氣候變遷的霧氣因子（2018 + 變化量）
  fog_plus_dif_pri <- fog_pri + fog_dif_pri
  
  sim_arr_with_cc <- fun_ca_process(
    in_arr = lu_1718_arr,
    it_step = combo$it_step,
    tg_value = 6,
    final_lu_map = lu_1718_map,
    rain_priority = rain_pri,
    fog_priority = fog_plus_dif_pri, # 使用有變化的霧氣資料
    rain_weight = combo$rain_weight,
    priority_influence = combo$priority_influence,
    threshold_method = combo$threshold_method,
    threshold_value = combo$threshold_value
  )
  sim_rast_with_cc <- rast(sim_arr_with_cc, extent = ext(lu_1718_map), crs = crs(lu_1718_map))
  
  # 存檔
  base_tag_with_cc <- paste0("with_cc_run", combo$run_id,
                             "_it", combo$it_step,
                             "_rw", sprintf("%.2f", combo$rain_weight),
                             "_pi", sprintf("%.2f", combo$priority_influence),
                             "_thr", sprintf("%.2f", combo$threshold_value))
  tif_path_with_cc <- file.path(out_dir_path, paste0(base_tag_with_cc, ".tif"))
  png_path_with_cc <- file.path(out_dir_path, paste0(base_tag_with_cc, ".png"))
  
  writeRaster(sim_rast_with_cc, tif_path_with_cc, overwrite = TRUE)
  png(png_path_with_cc, width = 1600, height = 1400, res = 200)
  plot(sim_rast_with_cc, col = cols_1to6, main = paste0("LU Simulation: With Climate Change (Run ", combo$run_id, ")"))
  dev.off()
  
  # 記錄 Grass 像元數
  grass_count_with_cc <- sum(values(sim_rast_with_cc) == 6, na.rm = TRUE)
  
  # 新增功能：找出新長出的 Grass，並將其標記為 7-pink
  new_grass_pixels_with_cc <- (lu_1718_arr != 6) & (sim_arr_with_cc == 6)
  
  # 記錄新長出的 Grass 像元數
  new_grass_count_with_cc <- sum(new_grass_pixels_with_cc, na.rm = TRUE)
  
  sim_arr_with_cc_pink <- sim_arr_with_cc
  sim_arr_with_cc_pink[new_grass_pixels_with_cc] <- 7
  
  sim_rast_with_cc_pink <- rast(sim_arr_with_cc_pink, extent = ext(lu_1718_map), crs = crs(lu_1718_map))
  
  # 存檔，檔名加上 _pink
  tif_path_with_cc_pink <- file.path(out_dir_path, paste0(base_tag_with_cc, "_pink.tif"))
  png_path_with_cc_pink <- file.path(out_dir_path, paste0(base_tag_with_cc, "_pink.png"))
  
  writeRaster(sim_rast_with_cc_pink, tif_path_with_cc_pink, overwrite = TRUE)
  png(png_path_with_cc_pink, width = 1600, height = 1400, res = 200)
  plot(sim_rast_with_cc_pink, col = cols_1to7, main = paste0("LU Simulation: With Climate Change (Run ", combo$run_id, ") - New Grass Highlighted"))
  dev.off()
  
  # 更新資料框，加入新長出的 Grass 數量
  grass_counts_df <- rbind(grass_counts_df, data.frame(
    run_id = combo$run_id,
    description = "With_Climate_Change",
    initial_grass_count = initial_grass_count_ref,
    simulated_grass_count = grass_count_with_cc,
    new_grass_pixels = new_grass_count_with_cc
  ))
}

# 輸出最終的 Grass 像元數比較表
output_csv_path <- file.path(out_dir_path, "grass_growth_comparison.csv")
write.csv(grass_counts_df, output_csv_path, row.names = FALSE)
cat(paste0("\n所有模擬結果與 Grass 像元數已儲存至:", output_csv_path, "\n"))
