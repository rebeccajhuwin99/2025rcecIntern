# =======================================================
# CA Simulation Validation: 2013 → 2023
# - For each shop point: crop 750m window
# - Run CA simulation (10 generations)
# - Compare with observed 2023
# - Output confusion map PNG + metrics CSV
# =======================================================
library(terra)
library(dplyr)

# ============= Conway CA（單步、累積活格可關） =============
fun_ca_process <- function(in_arr, it_step = 1, tg_value = 1, accumulate = TRUE) {
  if (is.null(in_arr)) stop("in_arr is NULL")
  if (length(in_arr) == 0) stop("in_arr has length 0")
  if (is.null(dim(in_arr))) stop("in_arr has no dim; expected a matrix")
  nx <- nrow(in_arr); ny <- ncol(in_arr)
  out_arr <- matrix(0L, nrow = nx, ncol = ny)
  tmp_arr <- out_arr
  
  mat <- in_arr
  mat[in_arr == tg_value] <- 1L
  mat[in_arr != tg_value] <- 0L
  temp_mat <- mat
  
  for (k in 1:it_step) {
    for (i in 1:nx) {
      for (j in 1:ny) {
        E <- ifelse(i+1 == nx+1, 1, i+1)
        W <- ifelse(i-1 == 0,    nx, i-1)
        N <- ifelse(j-1 == 0,    ny, j-1)
        S <- ifelse(j+1 == ny+1, 1,  j+1)
        numb_alive <- mat[W,N] + mat[i,N] + mat[E,N] +
                      mat[W,j]            + mat[E,j] +
                      mat[W,S] + mat[i,S] + mat[E,S]
        if (mat[i,j]==1 && numb_alive<2)              temp_mat[i,j] <- 0L
        if (mat[i,j]==1 && numb_alive>3)              temp_mat[i,j] <- 0L
        if (mat[i,j]==1 && (numb_alive==2 || numb_alive==3)) temp_mat[i,j] <- 1L
        if (mat[i,j]==0 && numb_alive==3)             temp_mat[i,j] <- 1L
      }
    }
    mat <- temp_mat
    if (accumulate) for (ii in 1:nx) for (jj in 1:ny) if (mat[ii,jj] == 1L) tmp_arr[ii,jj] <- 1L
  }
  if (accumulate) list(out_arr = tmp_arr) else list(out_arr = mat)
}

# ================= 檔案路徑 =================
tiff_2013_path <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/lu_maps/2013_all_grid_wgs84.tif"
tiff_2023_path <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/lu_maps/2023_all_grid_wgs84.tif"
shops_150_path <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/csv/Changhua_shops_group_0.0015_df_out.csv"
shops_500_path <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/csv/Changhua_shops_group_0.005_df_out.csv"
output_folder  <- "C:/Users/cwiev/Desktop/intern/august 13/ca_verification_2013to2023"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# ================= 載入資料 =================
r2013_all <- rast(tiff_2013_path)
r2023_all <- rast(tiff_2023_path)

shops_150 <- read.csv(shops_150_path)
shops_500 <- read.csv(shops_500_path)
shop_data <- bind_rows(shops_150, shops_500) %>%
  distinct(name, .keep_all = TRUE)
if (!"name" %in% names(shop_data)) shop_data$name <- paste0("pt_", seq_len(nrow(shop_data)))

# ================= 模擬參數 =================
n_years       <- 10          # 2013 → 2023 共 10 代
radius_m      <- 750
ALIVE_CLASSES <- c(100)      # 若 50 也算活格：c(100, 50)
set.seed(42); terraOptions(progress = 0)

# ================= 欄位名相容處理 =================
has_lonx <- "lon.x" %in% names(shop_data); has_latx <- "lat.x" %in% names(shop_data)
has_lon  <- "lon"   %in% names(shop_data); has_lat  <- "lat"   %in% names(shop_data)
stopifnot((has_lonx && has_latx) || (has_lon && has_lat))
get_lon <- function(df, idx) if (has_lonx) df$lon.x[idx] else df$lon[idx]
get_lat <- function(df, idx) if (has_latx) df$lat.x[idx] else df$lat[idx]

# ================= 主迴圈 =================
metrics_rows <- list()

for (i in 1:nrow(shop_data)) {
  shop_lon <- suppressWarnings(as.numeric(get_lon(shop_data, i)))
  shop_lat <- suppressWarnings(as.numeric(get_lat(shop_data, i)))
  if (is.na(shop_lon) || is.na(shop_lat)) {
    cat("⚠️ NA coords, skip. name=", shop_data$name[i], "\n"); next
  }
  
  # 裁切 750m 視窗
  delta_deg <- radius_m / 111000
  ext_crop <- ext(shop_lon - delta_deg, shop_lon + delta_deg,
                  shop_lat - delta_deg, shop_lat + delta_deg)
  r2013 <- try(crop(r2013_all, ext_crop), silent = TRUE)
  r2023 <- try(crop(r2023_all, ext_crop), silent = TRUE)
  if (inherits(r2013, "try-error") || inherits(r2023, "try-error") ||
      ncell(r2013)==0 || ncell(r2023)==0) {
    cat("⚠️ Crop failed/empty, skip. name=", shop_data$name[i], "\n"); next
  }
  cat("👉 Running:", shop_data$name[i], " (", nrow(r2013), "x", ncol(r2013), ")\n")
  
  # —— 2013 重編碼 → 初始二值 ——
  ca0 <- values(r2013)
  ca0[ca0 == 2] <- 100
  ca0[ca0 == 5] <- 50
  ca0[!(ca0 %in% c(100, 50))] <- 0
  nX <- nrow(r2013); nY <- ncol(r2013)
  ca0_mat <- matrix(ca0, nrow = nX, ncol = nY, byrow = TRUE)
  ini_bin <- matrix(as.integer(ca0_mat %in% ALIVE_CLASSES), nrow = nX, ncol = nY, byrow = TRUE)
  
  # —— 區帶與機率 —— 
  coords   <- crds(r2013)
  dists_m  <- sqrt((coords[,1] - shop_lon)^2 + (coords[,2] - shop_lat)^2) * 111000
  zone_cls <- ifelse(dists_m <= 150, 1, ifelse(dists_m <= 500, 2, 3))
  inv_score <- 1 / (dists_m + 1); inv_score <- inv_score / max(inv_score)
  
  # —— 初值 & 上限 ——
  current    <- ca0
  max_growth <- ca0
  max_growth[zone_cls == 1] <- pmin(100, ca0[zone_cls == 1] + 30)
  max_growth[zone_cls == 2] <- pmin(100, ca0[zone_cls == 2] + 20)
  
  # —— 模擬 n_years —— 
  curr_bin <- ini_bin
  for (t in 1:n_years) {
    ca_step  <- fun_ca_process(curr_bin, it_step = 1, tg_value = 1, accumulate = FALSE)
    cand_bin <- ca_step$out_arr
    cand_vec <- as.vector(t(cand_bin))
    
    max_t <- max_growth
    max_t[cand_vec == 0 | zone_cls == 3] <- current[cand_vec == 0 | zone_cls == 3]
    
    pickable <- (cand_vec == 1) & (zone_cls != 3)
    selected <- rep(FALSE, length(current))
    if (any(pickable)) selected[pickable] <- runif(sum(pickable)) < inv_score[pickable]
    
    growth <- rep(0, length(current))
    if (any(selected)) growth[selected] <- runif(sum(selected), 1, 100)
    
    proposed <- current + growth
    current  <- ifelse(proposed <= max_t, proposed, max_t)
    curr_bin <- matrix(as.integer(current > 0), nrow = nX, ncol = nY, byrow = TRUE)
  }
  
  # —— 2023 實際 → 二值 —— 
  real <- values(r2023)
  real[real == 2] <- 100
  real[real == 5] <- 50
  real[!(real %in% c(100, 50))] <- 0
  real_mat <- matrix(real, nrow = nX, ncol = nY, byrow = TRUE)
  real_bin <- ifelse(real_mat > 0, 1, 0)
  
  # —— 模擬結果二值 —— 
  sim_bin_vec <- ifelse(current > 0, 1, 0)
  sim_bin_mat <- matrix(sim_bin_vec, nrow = nX, ncol = nY, byrow = TRUE)
  
  # —— 混淆矩陣與指標 —— 
  compare <- matrix(NA_integer_, nrow = nX, ncol = nY)
  compare[sim_bin_mat == 1 & real_bin == 1] <- 1  # TP
  compare[sim_bin_mat == 1 & real_bin == 0] <- 2  # FP
  compare[sim_bin_mat == 0 & real_bin == 1] <- 3  # FN
  compare[sim_bin_mat == 0 & real_bin == 0] <- 0  # TN
  
  TP <- sum(compare == 1, na.rm = TRUE)
  FP <- sum(compare == 2, na.rm = TRUE)
  FN <- sum(compare == 3, na.rm = TRUE)
  
  precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  recall    <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  f1        <- ifelse((precision + recall) == 0, 0, 2 * precision * recall / (precision + recall))
  
  cat(sprintf("    📊 P=%.3f  R=%.3f  F1=%.3f    (TP=%d FP=%d FN=%d)\n", precision, recall, f1, TP, FP, FN))
  
  metrics_rows[[length(metrics_rows)+1]] <- data.frame(
    index = i,
    name  = as.character(shop_data$name[i]),
    lon   = shop_lon,
    lat   = shop_lat,
    P = precision, R = recall, F1 = f1,
    stringsAsFactors = FALSE
  )
  
  # —— 輸出比較圖 —— 
  result_r <- r2013; values(result_r) <- as.vector(t(compare))
  png(file.path(output_folder, sprintf("sim_plot_%03d_comparison.png", i)), width = 1000, height = 800)
  col_scheme <- c("white", "#A8C653", "#8bd2e8", "#e86d5f")  # TN, TP, FP, FN
  legend_labels <- c("TN","TP","FP","FN")
  plot(result_r, col = col_scheme,
       main = paste0("Comparison - ", shop_data$name[i],
                     "\nP=", round(precision,2),
                     "  R=", round(recall,2),
                     "  F1=", round(f1,2)),
       legend = FALSE)
  points(shop_lon, shop_lat, pch = 19, col = "black")
  symbols(shop_lon, shop_lat, circles = 150/111000, add = TRUE, inches = FALSE, fg = "blue", lty = 2)
  symbols(shop_lon, shop_lat, circles = 500/111000, add = TRUE, inches = FALSE, fg = "red",  lty = 2)
  legend("bottomleft", fill = col_scheme, legend = legend_labels, bty = "n")
  dev.off()
}

# 彙整輸出
metrics_df <- dplyr::bind_rows(metrics_rows)
if (nrow(metrics_df) > 0) {
  utils::write.csv(metrics_df, file.path(output_folder, "simulation_metrics_raw.csv"), row.names = FALSE)
  cat("📄 已輸出：", file.path(output_folder, "simulation_metrics_raw.csv"), "\n")
}

cat("✅ 所有點位模擬與比對完成！輸出資料夾：", output_folder, "\n")

