# =======================================================
# Predict 2023 from 2013 for ALL points
# Only save growth-diff PNGs
# =======================================================
library(terra)
library(dplyr)
library(readr)

# ---------- 路徑 ----------
tiff_2013_path <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/lu_maps/2013_all_grid_wgs84.tif"
shops_150_path <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/csv/Changhua_shops_group_0.0015_df_out.csv"
shops_500_path <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/csv/Changhua_shops_group_0.005_df_out.csv"

# 輸出資料夾（只留 growth_diff）
out_root <- "C:/Users/cwiev/Desktop/intern/august 13/ca_prediction_2013to2023"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir_growth_png  <- file.path(out_root, "growth_difference_png")
dir.create(dir_growth_png,  showWarnings = FALSE, recursive = TRUE)

# ---------- 參數 ----------
n_years       <- 10
radius_m      <- 750
ALIVE_CLASSES <- c(100)  # 若 50 也算都市：c(100, 50)
set.seed(42); terraOptions(progress = 0)

# ---------- CA：單步，不累積 ----------
fun_ca_process <- function(in_arr, it_step = 1, tg_value = 1, accumulate = FALSE) {
  if (is.null(in_arr) || is.null(dim(in_arr))) stop("in_arr invalid")
  nx <- nrow(in_arr); ny <- ncol(in_arr)
  out_arr <- matrix(0L, nx, ny); tmp_arr <- out_arr
  mat <- in_arr; mat[in_arr == tg_value] <- 1L; mat[in_arr != tg_value] <- 0L
  temp_mat <- mat
  for (k in 1:it_step) {
    for (i in 1:nx) for (j in 1:ny) {
      E <- ifelse(i+1 == nx+1, 1, i+1)
      W <- ifelse(i-1 == 0,     nx, i-1)
      N <- ifelse(j-1 == 0,     ny, j-1)
      S <- ifelse(j+1 == ny+1,  1,  j+1)
      numb_alive <- mat[W,N] + mat[i,N] + mat[E,N] +
        mat[W,j]            + mat[E,j] +
        mat[W,S] + mat[i,S] + mat[E,S]
      if (mat[i,j]==1 && numb_alive<2)                     temp_mat[i,j] <- 0L
      if (mat[i,j]==1 && numb_alive>3)                     temp_mat[i,j] <- 0L
      if (mat[i,j]==1 && (numb_alive==2 || numb_alive==3)) temp_mat[i,j] <- 1L
      if (mat[i,j]==0 && numb_alive==3)                    temp_mat[i,j] <- 1L
    }
    mat <- temp_mat
    if (accumulate) for (ii in 1:nx) for (jj in 1:ny) if (mat[ii,jj] == 1L) tmp_arr[ii,jj] <- 1L
  }
  if (accumulate) list(out_arr = tmp_arr) else list(out_arr = mat)
}

# ---------- 讀資料 ----------
r2013_all <- rast(tiff_2013_path)
shops_150  <- read.csv(shops_150_path, stringsAsFactors = FALSE)
shops_500  <- read.csv(shops_500_path, stringsAsFactors = FALSE)

# 合併並去重（以 name 唯一）
shop_data <- bind_rows(shops_150, shops_500) %>%
  distinct(name, .keep_all = TRUE)

# 經緯度欄位相容
has_lonx <- "lon.x" %in% names(shop_data); has_latx <- "lat.x" %in% names(shop_data)
has_lon  <- "lon"   %in% names(shop_data); has_lat  <- "lat"   %in% names(shop_data)
stopifnot((has_lonx && has_latx) || (has_lon && has_lat))
get_lon <- function(df, idx) if (has_lonx) df$lon.x[idx] else df$lon[idx]
get_lat <- function(df, idx) if (has_latx) df$lat.x[idx] else df$lat[idx]

# ---------- 全部點位迴圈 ----------
for (i in 1:nrow(shop_data)) {
  name_i  <- shop_data$name[i]
  shop_lon <- suppressWarnings(as.numeric(get_lon(shop_data, i)))
  shop_lat <- suppressWarnings(as.numeric(get_lat(shop_data, i)))
  if (is.na(shop_lon) || is.na(shop_lat)) {
    cat("⚠️ NA coords, skip. name=", name_i, "\n"); next
  }
  
  # 裁切 750m 視窗
  delta_deg <- radius_m / 111000
  ext_crop <- ext(shop_lon - delta_deg, shop_lon + delta_deg,
                  shop_lat - delta_deg, shop_lat + delta_deg)
  r2013 <- try(crop(r2013_all, ext_crop), silent = TRUE)
  if (inherits(r2013, "try-error") || ncell(r2013)==0) {
    cat("⚠️ Crop failed/empty, skip. name=", name_i, "\n"); next
  }
  cat("👉 Running:", name_i, "  (", nrow(r2013), "x", ncol(r2013), ")\n")
  
  # ---- 2013 重編碼 ----
  ca0 <- values(r2013)
  ca0[ca0 == 2] <- 100
  ca0[ca0 == 5] <- 50
  ca0[!(ca0 %in% c(100, 50))] <- 0
  nX <- nrow(r2013); nY <- ncol(r2013)
  ca0_mat <- matrix(ca0, nrow=nX, ncol=nY, byrow=TRUE)
  curr_bin <- matrix(as.integer(ca0_mat %in% ALIVE_CLASSES), nrow=nX, ncol=nY, byrow=TRUE)
  
  # ---- 距離分區 ----
  coords   <- crds(r2013)
  dists_m  <- sqrt((coords[,1] - shop_lon)^2 + (coords[,2] - shop_lat)^2) * 111000
  zone_cls <- ifelse(dists_m <= 150, 1, ifelse(dists_m <= 500, 2, 3))
  inv_score <- 1/(dists_m + 1); inv_score <- inv_score / max(inv_score)
  
  # ---- 連續層初值 & 上限 ----
  current    <- ca0
  max_growth <- ca0
  max_growth[zone_cls == 1] <- pmin(100, ca0[zone_cls == 1] + 30)
  max_growth[zone_cls == 2] <- pmin(100, ca0[zone_cls == 2] + 20)
  
  # ---- 每年演化 ----
  for (year in 1:n_years) {
    cand_bin <- fun_ca_process(curr_bin, it_step=1, tg_value=1, accumulate=FALSE)$out_arr
    cand_vec <- as.vector(t(cand_bin))
    
    max_t <- max_growth
    block <- (cand_vec == 0) | (zone_cls == 3)
    max_t[block] <- current[block]
    
    pickable <- (cand_vec == 1) & (zone_cls != 3)
    selected <- rep(FALSE, length(current))
    if (any(pickable)) selected[pickable] <- runif(sum(pickable)) < inv_score[pickable]
    
    growth <- rep(0, length(current))
    if (any(selected)) growth[selected] <- runif(sum(selected), 1, 100)
    
    proposed <- current + growth
    current  <- ifelse(proposed <= max_t, proposed, max_t)
    curr_bin <- matrix(as.integer(current > 0), nrow=nX, ncol=nY, byrow=TRUE)
  }
  
  # ---- 圖：成長差異 ----
  growth_diff <- current - ca0
  growth_raster <- r2013; values(growth_raster) <- growth_diff
  png(file.path(dir_growth_png, sprintf("growth_diff_%s.png", as.character(name_i))),
      width = 900, height = 700)
  plot(growth_raster, col=colorRampPalette(c("white","orange","red"))(100),
       main = paste0("Urban Growth Difference (2013→2023) - ", name_i))
  points(shop_lon, shop_lat, col="darkblue", pch=19)
  symbols(shop_lon, shop_lat, circles=150/111000, add=TRUE, inches=FALSE, fg="blue", lty=2)
  symbols(shop_lon, shop_lat, circles=500/111000, add=TRUE, inches=FALSE, fg="red",  lty=2)
  mtext(paste0("Lon: ", round(shop_lon,5), "  Lat: ", round(shop_lat,5)), side=1, line=3, cex=0.9)
  dev.off()
}

cat("🎉 Done. Outputs in:\n",
    "- Growth difference PNGs: ", dir_growth_png, "\n", sep = "")
