library(terra)
library(sf)
library(dplyr)

## ========= Utils =========
fun_ca_step <- function(in_mat, tg_value = 1, accumulate = FALSE) {
  if (is.null(dim(in_mat))) stop("in_mat must be a matrix")
  nx <- nrow(in_mat); ny <- ncol(in_mat)
  mat <- ifelse(in_mat == tg_value, 1L, 0L); temp <- mat; out_acc <- matrix(0L, nx, ny)
  for (i in 1:nx) for (j in 1:ny) {
    E <- ifelse(i+1 == nx+1, 1, i+1); W <- ifelse(i-1 == 0, nx, i-1)
    N <- ifelse(j-1 == 0, ny, j-1);  S <- ifelse(j+1 == ny+1, 1, j+1)
    nlive <- mat[W,N]+mat[i,N]+mat[E,N]+mat[W,j]+mat[E,j]+mat[W,S]+mat[i,S]+mat[E,S]
    if (mat[i,j]==1 && nlive<2) temp[i,j] <- 0L
    if (mat[i,j]==1 && nlive>3) temp[i,j] <- 0L
    if (mat[i,j]==1 && (nlive==2 || nlive==3)) temp[i,j] <- 1L
    if (mat[i,j]==0 && nlive==3) temp[i,j] <- 1L
  }
  if (accumulate) { out_acc[temp==1L] <- 1L; return(list(out_arr=out_acc)) } else { return(list(out_arr=temp)) }
}

## ========= 路徑 =========
BASE_R_PATH  <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/lu_maps/2023_all_grid_wgs84.tif"
STORES_CSV   <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/csv/all_convenience_stores.csv"
OUT_DIR      <- "C:/Users/cwiev/Desktop/intern/august 13/future_20yrs_CA"
JUNCTIONS_CSV<- file.path(OUT_DIR, "junctions_extracted.csv") # 已先萃取好的路口點 CSV（含 name,lon,lat）

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "growth_diff"), showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "new_growth"),  showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "geojson"),     showWarnings=FALSE)

## ========= 共用參數 =========
RADIUS_M  <- 750
ALIVE_CLASSES <- c(100)     # 100=都市；50=次級；其他=0
GROWTH_MINMAX <- c(1,30)
ZONE_CAPS     <- c(30,20)   # 150m/500m 每年上限的累進空間
STORE_INFLUENCE_RADIUS_M <- 150
STORE_PROB_BASE   <- 0.50
STORE_PROB_TARGET <- 0.65   # 降低一點，避免過度樂觀
STORE_PROB_MULT   <- STORE_PROB_TARGET/STORE_PROB_BASE
EXIST_MIN_DIST_M  <- 400    # 與既有店至少相距（收嚴）
JUNC_MIN_SPACING_M<- 150    # 候選路口先抽稀（網格法）

set.seed(42); terraOptions(progress=0)

## ========= Stage-B0 快篩門檻 (預處理) =========
MIN_AVAILABLE_150_B0 <- 60  # 快篩前先檢查可用空間

## ========= Stage-B 快篩門檻（輕模擬） =========
N_YEARS_B          <- 6
FEEDBACK_THRESH_B  <- 40
MIN_RATE_150_B     <- 0.25
MIN_NEW_150_B      <- 20

## ========= Stage-A 深模擬門檻（嚴格） =========
N_YEARS_A          <- 20
FEEDBACK_THRESH_A  <- 60
MIN_AVAILABLE_150_A<- 80
MIN_RATE_150_A     <- 0.35
MIN_NEW_150_A      <- 40

## ========= Stage-C 決策優化（參數保留，實作將改為 base-R）=========
MAX_NEW_SITES       <- 388   # 最多選幾家（保留參數名以相容）
SITES_MIN_SPACING_M <- 500   # 彼此最小間距

## ========= 載入資料 =========
r_base <- rast(BASE_R_PATH)

existing_sf <- NULL
if (file.exists(STORES_CSV)) {
  s <- read.csv(STORES_CSV)
  if (!("lon" %in% names(s)) && "lon.x" %in% names(s)) s$lon <- s$lon.x
  if (!("lat" %in% names(s)) && "lat.x" %in% names(s)) s$lat <- s$lat.x
  s <- s[!is.na(s$lon)&!is.na(s$lat),]
  if (nrow(s)>0) existing_sf <- st_as_sf(s, coords=c("lon","lat"), crs=4326, remove=FALSE)
}

stopifnot(file.exists(JUNCTIONS_CSV))
poi_df <- read.csv(JUNCTIONS_CSV)
if (!("name" %in% names(poi_df))) poi_df$name <- paste0("jx_", seq_len(nrow(poi_df)))
stopifnot(all(c("lon","lat") %in% names(poi_df)))

# 先建立 poi_sf，再取樣（修正一行順序，避免未定義變數）
poi_sf <- st_as_sf(poi_df, coords=c("lon","lat"), crs=4326, remove=FALSE)

## 排除太靠近既有店
if (!is.null(existing_sf) && nrow(existing_sf)>0) {
  buf <- st_buffer(st_transform(existing_sf, 3826), EXIST_MIN_DIST_M)
  keep <- lengths(st_intersects(st_transform(poi_sf,3826), buf))==0
  poi_sf <- poi_sf[keep,]
  cat(sprintf("排除與既有店 < %dm：剩 %d 點\n", EXIST_MIN_DIST_M, nrow(poi_sf)))
}
## 候選點抽稀
poi_m <- st_transform(poi_sf, 3826)
grid  <- st_make_grid(poi_m, cellsize = JUNC_MIN_SPACING_M)
join  <- st_join(poi_m, st_sf(gid = seq_along(grid), geometry = grid), join = st_within)
keep_idx <- tapply(seq_len(nrow(poi_m)), join$gid, function(v) v[1])
poi_sf <- st_transform(poi_m[unlist(keep_idx), ], 4326)
cat(sprintf("抽稀後候選點：%d\n", nrow(poi_sf)))

# 抽稀後再隨機抽 100（如不足 100 就全取）
set.seed(42)
poi_sf <- poi_sf[sample(seq_len(nrow(poi_sf)), size = min(100, nrow(poi_sf))), ]
cat(sprintf("抽稀後取樣 100：%d 點\n", nrow(poi_sf)))

poi_data <- poi_sf |> st_drop_geometry() |> dplyr::select(name,lon,lat) |>
  filter(!is.na(lon), !is.na(lat))
poi_data$lon <- as.numeric(poi_data$lon); poi_data$lat <- as.numeric(poi_data$lat)
stopifnot(nrow(poi_data)>0)

## ========= 共用模擬函數 =========
run_sim <- function(shop_lon, shop_lat, N_YEARS, FEEDBACK_THRESH,
                    save_outputs=FALSE, sim_idx=NULL, name_for_title="") {
  delta_deg <- RADIUS_M/111000
  ext_crop <- ext(shop_lon - delta_deg, shop_lon + delta_deg,
                  shop_lat - delta_deg, shop_lat + delta_deg)
  r0 <- try(crop(r_base, ext_crop), silent=TRUE)
  if (inherits(r0,"try-error") || ncell(r0)==0) return(NULL)
  vals0 <- values(r0); vals0[vals0==2] <- 100; vals0[vals0==5] <- 50; vals0[!(vals0 %in% c(100,50))] <- 0
  nX <- nrow(r0); nY <- ncol(r0)
  mat0 <- matrix(vals0, nrow=nX, ncol=nY, byrow=TRUE)
  bin0 <- matrix(as.integer(mat0 %in% ALIVE_CLASSES), nrow=nX, ncol=nY, byrow=TRUE)
  coords <- crds(r0)
  d_m  <- sqrt((coords[,1]-shop_lon)^2 + (coords[,2]-shop_lat)^2)*111000
  zone <- ifelse(d_m<=150,1, ifelse(d_m<=500,2,3))
  p_inv<- 1/(d_m+1); p_inv <- p_inv/max(p_inv)
  within <- d_m <= STORE_INFLUENCE_RADIUS_M
  p_use <- p_inv; p_use[within] <- pmin(1, p_inv[within]*STORE_PROB_MULT)
  current <- vals0; cap <- vals0
  cap[zone==1] <- pmin(100, vals0[zone==1] + ZONE_CAPS[1])
  cap[zone==2] <- pmin(100, vals0[zone==2] + ZONE_CAPS[2])
  curr_bin <- bin0
  for (yr in 1:N_YEARS) {
    ca_one   <- fun_ca_step(curr_bin, tg_value=1, accumulate=FALSE)
    cand_vec <- as.vector(t(ca_one$out_arr))
    cap_t <- cap; cap_t[cand_vec==0 | zone==3] <- current[cand_vec==0 | zone==3]
    pickable <- (cand_vec==1) & (zone!=3)
    selected <- rep(FALSE, length(current))
    if (any(pickable)) selected[pickable] <- runif(sum(pickable)) < p_use[pickable]
    growth <- rep(0, length(current))
    if (any(selected)) growth[selected] <- runif(sum(selected), GROWTH_MINMAX[1], GROWTH_MINMAX[2])
    proposed <- current + growth
    current  <- ifelse(proposed<=cap_t, proposed, cap_t)
    curr_bin <- matrix(as.integer(current > FEEDBACK_THRESH), nrow=nX, ncol=nY, byrow=TRUE)
  }
  future_bin <- matrix(as.integer(current > FEEDBACK_THRESH), nrow=nX, ncol=nY, byrow=TRUE)
  new_bin    <- ifelse(bin0==0 & future_bin==1, 1, NA)
  new_vec    <- as.vector(t(new_bin))
  bin0_vec   <- as.vector(t(bin0))
  new_total  <- sum(new_vec==1, na.rm=TRUE)
  avail150   <- sum(bin0_vec==0 & zone==1, na.rm=TRUE)
  new150     <- sum(new_vec==1  & zone==1, na.rm=TRUE)
  rate150    <- ifelse(avail150==0, 0, new150/avail150)
  avail500   <- sum(bin0_vec==0 & zone==2, na.rm=TRUE)
  new500     <- sum(new_vec==1  & zone==2, na.rm=TRUE)
  rate500    <- ifelse(avail500==0, 0, new500/avail500)
  if (isTRUE(save_outputs)) {
    diff_vec <- current - vals0
    r_diff <- r0; values(r_diff) <- diff_vec
    png(file.path(OUT_DIR, "growth_diff", sprintf("sim_plot_%03d_growth_diff.png", sim_idx)), 800, 600)
    plot(r_diff, col=colorRampPalette(c("white","orange","red"))(100),
         main=sprintf("Growth (%dy) - %s | 150m %d/%d (%.0f%%)", N_YEARS, name_for_title, new150, avail150, 100*rate150))
    points(shop_lon, shop_lat, col="darkblue", pch=19)
    symbols(shop_lon, shop_lat, circles=150/111000, add=TRUE, inches=FALSE, fg="blue", lty=2)
    symbols(shop_lon, shop_lat, circles=500/111000, add=TRUE, inches=FALSE, fg="red",  lty=2)
    dev.off()
    r_new <- r0; values(r_new) <- as.vector(t(new_bin))
    png(file.path(OUT_DIR, "new_growth", sprintf("sim_plot_%03d_new_growth.png", sim_idx)), 800, 600)
    plot(r_new, col="red", legend=FALSE,
         main=sprintf("New Urban (0→1) - %s | 150m %d/%d (%.0f%%)", name_for_title, new150, avail150, 100*rate150))
    points(shop_lon, shop_lat, col="darkblue", pch=19)
    symbols(shop_lon, shop_lat, circles=150/111000, add=TRUE, inches=FALSE, fg="blue", lty=2)
    symbols(shop_lon, shop_lat, circles=500/111000, add=TRUE, inches=FALSE, fg="red",  lty=2)
    dev.off()
    r_future <- r0; values(r_future) <- current
    writeRaster(r_future, file.path(OUT_DIR, sprintf("future_continuous_%03d.tif", sim_idx)), overwrite=TRUE)
  }
  list(avail150=avail150, new150=new150, rate150=rate150,
       avail500=avail500, new500=new500, rate500=rate500,
       new_total=new_total)
}


## == Stage-B0 可用空間預篩選 ==
cat("\n== Stage-B0 可用空間預篩選開始 ==\n")
idx_keep <- logical(nrow(poi_data))
for (i in seq_len(nrow(poi_data))) {
  shop_lon <- poi_data$lon[i]; shop_lat <- poi_data$lat[i]
  delta_deg <- RADIUS_M / 111000
  ext_crop <- ext(shop_lon - delta_deg, shop_lon + delta_deg,
                  shop_lat - delta_deg, shop_lat + delta_deg)
  r0 <- try(crop(r_base, ext_crop), silent = TRUE)
  if (inherits(r0,"try-error") || ncell(r0)==0) { idx_keep[i] <- FALSE; next }
  
  nX <- nrow(r0); nY <- ncol(r0)
  if (is.null(nX) || is.null(nY) || nX==0 || nY==0) { idx_keep[i] <- FALSE; next }
  
  vals0 <- values(r0)
  vals0[vals0==2] <- 100
  vals0[vals0==5] <- 50
  vals0[!(vals0 %in% c(100,50))] <- 0
  
  # 與主程式一致：先二值化（100 為都市=1，其餘=0）
  bin0 <- matrix(as.integer(matrix(vals0, nrow=nX, ncol=nY, byrow=TRUE) %in% ALIVE_CLASSES),
                 nrow=nX, ncol=nY, byrow=TRUE)
  
  coords <- crds(r0)
  d_m <- sqrt((coords[,1]-shop_lon)^2 + (coords[,2]-shop_lat)^2) * 111000
  zone <- ifelse(d_m<=150,1, ifelse(d_m<=500,2,3))
  
  bin0_vec <- as.vector(t(bin0))
  available_150 <- sum(bin0_vec==0 & zone==1, na.rm=TRUE)
  
  idx_keep[i] <- available_150 >= MIN_AVAILABLE_150_B0
  if (i %% 500 == 0) cat(sprintf("...B0 預篩到第 %d / %d 點\n", i, nrow(poi_data)))
}
poi_data <- poi_data[idx_keep, ]
cat(sprintf("預篩選後剩餘點數：%d\n", nrow(poi_data)))



## ========= Stage-B 快篩開始 =========
cat("\n== Stage-B 快篩開始 ==\n")
recB <- data.frame(
  name = character(), lon = numeric(), lat = numeric(),
  available_150 = numeric(), new_in150 = numeric(), new_rate_150 = numeric(),
  available_500 = numeric(), new_in500 = numeric(), new_rate_500 = numeric(),
  new_total = numeric(),
  stringsAsFactors = FALSE
)
if (nrow(poi_data) > 0) {
  for (i in seq_len(nrow(poi_data))) {
    nm <- poi_data$name[i]; lon <- poi_data$lon[i]; lat <- poi_data$lat[i]
    res <- run_sim(lon, lat, N_YEARS_B, FEEDBACK_THRESH_B, save_outputs=FALSE)
    if (is.null(res)) next
    passB <- (res$avail150 >= MIN_AVAILABLE_150_B) &&
      (res$rate150  >= MIN_RATE_150_B)     &&
      (res$new150   >= MIN_NEW_150_B)
    if (passB) {
      new_row <- data.frame(
        name = nm, lon = lon, lat = lat,
        available_150 = res$avail150, new_in150 = res$new150, new_rate_150 = res$rate150,
        available_500 = res$avail500, new_in500 = res$new500, new_rate_500 = res$rate500,
        new_total = res$new_total,
        stringsAsFactors = FALSE
      )
      recB <- rbind(recB, new_row)
    }
    if (i %% 100 == 0) cat(sprintf("...B 篩到第 %d / %d 點；目前保留 %d\n", i, nrow(poi_data), nrow(recB)))
  }
}
if (nrow(recB)>0) {
  write.csv(recB, file.path(OUT_DIR, "stageB_candidates.csv"), row.names=FALSE)
}
cat(sprintf("Stage-B 通過：%d / %d\n", nrow(recB), nrow(poi_data)))


## == Stage-A 深模擬開始 ==
cat("\n== Stage-A 深模擬開始 ==\n")
recA <- data.frame(
  name=character(), lon=numeric(), lat=numeric(),
  available_150=numeric(), new_in150=numeric(), new_rate_150=numeric(),
  available_500=numeric(), new_in500=numeric(), new_rate_500=numeric(),
  new_total=numeric(),
  stringsAsFactors = FALSE
)

sim_idx <- 1
if (nrow(recB)>0) {
  for (i in seq_len(nrow(recB))) {
    nm  <- recB$name[i]; lon <- as.numeric(recB$lon[i]); lat <- as.numeric(recB$lat[i])
    res <- run_sim(lon, lat, N_YEARS_A, FEEDBACK_THRESH_A,
                   save_outputs=TRUE, sim_idx=sim_idx, name_for_title=nm)
    if (is.null(res)) next
    
    passA <- (res$avail150 >= MIN_AVAILABLE_150_A) &&
      (res$rate150  >= MIN_RATE_150_A)      &&
      (res$new150   >= MIN_NEW_150_A)
    
    if (passA) {
      new_row <- data.frame(
        name = nm, lon = lon, lat = lat,
        available_150 = as.numeric(res$avail150),
        new_in150     = as.numeric(res$new150),
        new_rate_150  = as.numeric(res$rate150),
        available_500 = as.numeric(res$avail500),
        new_in500     = as.numeric(res$new500),
        new_rate_500  = as.numeric(res$rate500),
        new_total     = as.numeric(res$new_total),
        stringsAsFactors = FALSE
      )
      recA <- rbind(recA, new_row)
      sim_idx <- sim_idx + 1
    }
    
    if (i %% 100 == 0) cat(sprintf("...A 深模擬到第 %d / %d；目前保留 %d\n", i, nrow(recB), nrow(recA)))
    gc()
  }
}

if (nrow(recA) > 0) {
  recA <- recA |>
    arrange(desc(new_rate_150), desc(new_in150), desc(new_rate_500), desc(new_total))
  write.csv(recA, file.path(OUT_DIR, "stageA_strong.csv"), row.names=FALSE)
}
cat(sprintf("Stage-A 通過：%d\n", nrow(recA)))


## =======================
## Stage-C：base-R 版（無 sf/terra）
## 先全掃距離（不中斷），再依表現排序取前 388
## =======================
options(stringsAsFactors = FALSE)

STAGEA_CSV          <- file.path(OUT_DIR, "stageA_strong.csv")
SITES_MIN_SPACING_M <- SITES_MIN_SPACING_M  # 已於前段定義為 500
TOP_N_FINAL         <- MAX_NEW_SITES        # 388

stopifnot(file.exists(STAGEA_CSV))
dir.create(file.path(OUT_DIR, "geojson"), recursive = TRUE, showWarnings = FALSE)

recA2 <- read.csv(STAGEA_CSV, check.names = FALSE)
cat(sprintf("已載入 Stage-A 篩選結果，共 %d 個點位。\n", nrow(recA2)))
if (nrow(recA2) == 0) {
  cat("Stage-A 無點位，結束。\n")
} else {
  # 欄位清理
  stopifnot(all(c("lon","lat") %in% names(recA2)))
  recA2$lon <- as.numeric(recA2$lon)
  recA2$lat <- as.numeric(recA2$lat)
  recA2 <- recA2[!is.na(recA2$lon) & !is.na(recA2$lat), , drop = FALSE]
  
  # 定義排序鍵
  score_cols <- c("new_rate_150","new_in150","new_rate_500","new_total")
  has_cols <- score_cols %in% names(recA2)
  if (all(has_cols)) {
    base_order <- order(-recA2$new_rate_150, -recA2$new_in150, -recA2$new_rate_500, -recA2$new_total,
                        na.last = TRUE, decreasing = FALSE)
  } else {
    base_order <- seq_len(nrow(recA2))
  }
  
  # Haversine 距離（公尺）
  .rad <- pi / 180
  haversine_m <- function(lat1, lon1, lat2, lon2) {
    dlat <- (lat2 - lat1) * .rad
    dlon <- (lon2 - lon1) * .rad
    a <- sin(dlat/2)^2 + cos(lat1*.rad) * cos(lat2*.rad) * sin(dlon/2)^2
    2 * 6371000 * asin(pmin(1, sqrt(a)))
  }
  
  # 全掃描 spacing
  n <- nrow(recA2)
  chosen_idx <- logical(n)
  blocked_by_spacing <- logical(n)
  sel_lon <- numeric(0); sel_lat <- numeric(0)
  
  cat("\n== Stage-C 全掃描（先檢查間距、不中斷）==\n")
  for (k in seq_along(base_order)) {
    i <- base_order[k]
    lon <- recA2$lon[i]; lat <- recA2$lat[i]
    
    if (!length(sel_lon)) {
      chosen_idx[i] <- TRUE
      sel_lon <- c(sel_lon, lon); sel_lat <- c(sel_lat, lat)
    } else {
      d <- haversine_m(lat, lon, sel_lat, sel_lon)
      if (min(d) >= SITES_MIN_SPACING_M) {
        chosen_idx[i] <- TRUE
        sel_lon <- c(sel_lon, lon); sel_lat <- c(sel_lat, lat)
      } else {
        blocked_by_spacing[i] <- TRUE
      }
    }
    if (k %% 200 == 0) cat(sprintf("...掃描 %d/%d，目前可行候選：%d\n", k, n, sum(chosen_idx)))
  }
  
  candidates <- recA2[chosen_idx, , drop = FALSE]
  cat(sprintf("\n通過 %dm 間距的候選總數：%d\n", SITES_MIN_SPACING_M, nrow(candidates)))
  
  # 在候選集合上再排序，最後取前 388
  if (nrow(candidates) > 0) {
    if (all(score_cols %in% names(candidates))) {
      ord <- order(-candidates$new_rate_150, -candidates$new_in150,
                   -candidates$new_rate_500, -candidates$new_total,
                   na.last = TRUE)
    } else {
      ord <- order(match(rownames(candidates), rownames(recA2)))
    }
    candidates_sorted <- candidates[ord, , drop = FALSE]
    selected <- head(candidates_sorted, TOP_N_FINAL)
  } else {
    selected <- candidates
  }
  
  # 落盤：診斷、候選、最終
  diag <- data.frame(
    idx = seq_len(nrow(recA2)),
    name = if ("name" %in% names(recA2)) recA2$name else NA_character_,
    lon = recA2$lon, lat = recA2$lat,
    chosen_by_spacing = chosen_idx,
    blocked_by_spacing = blocked_by_spacing
  )
  write.csv(diag,       file.path(OUT_DIR, "stageC_scan_diag.csv"), row.names = FALSE)
  write.csv(candidates, file.path(OUT_DIR, "stageC_candidates_spacing.csv"), row.names = FALSE)
  write.csv(selected,   file.path(OUT_DIR, "selected_sites.csv"), row.names = FALSE)
  
  # 輸出 GeoJSON（純 base R）
  geojson_path <- file.path(OUT_DIR, "geojson", "selected_sites.geojson")
  json_escape <- function(x) { x <- gsub("\\\\","\\\\\\\\", x); gsub("\"","\\\\\"", x) }
  fmt_val <- function(v) {
    if (is.na(v)) return("null")
    if (is.numeric(v)) return(sub(",", "", format(v, scientific = FALSE, trim = TRUE)))
    paste0("\"", json_escape(as.character(v)), "\"")
  }
  prop_keys <- setdiff(names(selected), c("lon","lat"))
  feature_strings <- character(nrow(selected))
  for (i in seq_len(nrow(selected))) {
    props <- paste0(
      mapply(function(k, val) paste0("\"", json_escape(k), "\":", fmt_val(val)),
             prop_keys, selected[i, prop_keys, drop = TRUE]),
      collapse = ","
    )
    geom <- paste0("\"type\":\"Point\",\"coordinates\":[",
                   fmt_val(selected$lon[i]), ",", fmt_val(selected$lat[i]), "]")
    feature_strings[i] <- paste0(
      "{\"type\":\"Feature\",\"geometry\":{", geom, "},\"properties\":{", props, "}}"
    )
  }
  geojson_text <- paste0(
    "{\"type\":\"FeatureCollection\",\"features\":[",
    paste(feature_strings, collapse = ","),
    "]}"
  )
  writeLines(geojson_text, geojson_path)
  
  cat(sprintf("\n✅ 最終選址：%d 個（間距 ≥ %dm；候選 %d → 取前 %d）\n",
              nrow(selected), SITES_MIN_SPACING_M, nrow(candidates), TOP_N_FINAL))
  cat("輸出：stageC_scan_diag.csv、stageC_candidates_spacing.csv、selected_sites.csv + GeoJSON\n")
}
