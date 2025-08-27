library(sf)
library(dplyr)

# ================== 參數設定 (與主程式相同) ==================
ROADS_SHP_PATH  <- "C:\\Users\\cwiev\\Desktop\\intern\\august 13\\changhua_roads\\changhua_roads.shp"
OUT_DIR         <- "C:/Users/cwiev/Desktop/intern/august 13/future_20yrs_CA"
MIN_DIST_FROM_EXISTING_M <- 80
STORES_CSV_PATH <- "C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops/csv/all_convenience_stores.csv"
SCENARIO <- "junctions_only_new" # or "co_locate"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ================== 萃取與輸出路口點 ==================
cat("正在從路網圖中自動萃取交叉路口...\n")
roads <- suppressWarnings(st_read(ROADS_SHP_PATH, quiet = TRUE))
if (is.na(st_crs(roads))) st_crs(roads) <- 4326
roads_m <- st_transform(roads, 3826)

inter_union <- st_union(roads_m)
inter_geom  <- st_intersection(inter_union)
inter_pts   <- st_cast(inter_geom, "POINT")
if (length(inter_pts) == 0 || all(st_is_empty(inter_pts))) {
  stop("無法從路網取得交叉點，請檢查 shp 或其幾何類型。")
}

cat("正在對交叉點進行抽稀...\n")
grid <- st_make_grid(inter_pts, cellsize = 30)
join <- st_join(
  st_as_sf(data.frame(id = seq_along(inter_pts)), geometry = inter_pts),
  st_sf(gid = seq_along(grid), geometry = grid), join = st_within
)
agg <- suppressWarnings(st_centroid(st_geometry(grid[na.omit(join$gid), ])))
poi_sf <- st_as_sf(data.frame(name = paste0("jx_", seq_along(agg))), geometry = agg, crs = st_crs(inter_pts))
poi_sf <- st_transform(poi_sf, 4326)

# 排除太靠近既有店的候選點 (與主程式相同，確保一致性)
existing_sf <- NULL
if (file.exists(STORES_CSV_PATH)) {
  stores_df <- read.csv(STORES_CSV_PATH)
  if (!("lon" %in% names(stores_df)) && "lon.x" %in% names(stores_df)) stores_df$lon <- stores_df$lon.x
  if (!("lat" %in% names(stores_df)) && "lat.x" %in% names(stores_df)) stores_df$lat <- stores_df$lat.x
  stores_df <- stores_df[!is.na(stores_df$lon) & !is.na(stores_df$lat), ]
  if (nrow(stores_df) > 0) {
    existing_sf <- st_as_sf(stores_df, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  }
}
if (!is.null(existing_sf) && nrow(existing_sf) > 0 && SCENARIO == "junctions_only_new") {
  buf <- st_buffer(st_transform(existing_sf, 3826), MIN_DIST_FROM_EXISTING_M)
  poi_keep <- lengths(st_intersects(st_transform(poi_sf, 3826), buf)) == 0
  poi_sf <- poi_sf[poi_keep, ]
  cat(sprintf("已篩選掉與現有超商距離 < %dm 的路口點，剩餘 %d 個。\n", MIN_DIST_FROM_EXISTING_M, nrow(poi_sf)))
}

poi_data <- cbind(st_drop_geometry(st_transform(poi_sf, 4326)), st_coordinates(st_transform(poi_sf, 4326)))
names(poi_data)[(ncol(poi_data)-1):ncol(poi_data)] <- c("lon","lat")
if (!"name" %in% names(poi_data)) poi_data$name <- paste0("jx_", seq_len(nrow(poi_data)))
if (nrow(poi_data) == 0) stop("沒有路口點可供模擬，請檢查 shp 或篩選條件。")

JUNCTIONS_PATH <- file.path(OUT_DIR, "junctions_extracted.csv")
write.csv(poi_data, JUNCTIONS_PATH, row.names = FALSE)
cat("\n🎉 路口點已成功輸出到：", JUNCTIONS_PATH, "\n")

