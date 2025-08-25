# 載入必要的套件
library(terra)

# --- 設定輸入檔案路徑 ---
input_tif_path <- "C:/Users/user/Downloads/1718_229_198_tree20.tif"
output_path_png <- "output_image.png"

# --- 載入資料 ---
lu_map <- rast(input_tif_path)

# --- 定義 Land Use (LU) 分類顏色配置 ---
categories_df <- data.frame(
  value = c(1, 2, 3, 4, 5, 6),
  label = c("1-Forest", "2-Built-up", "3-Water", "4-Agri", "5-Unkn", "6-Tea")
)

# 設置顏色
custom_colors <- c(
  "darkgreen",  # 1-Forest
  "brown",      # 2-Built-up
  "lightblue",  # 3-Water
  "yellow",     # 4-Agri
  "gray",       # 5-Unkn
  "yellowgreen" # 6-Grass
)

# 設置類別，這是為了讓 `plot` 函式知道如何處理這些數值
levels(lu_map) <- categories_df

# --- 繪製並儲存 PNG 檔案 ---
png(output_path_png, width = 800, height = 600, units = "px", res = 100)
plot(lu_map,
     main = "Land Use Map for 2018",
     # 使用 breaks 參數明確定義每個類別的範圍，並使用 col 參數綁定顏色
     breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
     col = custom_colors,
     legend = TRUE,
     # 修正後的圖例參數
     plg = list(
       title = "Land Use", 
       legend = categories_df$label, 
       cex = 1.0, 
       yjust = 1.05
     ),
     axes = TRUE)
dev.off()

cat(paste0("PNG 地圖已成功儲存至: ", output_path_png, "\n"))
