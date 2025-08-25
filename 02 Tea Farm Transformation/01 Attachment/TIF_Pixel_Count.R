#install.packages("raster")
library(raster)

# 載入 raster
# filename <- basename("C:/Users/user/Downloads/2024_229_198_tree20.tif")
# filepath <- file.path("C:/Users/user/Downloads", filename)

filepath <- "C:/Users/user/Downloads/1718_229_198_tree20.tif"
aa <- raster(filepath)

# 取出像素值並去除 NA
vals <- getValues(aa)
vals <- vals[!is.na(vals)]
print(length(vals))

# 計算每個像素值的出現次數 (count)
counts <- table(vals)
df_counts <- as.data.frame(counts)
colnames(df_counts) <- c("value", "count")
print(counts)

# 加入百分比欄位
df_counts$percentage <- round(100 * df_counts$count / length(vals), 2)  # 四捨五入到小數第 2 位

# 圖標題與 Y 軸上限
year <- substr(filename, 1, 4)
# title_text <- paste0("Pixel Value Counts for  ", year," (Grass to Forest)")
title_text <- paste0("Pixel Value Counts for 1718")

max_count <- max(df_counts$count)

# 繪製柱狀圖
bp <- barplot(
  height = df_counts$count,
  names.arg = df_counts$value,
  main = title_text,
  xlab = "Pixel Value",
  ylab = "Count",
  ylim = c(0, max_count * 1.2)  # 稍微多一點空間放上標籤
)

# 顯示數量 + 百分比
text(
  x = bp,
  y = df_counts$count * 1.02,
  #labels = df_counts$count,
  # labels = df_counts$count, "\n(", df_counts$percentage, "%)"),
  labels = paste0(df_counts$percentage, "%"),
  cex = 0.6,
  col = "blue",
  pos = 3
)
