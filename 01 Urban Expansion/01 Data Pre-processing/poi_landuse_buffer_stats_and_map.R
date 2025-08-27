# convert the address to latitude and longitude

setwd("C:/Users/cwiev/Desktop/intern/july 22/Changhua_Shops")

# 載入套件
library(httr)
library(jsonlite)
library(ggmap)
library(sp)
library(sf) 
library(stars)
library(dplyr)

#地址轉換函數
addr_to_lalo <- function(addr = c("中央研究院"))
{
  # 設定 Google Maps API 金鑰
  api_key <- read.csv("./api.key.txt",header=FALSE)$V1
  
  # 將地址進行 URL 編碼
  encoded_address <- URLencode(addr)
  
  # 建構 API 請求 URL，包含中文語言參數
  url <- paste0("https://maps.googleapis.com/maps/api/geocode/json?address=", 
                encoded_address,"&key=", api_key, "&language=zh-TW")
  
  # 發送 GET 請求
  response <- GET(url)
  
  # 檢查回應狀態
  if (status_code(response) != 200) {
    stop("API 請求失敗")
  }
  
  # 解析 JSON 回傳資料
  data <- content(response, as = "text", encoding = "UTF-8")
  json_data <- fromJSON(data)
  print(json_data)
  # 提取經緯度
  if (json_data$status == "OK") {
    lat <- json_data$results[[3]]$location$lat
    lon <- json_data$results[[3]]$location$lng
    print(paste("Coordinate:",lat,",",lon,sep=""))
    #Convert coordinate reference system    
    #create a sp object with la lo and ref coordinate stsyetm)
    cod.wgs84 <- SpatialPoints( coords = cbind(lon, lat),proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs') )   
    cod.twd97 <- spTransform(cod.wgs84, sp::CRS("+init=epsg:3826"))
    # 顯示結果
    print(paste("Address:",addr," Coordinates(WGS84)/", "Lng:",cod.wgs84@coords[1], "Lat:", cod.wgs84@coords[2], sep=" "))
    # 顯示結果
    print(paste("Address:",addr," Coordinates(TWD97)/", "X:",cod.twd97@coords[1], "Y:", cod.twd97@coords[2], sep=" "))
    
    return(data.frame(lat = lat, lon = lon, x=cod.twd97@coords[1], y= cod.twd97@coords[2]))
  } else {
    stop("地址無法解析")
  }#end if  
}#end of addr2lalo

# 顯示地圖的函數 輸入台灣的地址
addr_map <- function(lalo=FALSE, lat=23.5, lon=121.0,  addr=c("台北101"), zoom = 14, rad = 0.01, year=2023, map_id=2, plot_ld=TRUE)
{
  # INIT
  #  addr=c("台北101")
  #  zoom=14
  #  rad=0.01
  # 設定 Google Maps API 金鑰
  api_key <- read.csv("./api.key.txt",header=FALSE)$V1
  ggmap::register_google(key=api_key)
  #利用Google Sataic Map API 功能 進行地址傳換機緯度的變更
  # 取得中心位置的地圖
  
  if (lalo==FALSE) {
    #conver the address to latitude and longitude 
    cod <- addr_to_lalo(addr=addr)
  }else{
    #load the latitude and longitude 
    cod.wgs84 <- SpatialPoints( coords = cbind(lon, lat),proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs') )
    cod.twd97 <- spTransform(cod.wgs84, sp::CRS("+init=epsg:3826"))
    cod <- data.frame(lon = lon, lat = lat, x=cod.twd97@coords[1], y= cod.twd97@coords[2] )
  }
  map_types=c("roadmap", "satellite", "hybrid", "terrain")
  used_map = map_types[map_id]
  
  #圓框線
  # Convert the location to an sf object
  center_point <- st_point(c(cod$lon, cod$lat)) %>%
    st_sfc(crs = 3857)
  # Define the radius of the circle (in meters)
  radius <- rad # 1000 m
  # Create a circular buffer around the center point
  circle <- st_buffer(center_point, dist = radius)
  circle_sf <- st_transform(circle, crs = 3857) 
  
  #讀取GEOTIFF 檔案並傳換座標
  tiff_file <- read_stars(paste("./lu_maps/",year,"_all_grid_wgs84.tif",sep=""))
  st_crs(tiff_file) <- st_crs(3857)  #TWD97/WGS80
  
  # 使用圓形遮罩裁剪 TIFF 檔案
  tiff_cropped <-st_crop(tiff_file, circle_sf) 
  
  #計算比例
  ## Convert the raster values to a data frame to analyze value frequencies
  raster_values <- as.data.frame(tiff_cropped, as_points = FALSE)
  # Rename the value column for easier access
  names(raster_values)[3] <- "value"
  # Remove NA values
  raster_values <- raster_values %>% filter(!is.na(value)) 
  # Calculate the count and percentage of each unique value
  value_counts <- raster_values %>%
    dplyr::group_by(value) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::mutate(percentage = (count / sum(count)) * 100)
  
  # 給定調色盤
  if (value_counts[1]$value[1] == 0) {
    colors <- c("0"="lightgray","1"="#439c6e","2"="#e86d5f","3"="#8bd2e8","4"="#f0d86e","5"="#999999","6"="#99ad50")
    print(value_counts)
    rad_for = round(value_counts[2,3],digits=1)
    rad_urb = round(value_counts[3,3],digits=1)+round(value_counts[6,3],digits=1)
    rad_wat = round(value_counts[4,3],digits=1)
    rad_agr = round(value_counts[5,3],digits=1)
    rad_gra = round(value_counts[7,3],digits=1)
    print(paste("森林:", rad_for,"%","建物:", rad_urb,"%", "農田:", rad_agr,"%","綠地:", rad_gra,"%","水體:", rad_wat,"%",sep=" "))
  }else{
    colors <- c("1"="#439c6e","2"="#e86d5f","3"="#8bd2e8","4"="#f0d86e","5"="#999999","6"="#99ad50")
    # Print the result
    print(value_counts)
    rad_for = round(value_counts[1,3],digits=1)
    rad_urb = round(value_counts[2,3],digits=1)+round(value_counts[5,3],digits=1)
    rad_wat = round(value_counts[3,3],digits=1)
    rad_agr = round(value_counts[4,3],digits=1)
    rad_gra = round(value_counts[6,3],digits=1)
    print(paste("森林:", rad_for,"%","建物:", rad_urb,"%", "農田:", rad_agr,"%","綠地:", rad_gra,"%","水體:", rad_wat,"%",sep=" "))
  }
  
  # get map from google
  if (plot_ld == T) {
    map <- get_map(location = c(lon = cod$lon, lat = cod$lat), zoom = zoom, source = "google", maptype = used_map)
    base_map <- ggmap(map)
    
    # Plot the circle on the map
    map_with_circle <- base_map +
      geom_stars(data = tiff_cropped, alpha = 0.7) +
      geom_sf(data = circle_sf, fill = NA, color = "gray",  size=1, linewidth=1,  inherit.aes = FALSE) + #顯示圓邊界 
      scale_fill_gradientn(colours = colors[] ) +
      geom_point(aes(x = lon, y = lat), data=cod, color = "yellow", size=1, shape=21, stroke=2) +
      labs(title =  paste("森林(F):", rad_for,"%", ", 農田/綠地(A/G):", rad_agr+rad_gra,"%",", 水體(W):", rad_wat,"%","裸土/建物(B/B):", rad_urb,"%",sep=""),
           caption = paste("查詢地址:",addr,",圖資年份:",year,", TWD97-參考座標(",round(x=cod$x,digits=1),", ", round(x=cod$y,digits=1),")",sep=""),
           x = "經度(WGS84)",  y = "緯度(WGS84)") +
      theme(plot.caption = element_text(size = 10), legend.position ="none") +# Adjust size here & 移除legend 
      coord_sf(datum = st_crs(circle_sf))
    
    # plot the result
    print(map_with_circle)
  } #end if
  
  # return a list of share of land 
  return(list(year=year,forest=rad_for[[1]],agri=rad_agr[[1]], grass=rad_gra[[1]], water=rad_wat[[1]], urban=rad_urb[[1]], radius=rad)) 
}#addr_map

# import the table data for the analysis
in_table <- read.csv("./Changhua_shops/彰化便利店_XY.txt")

rad = c(0.00500)

for ( irun in 1:1) {
  xydel = rad[irun]
  
  lu_data <- data.frame(name=character(), year_0 = numeric(),lon = numeric(), lat =numeric(), x = numeric(), y = numeric(),
                        forest = numeric(), agri= numeric(),year_1=numeric(), year_2=numeric(),status=numeric(),
                        grass = numeric(), water = numeric(), urban = numeric(), radius = numeric() ) 
  
  #loop for year for check the dynamics of land use change
  yrs <- c(2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023)
  
  # sample loop
  for (idf in 1:775 ) {
    for (iyr in yrs)
    {
      plot_ld = FALSE
      addr_value <- addr_map(lalo=TRUE, lat= in_table$lat[idf] , lon=in_table$lon[idf] , zoom = 16,
                             year=iyr, addr=in_table$分公司地址[idf] , plot_ld=plot_ld, rad=xydel)
      addr_value$year_0 <- iyr
      addr_value$year_1 <- in_table$分公司核准設立日期[idf] 
      addr_value$year_1 <- (round(addr_value$year_1/10000)+1911)*10000 + addr_value$year_1%%10000
      addr_value$year_2 <- in_table$分公司最後核准變更日期[idf] 
      addr_value$year_2 <- (round(addr_value$year_2/10000)+1911)*10000 + addr_value$year_2%%10000
      if (is.na(addr_value$year_1) ) addr_value$year_1 <- addr_value$year_2 
      if (is.na(addr_value$year_2) ) addr_value$year_2 <- addr_value$year_1
      
      year_0 <- iyr 
      year_1 <- round(addr_value$year_1/10000)
      year_2 <- round(addr_value$year_2/10000)
      
      print(paste("year_0:", year_0 , "year_1:", year_1, "year_2:", year_2,sep=", "))
      
      if ( (year_0 > year_1) & (year_0 > year_2) ) {
        addr_value$status <- 1.0
      } else if ( year_0 == year_2  ) {
        addr_value$status <- 1.0
      } else {
        addr_value$status <- 0.
      }
      
      addr_value$x <- in_table$x[idf]
      addr_value$y <- in_table$y[idf]
      addr_value$lon <- in_table$lon[idf]
      addr_value$lat <- in_table$lat[idf]
      
      #combine the results 
      lu_data <- rbind(lu_data, cbind(name= in_table$分公司統一編號[idf] ,as.data.frame(addr_value,  stringsAsFactors = FALSE)))
    }# end of iyr
  }# end of idf 
  
  # grop the data & save the group table 
  ld_go <- TRUE 
  if (ld_go) {
    #group the table
    #Group by site_name and count the number of unique statuses for each site
    site_status_counts <- lu_data %>%
      group_by(name) %>%
      summarise(
        unique_status_count = n_distinct(status) # Count distinct statuses
      ) %>%
      ungroup()
    
    print("--- Unique Status Count per Site ---")
    print(site_status_counts)
    cat("\n")
    
    #Filter out sites that have only one unique status
    sites_to_keep <- site_status_counts %>%
      filter(unique_status_count > 1) %>% # Keep sites with more than 1 unique status
      pull(name) # Extract just the site names
    
    print("--- Sites to Keep (those with multiple statuses) ---")
    print(sites_to_keep)
    cat("\n")
    
    #Filter the original data to keep only the identified sites
    filtered_data <- lu_data %>%
      filter(name %in% sites_to_keep)
    
    print("--- Filtered Data (Sites with single status removed) ---")
    print(filtered_data)
    cat("\n")
    
    group_table <- filtered_data %>%
      group_by(name, status) %>%
      summarise(
        urban_median = median(urban, na.rm = TRUE),
        urban_max = max(urban, na.rm = TRUE),
        urban_min = min(urban, na.rm = TRUE),
        lon = mean(lon, na.rm = TRUE),
        lat = mean(lat, na.rm = TRUE),
        x = mean(x, na.rm = TRUE),
        y = mean(y, na.rm = TRUE)
      ) %>%
      ungroup()
    
    #subset the table with status and then combine it together
    sta_0 <- subset(group_table, status==0)
    sta_1 <- subset(group_table, status==1)
    group_table <- left_join(sta_0, sta_1, by="name")
    group_table <- subset(group_table, (urban_median.y - urban_median.x) > 5)  
    
    #plot the results
    boxplot(group_table$urban_min.x, group_table$urban_max.y)
    
    write.table(x=group_table,file=paste("Changhua_shops_group_",xydel,"_df_out_test500m.txt",sep=""),sep=",", row.names=FALSE)
  } # end ld_go
  
  write.table(x=lu_data,file=paste("Changhua_shops_",xydel,"_df_out_test500m.txt",sep=""),sep=",", row.names=FALSE)
  print(paste("Working on grid:",idf,"/Totoal:",length(in_table$x),sep=""))
}# END IRUN
