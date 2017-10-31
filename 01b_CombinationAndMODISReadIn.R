library(dplyr)
library(lubridate)

# Full 
sites <- readRDS("sites_center_all.rds") 
sites <- sites[which(sites$Max_linear_extent <= quantile(sites$Max_linear_extent,.99,na.rm=T)),] # Get only sites within the 99 quantile

length(unique(sites$SSBS))
sites$midyear <- year(sites$Sample_midpoint)
sites$StudyLength <- (sites$Sample_end_latest - sites$Sample_start_earliest) # Study duration
sites <- subset(sites,StudyLength <= 365)  # Kick out studies longer that a 1 year
length(unique(sites$SSBS))
saveRDS(sites, "AllPredictsSites_full.rds" )


library(rgdal)
sp <- sites[,c("SS","SSBS","Longitude","Latitude")]
coordinates(sp) <- ~Longitude+Latitude
proj4string(sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
writeOGR(sp,"PREDICTS_centerCoordinates.kml","PREDICTS_buffer",driver = "KML",overwrite_layer = T)

#### Load and format for all PREDICTS sites - VERSION 2 - geoJSON ####
myLog <- function(...) {
  cat(paste0("[Extract] ", Sys.time(), " | ", ..., "\n"))
}
library(sp)
library(rgeos)
library("plotKML")
library(raster)
library(rgdal)
library(stringr)
library(data.table)
library(assertthat)
library(jsonlite)
# Load in al MODIS BRDF Bands
readInFormatv2 <- function(x,idv = "SSBS"){
  val =  as.data.frame(fromJSON(x)$features$properties) %>% dplyr::select(-latitude, -longitude) %>% 
    reshape2::melt(., id.vars = idv) %>%
    distinct() %>% mutate(variable = as.character(variable)) %>% 
    mutate(year = as.numeric( str_sub(variable,13,16)), # Get year out of file name
           month = as.numeric( str_sub(variable,18,19) ), # get month
           day = as.numeric( str_sub(variable,21,22) )) %>%  # get day
    # Make a date column
    mutate(date = ymd(paste(year,month,day,sep="-")))
  # Reorder
  val <- val[order(val$date,decreasing = F),]
  return(val)
}

#x = as.data.frame(jsonlite::fromJSON("Extracts/PREDICTS_center_MCD43A4_Band1.geojson")$features$properties)

myLog("Get all BRDF bands files in a formatted way")
b1 <- readInFormatv2("Extracts/PREDICTS_center_MCD43A4_Band1.geojson",idv = "SSBS") %>% mutate(Band = 1)
saveRDS(b1,"Extracts/PREDICTS_center_MCD43A4_1.rds");rm(b1);gc()
b2 <- readInFormatv2("Extracts/PREDICTS_center_MCD43A4_Band2.geojson",idv = "SSBS") %>% mutate(Band = 2)
saveRDS(b2,"Extracts/PREDICTS_center_MCD43A4_2.rds");rm(b2);gc()
b3 <- readInFormatv2("Extracts/PREDICTS_center_MCD43A4_Band3.geojson",idv = "SSBS") %>% mutate(Band = 3)
saveRDS(b3,"Extracts/PREDICTS_center_MCD43A4_3.rds");rm(b3);gc()
b4 <- readInFormatv2("Extracts/PREDICTS_center_MCD43A4_Band4.geojson",idv = "SSBS") %>% mutate(Band = 4)
saveRDS(b4,"Extracts/PREDICTS_center_MCD43A4_4.rds");rm(b4);gc()
b5 <- readInFormatv2("Extracts/PREDICTS_center_MCD43A4_Band5.geojson",idv = "SSBS") %>% mutate(Band = 5)
saveRDS(b5,"Extracts/PREDICTS_center_MCD43A4_5.rds");rm(b5);gc()
b6 <- readInFormatv2("Extracts/PREDICTS_center_MCD43A4_Band6.geojson",idv = "SSBS") %>% mutate(Band = 6)
saveRDS(b6,"Extracts/PREDICTS_center_MCD43A4_6.rds");rm(b6);gc()
b7 <- readInFormatv2("Extracts/PREDICTS_center_MCD43A4_Band7.geojson",idv = "SSBS") %>% mutate(Band = 7)
saveRDS(b7,"Extracts/PREDICTS_center_MCD43A4_7.rds");rm(b7);gc()
myLog("DONE")

