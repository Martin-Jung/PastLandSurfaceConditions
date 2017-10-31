#### IDEAL Time series ####
library(rvest)
library(zoo)
URL <- "http://e4ftl01.cr.usgs.gov/MOTA/MCD43A4.005/" # Load MCD43 product
pg <- read_html(URL)
# Scrap all dates
dates <- html_attr(html_nodes(pg, "a"), "href")
# Remove the everything up to and including to mota
dates <- dates[seq(1,grep(pattern = "/MOTA/",dates))*-1]
# Remove /
dates <- unlist(lapply(dates, function(x) str_replace(x,"/","") ) )
# This should be sorted now, so create a TS and save
ideal <- zooreg(NA,order.by = ymd(dates))
saveRDS(ideal,"MODIS_MCD43A4_IdealTimeseries.rds")
