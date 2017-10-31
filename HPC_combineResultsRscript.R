#### Aggregate results ####
library(lubridate)
library(zoo)
library(stringr)
# Analysis packages
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
#library(greenbrown)
library(zoo)
library(xts)
library(forecast) # additional smoothers
library(signal) # sg-filter
library(imputeTS) # kalman smooth auto function
myLog <- function(...) {
  cat(paste0("[TS_Prepare] ", Sys.time(), " | ", ..., "\n"))
}
co.var <- function(x) ( sd(x, na.rm = T)/mean(x, na.rm = T))
outlier = TRUE
# Ideal time series
ideal <- readRDS("MODIS_MCD43A4_IdealTimeseries.rds")

r1 <- readRDS("Center_FullTimeSeriesSmooth_Final_2001.rds")
r2 <- readRDS("Center_FullTimeSeriesSmooth_Final_2002.rds")
r3 <- readRDS("Center_FullTimeSeriesSmooth_Final_2003.rds")
r4 <- readRDS("Center_FullTimeSeriesSmooth_Final_2004.rds")
r5 <- readRDS("Center_FullTimeSeriesSmooth_Final_2005.rds")
r6 <- readRDS("Center_FullTimeSeriesSmooth_Final_2006.rds")
r7 <- readRDS("Center_FullTimeSeriesSmooth_Final_2007.rds")
r8 <- readRDS("Center_FullTimeSeriesSmooth_Final_2008.rds")
r9 <- readRDS("Center_FullTimeSeriesSmooth_Final_2009.rds")
r10 <- readRDS("Center_FullTimeSeriesSmooth_Final_2010.rds")
r11 <- readRDS("Center_FullTimeSeriesSmooth_Final_2011.rds")
r12 <- readRDS("Center_FullTimeSeriesSmooth_Final_2012.rds")
r13 <- readRDS("Center_FullTimeSeriesSmooth_Final_2013.rds")


# Combine all
r.full <- c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)
rm(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)
saveRDS(r.full,"Center_FullTimeSeriesSmooth_Final.rds")


sites <- readRDS("AllPredictsSites_full.rds")
# Which ones are still missing ?
(missing = sites$SSBS[which(!(sites$SSBS %in% names(r.full)))])
sites <- subset(sites,SSBS %in% missing)

# Reclassifying and running again
md <- data.frame(missing,newid=NA)
md$newid <- str_replace_all(md$missing,pattern = "ß","?")
md$newid <- str_replace_all(md$newid,pattern = "ü","?")
md$newid <- str_replace_all(md$newid,pattern = "é","?")

gc()
myLog("Loading Files")
b1 <- (readRDS("Extracts/PREDICTS_center_MCD43A4_1.rds"))
b1 <- subset(b1,SSBS %in% md$newid)
b2 <- (readRDS("Extracts/PREDICTS_center_MCD43A4_2.rds"))
b2 <- subset(b2,SSBS %in% md$newid)
b3 <- (readRDS("Extracts/PREDICTS_center_MCD43A4_3.rds"))
b3 <- subset(b3,SSBS %in% md$newid)
b4 <- (readRDS("Extracts/PREDICTS_center_MCD43A4_4.rds"))
b4 <- subset(b4,SSBS %in% md$newid)
b5 <- (readRDS("Extracts/PREDICTS_center_MCD43A4_5.rds"))
b5 <- subset(b5,SSBS %in% md$newid)
b6 <- (readRDS("Extracts/PREDICTS_center_MCD43A4_6.rds"))
b6 <- subset(b6,SSBS %in% md$newid)
b7 <- (readRDS("Extracts/PREDICTS_center_MCD43A4_7.rds"))
b7 <- subset(b7,SSBS %in% md$newid)

# Reclassify back to correct names
b1$SSBS <- md$missing[match(b1$SSBS,md$newid)]
b2$SSBS <- md$missing[match(b2$SSBS,md$newid)]
b3$SSBS <- md$missing[match(b3$SSBS,md$newid)]
b4$SSBS <- md$missing[match(b4$SSBS,md$newid)]
b5$SSBS <- md$missing[match(b5$SSBS,md$newid)]
b6$SSBS <- md$missing[match(b6$SSBS,md$newid)]
b7$SSBS <- md$missing[match(b7$SSBS,md$newid)]

#### Construct history ####
# SSBS |- Bands ... and |- NDVI ...
result <- list() # final result list. Appears as if saving does not work in jobs / R 3.2
# Do the thing, Julie
for(id in unique(sites$SSBS)){
  myLog("Processing ",id)
  # Create a subset
  s <- subset(sites,SSBS == id)
  if(nrow(s)==0) stop("Something went wrong. Encoding?")
  myLog(unique(s$Source_ID)," | Start: ",unique(s$Sample_end_latest))
  
  # Construct the interval of interest
  # Subset bands to respective ids
  sub_b1 <- b1[which(b1$SSBS == id),]
  sub_b2 <- b2[which(b2$SSBS == id),]
  sub_b3 <- b3[which(b3$SSBS == id),]
  sub_b4 <- b4[which(b4$SSBS == id),]
  sub_b5 <- b5[which(b5$SSBS == id),]
  sub_b6 <- b6[which(b6$SSBS == id),]
  sub_b7 <- b7[which(b7$SSBS == id),]
  # Finally aggregate all
  sub <- rbind(sub_b1,sub_b2)
  sub <- rbind(sub,sub_b3)
  sub <- rbind(sub,sub_b4)
  sub <- rbind(sub,sub_b5)
  sub <- rbind(sub,sub_b6)
  sub <- rbind(sub,sub_b7)
  rm(sub_b1,sub_b2,sub_b3,sub_b4,sub_b5,sub_b6,sub_b7) # clean up
  
  if(nrow(sub)==0){
    myLog("NO values found! Started in ",unique(s$Sample_start_earliest))
    next()
  }
  
  # Final preperations
  sub$value <- sub$value * 0.0001# Correct
  sub$Band <- as.factor(sub$Band)
  sub$date <- as.Date(ymd(sub$date),origin="2000-01-01")
  
  ## Now construct time-series for all bands
  for(b in seq(1:7)){
    myLog("       |- Band ",b)
    # Get band
    sub1 = sub %>% dplyr::filter(Band == b) %>% dplyr::select(value,date)
    
    # Construct zoo time series
    sub1 = zooreg(sub1$value,order.by = sub1$date)
    # Merge with ideal time series and use the merge product as new ts
    sub1 = merge(ideal,sub1,all = T)[,2]
    sub1.ori = sub1
    
    # Detect and remove outliers based on MAD defined threshold
    # http://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/
    if(outlier){
      x = abs(sub1 - median(sub1,na.rm = T)) / mad(sub1,na.rm = T)
      tr = ifelse(quantile(x,.99,na.rm=T)>2,quantile(x,.99,na.rm=T),NA) # Determine threshold
      if(!is.na(tr)){
        sub1[ union(which(x > tr),which(x > 10))] <- NA
      }
    }
    
    ## Convert to ts
    # # Split per year
    # sub1.split <- split.xts(ideal,f = "years")
    # 
    # # Get 8 days endpoints of ideal time series
    # ep <- endpoints(ideal,'days',k=8)
    # sub1_long <- period.apply(x=ideal,ep,FUN=function(x) ifelse(all(is.na(x)), NA, max(x,na.rm=T) )) # MVC
    # rm(ideal,tt,ep)
    
    # Convert to TS
    #sub1_long.ts <- as.ts(ideal)
    
    # Next fill all gaps using a kalman smoother of state-space model
    sub1_i <- try(
      na.kalman(coredata(sub1),model = "auto.arima",smooth = TRUE),
      silent=T)
    
    # The time series of NDVI was smoothed using the Savitzky–Golay filter with a length of 5 to reduce noise 
    # but still expose abrupt change events that might occur in the series (Jönsson & Eklundh, 2004)
    # Note this queries the signal package sgolayfilter function, which has other default parameters
    # that might defer from the default as set in the TIMESAT software
    sub1_i <- signal::sgolayfilt(sub1_i, n = 5)
    
    # Overwrite previous values
    coredata(sub1) <- sub1_i
    
    # Reset values
    s.na <- na.approx(sub1.ori,maxgap=5)
    sub1_i2 <- merge(sub1,s.na) # prev. cbind
    # Equalize gaps where continious gap length too long
    sub1_i2[which(is.na(apply(sub1_i2,1,mean))),1] <- NA
    subi_1 <- sub1_i2[,1];rm(sub1_i2,s.na)
    
    # Create window until sampling end + 1 day
    out <- window(subi_1,start = ymd(paste0("2000-01-01")),end = ymd(paste0( year(unique(s$Sample_end_latest)),"-12-31"))+1 )
    
    # Save 
    result[[id]][[paste0("Band",b)]] <- out
    
    # if(class(sub1_i)[1]=="try-error") { 
    #   #pa_flag = TRUE # Set to True
    #   next()
    # } else { 
    #   ## Ensure that there is a monthly values for every year since start ##
    #   # Generate all dates that should be present per month and merge with original data
    #   alldates <- data.table(date=seq.Date(as.Date(min(time(sub1_i))),as.Date(max(time(sub1_i))), by="months"))
    #   alldates$date <- as.yearmon(alldates$date)
    #   # Aggregate maximum per month
    #   df <- data.frame(date=as.POSIXct(time(sub1_i)),value=coredata(sub1_i))
    #   df$date <- as.yearmon(df$date)
    #   # Merge to find missing months since start
    #   dt <- merge(df, alldates, by="date", all=TRUE)
    #   
    #   # Recreate the zoo layer
    #   out <- zoo(dt$value,dt$date)
    #   out <- na.approx(na.aggregate(na.approx(out,maxgap=3),months,FUN=median))
    #   
    #   # Create a window for the sampling year
    #   #interval("2000-01-01",unique(s$Sample_start_earliest))
    #   out <- window(out,start = as.yearmon(ymd(paste0("2000-01-01"))),end = as.yearmon(ymd(paste0( year(unique(s$Sample_midpoint))-1,"-12-31"))) )
    #   
    #   if(length(out)>0){
    #     result[[id]][[paste0("Band",b)]] <- out
    #   }
    # }
    rm(b,sub1,sub1_i)
  }
  #if( length(result[[id]]) != 7) {
  #  result[[id]] <- NULL # remove from analyisis and skip
  #  next()
  #}
  ## Now construct vegetation indices for all 
  B1 = result[[id]][["Band1"]]
  B2 = result[[id]][["Band2"]]
  B3 = result[[id]][["Band3"]]
  B4 = result[[id]][["Band4"]]
  B5 = result[[id]][["Band5"]]
  B6 = result[[id]][["Band6"]]
  B7 = result[[id]][["Band7"]]
  
  #NDVI
  # (2-1) / (2+1)
  result[[id]][["NDVI"]] <- (B2 - B1) / (B2 + B1)
  # EVI = G * (NIR – RED)/(NIR + C1*RED - C2*BLUE + L))
  #G – Gain factor
  #L – Factor for canopy background adjustment
  #C1, C2: Coefficients for correcting aerosol influences from RED using BLUE
  #MODIS EVI algorithm: L = 1, G = 2.5, C1 = 6, C2 = 7.5
  result[[id]][["EVI"]] = 2.5 * ((B2 - B1) / (B2 + 6.0 * B1 - 7.5 * B3 + 1.0))
  # EVI2
  # Using only two bands without blue
  result[[id]][["EVI2"]] = 2.5 * ((B2 - B1) / (B2 + B1 + 1))
  # SAVI
  # SAVI = (1 + L) * (NIR – RED)/(NIR + RED + L)
  result[[id]][["SAVI"]] = (1 + 0.5) * (B2 - B1) / (B2 + B1 + 0.5)
  # NDMI
  # Gao 1996 (NIR - Swir)
  result[[id]][["NDMI"]] = (B2 - B5) / (B2 + B5)
  # GVMI 
  # Global Vegetation Moisture Index 
  #( 5 + 0.1 ) - ( 7 + 0.02 ) ( 5 + 0.1 ) + ( 7 + 0.02 )
  result[[id]][["GVMI"]] = (B5 + 0.1) - (B7 + 0.02) / (B5 + 0.1) + (B7 + 0.02)
  # NBRI -  Normalised Burn Ratio Index 
  # NDSWIR
  #  (nir - swir2)/(nir + swir2) 
  result[[id]][["NBRI"]] = (B2 - B7) / (B2 + B7)
  # Heterogeneity
  # CV
  zz = try(merge(B1,B2,B3,B4,B5,B6,B7),silent = T)
  sh = try( zoo(apply(zz,1,co.var),time(zz)),silent = T)
  if( class(sh)[1]!="try-error") result[[id]][["CV"]] <- sh
  
  d <- list(NULL, c("brightness", "greenness", "wetness"))
  tass = matrix(c(
    #Lobser & Cohen (2007) 
    0.4395, 0.5945, 0.2460, 0.3918, 0.3506, 0.2136, 0.2678, 
    -0.4064, 0.5129,-0.2744,-0.2893, 0.4882,-0.0036,-0.4169,
    0.1147, 0.2489, 0.2408, 0.3132,-0.3122,-0.6416,-0.5087), ncol = 3, dimnames = d)
  # Calculate tasseled cap transformed values
  result[[id]][["Brightness"]] = tass[1,1]*B1+tass[2,1]*B2+tass[3,1]*B3+tass[4,1]*B4+tass[5,1]*B5+tass[6,1]*B6+tass[7,1]*B7
  result[[id]][["Greenness"]] = tass[1,2]*B1+tass[2,2]*B2+tass[3,2]*B3+tass[4,2]*B4+tass[5,2]*B5+tass[6,2]*B6+tass[7,2]*B7
  result[[id]][["Wetness"]] = tass[1,3]*B1+tass[2,3]*B2+tass[3,3]*B3+tass[4,3]*B4+tass[5,3]*B5+tass[6,3]*B6+tass[7,3]*B7
  
}

myLog("Now aggregating again and creating final output!")
r.complete <- c(r.full,result)
saveRDS(r.complete,"Center_FullTimeSeriesSmooth_Final.rds")
print("Done!")

