# Load packages
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
library(greenbrown)
library(zoo)
library(xts)
library(signal) # savitzky solaay filter
library(forecast) # smoothers
library(bfast)
library(data.table)
myLog <- function(...) {
  cat(paste0("[TS_Prepare] ", Sys.time(), " | ", ..., "\n"))
}

#### Get reference years data.frame ####
# Calculate annual metrics
# Site container, band, window
# winYear = year
getAnnualMetric <- function(sites,container,band,winYear = NULL ){
  require(greenbrown)
  # Output data.frame
  result <- data.frame()
  # Save the container
  ref <- container
  # Loop through
  for(id in unique(sites$SSBS)){
    myLog(band,"---",id)
    # Get annual time series
    x = ref[[id]][[band]]
    if(length(x)<3) next()
    # Check if windows was defined
    if(!is.null(winYear)){
      # Check if this is actually relevant
      if(winYear %in% unique(year(index(x)))) 
        x = window(x,start=as.Date(paste0(winYear,"-01-01")),end=as.Date(paste0(winYear,"-12-31"))) else
        next()
    }
    # Smooth with Savitzky-Golay
    x <- zoo(sgolayfilt(x),time(x))
    # Convert to TS
    xx <- ts(x,start = c(year(time(x)[1]),1),frequency = 365)
    ## Neutral phenology curve
    # Take the pot position (minimum)
    pot <- (PhenoDeriv(xx))["pot"]
    # Rescale the time series using the start of season as new start of season
    t = zoo(c(as.vector(x[pot:length(x)]),as.vector(x[1:pot-1])),1:46)
    r <- PhenoDeriv(ts(t,frequency = 365))
    # Integrated area under the curve total
    intg_total = (sum(coredata(t)) - as.vector(r["trough"])*length(t)) / length(t) #(((365*length(years))-16)*365) # average annual yield  (i.e. work out daily yield * 365
    # Integrated area of baseline (average) and peak
    # Calculate average
    avg <- mean(coredata(t),na.rm=T)
    # Extract window of interest woi
    woi <- window(t,start = min(which(t>=avg)),end=max(which(t>=avg)))
    intg_aavg = (sum(coredata(woi)) - as.vector(min(woi))*length(woi)) / length(woi) #(((365*length(years))-16)*365) # average annual yield  (i.e. work out daily yield * 365
    woi_f1 <- window(woi,start = index(woi)[1],end=which.max(woi))
    intg_grow = (sum(coredata(woi_f1)) - as.vector(min(woi_f1))*length(woi_f1)) / length(woi_f1)
    woi_f2 <- window(woi,start = which.max(woi)+1,end=max(index(woi)))
    intg_sensec = (sum(coredata(woi_f2)) - as.vector(min(woi_f2))*length(woi_f2)) / length(woi_f2)
    # Define output and save
    d <- data.frame(SSBS = id,band=band,year = year(time(x))[1],
                    avg = avg, intg_total = intg_total, intg_aavg = intg_aavg,intg_grow = intg_grow,intg_sensec = intg_sensec,
                    peak = r["peak"], through = r["trough"], mgs = r["mgs"],rsp = r["rsp"],rau = r["rau"])
    
    result <- rbind(result,d)
    rm(d)
  }
  return(result)
}

buildRefMetrics <- function(){
  # Load all precompiled data (Prep scripts 01-02b)
  sites <- readRDS("PlantStudiesOnFocus.rds")
  refy <- readRDS("ReferenceYearUnSmoothed.rds")

  # Calculate for every f****** band
  res <- rbind_all( lapply(names(refy[[1]]), function(x) getAnnualMetric(sites,refy,x) ) )
  saveRDS(res,"ReferenceYearMetrics.rds")
}

buildAllMetrics <- function(){
  # Load all precompiled data (Prep scripts 01-02b)
  sites <- readRDS("PlantStudiesOnFocus.rds")
  specPast <- readRDS("FullTimeSeriesUnSmoothed.rds")
  
  # Calculate for every f****** band
  res1 <- rbind_all( lapply(names(specPast[[1]]), function(x) getAnnualMetric(sites,specPast,band = x) ) )
  saveRDS(res,"AllAnnualMetrics.rds")
}

annual_metric <- function(x,id,band,winYear=NULL){
  require("greenbrown")
  myLog(id, "-" ,winYear)
  ## require a annual metric
  if(length(x)<3) next()
  # Check if windows was defined
  if(!is.null(winYear)){
    # Check if this is actually relevant
    if(winYear %in% unique(year(index(x)))) 
      x = window(x,start=as.Date(paste0(winYear,"-01-01")),end=as.Date(paste0(winYear,"-12-31"))) else
        next()
  }
  # Smooth with Savitzky-Golay
  x <- zoo(sgolayfilt(x),time(x))
  # Convert to TS
  xx <- ts(x,start = c(year(time(x)[1]),1),frequency = 365)
  ## Neutral phenology curve
  # Take the pot position (minimum)
  pot <- (PhenoDeriv(xx))["pot"]
  # Rescale the time series using the start of season as new start of season
  t = zoo(c(as.vector(x[pot:length(x)]),as.vector(x[1:pot-1])),1:46)
  r <- PhenoDeriv(ts(t,frequency = 365))
  # Integrated area under the curve total
  intg_total = (sum(coredata(t)) - as.vector(r["trough"])*length(t)) / length(t) #(((365*length(years))-16)*365) # average annual yield  (i.e. work out daily yield * 365
  # Integrated area of baseline (average) and peak
  # Calculate average
  avg <- mean(coredata(t),na.rm=T)
  # Extract window of interest woi
  woi <- window(t,start = min(which(t>=avg)),end=max(which(t>=avg)))
  intg_aavg = (sum(coredata(woi)) - as.vector(min(woi))*length(woi)) / length(woi) #(((365*length(years))-16)*365) # average annual yield  (i.e. work out daily yield * 365
  woi_f1 <- window(woi,start = index(woi)[1],end=which.max(woi))
  intg_grow = (sum(coredata(woi_f1)) - as.vector(min(woi_f1))*length(woi_f1)) / length(woi_f1)
  woi_f2 <- window(woi,start = which.max(woi)+1,end=max(index(woi)))
  intg_sensec = (sum(coredata(woi_f2)) - as.vector(min(woi_f2))*length(woi_f2)) / length(woi_f2)
  # Define output and save
  d <- data.frame(SSBS = id,band=band,year = year(time(x))[1],
                  avg = avg, intg_total = intg_total, intg_aavg = intg_aavg,intg_grow = intg_grow,intg_sensec = intg_sensec,
                  peak = r["peak"], through = r["trough"], mgs = r["mgs"],rsp = r["rsp"],rau = r["rau"])
  
  return(d)  
}