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
library(bfast)
library(data.table)
myLog <- function(...) {
  cat(paste0("[TS_Prepare] ", Sys.time(), " | ", ..., "\n"))
}
co.var <- function(x) ( sd(x, na.rm = T)/mean(x, na.rm = T))
# All Predicts sites
sites <- readRDS("AllPredictsSites.rds")
sites$midyear = year(sites$Sample_midpoint)
sites <- subset(sites,midyear > 2005)


#### Best possible gap-filling approach experiment ####
# Requires non-smoothed original time-series
# Test for best observation technique
# Use complete time-series
# take 49% out and approx
# than plot approx against real

# Quick visualization
r <- readRDS("FullTimeSeriesSmoothedCorrected.rds")
missing = unlist( lapply(r,function(x) length(which(is.na(x[["EVI2"]]))) / length(x[["EVI2"]])) )
id <- names(missing)[which(missing==0)]

# Experiment start
seed = 1024 # starting seed for reproduciability
run = 6 # Runs
propRange = round(seq(.05,.60,length.out = run),2) # Take from 5% to 45% out
res <- data.frame() # output 

for(i in 1:run) {
  set.seed(seed+i)
  prop = propRange[i]
  for(tid in id) {
    myLog(i," - ",tid)
    x  = r[[tid]]$EVI2
    # Make a proportional sample
    sam <- sample(time(x),round( prop * length(x) ))
    
    # Now replace with NA
    x_ori <- x
    x[sam] <- NA
    
    # Pass if time-series is shorter than 10
    if( length(x) < 10) next()
    
    # Define window to analyse
    win <- window(x_ori,start = time(na.trim(x))[1],end= time(na.trim(x))[length(na.trim(x))])
    
    obs <- win[sam] # Get observed values
    x <- na.trim(x,'both')
    
    # And fill
    # Linear
    pred_l <- na.approx(x)[sam]
    # Spline
    pred_sl <- na.spline(x)[sam]
    # Median of annual values and then overall
    pred_ml <- na.approx( na.aggregate(x,by = lubridate::yday,na.rm = T,FUN = median))[sam]
    # Kalman smoother
    fit <- (StructTS(ts(x,frequency = 8),type = "BSM"))
    ks <- KalmanSmooth(coredata(x),fit$model)$smooth[,1]
    pred_ks <- zoo(ks,time(win))[sam]
    
    # Current approach
    xx <- na.approx(na.aggregate(na.approx(x,maxgap=6),by = lubridate::yday,na.rm = T,FUN = median) )
    pred_lml <- xx[sam]
    # Save output if suitable
    if(length(obs)==0 | length(obs) == 1) next()
    res <- rbind(res,data.frame(
      SSBS = tid,
      prop = prop,
      Run = paste0("Run: ",i," \n Missing Data: ",prop),
      TL = length(x),
      rsq_linear = summary(lm(obs~pred_l))$r.squared,
      rsq_spline = summary(lm(obs~pred_sl))$r.squared,
      rsq_medannual = summary(lm(obs~pred_ml))$r.squared,
      rsq_kalman = summary(lm(obs~pred_ks))$r.squared,
      rsq_mixedcur = summary(lm(obs~pred_lml))$r.squared
    ))
  }
}

o <- res %>% dplyr::select(SSBS,prop,rsq_linear:rsq_mixedcur) %>% reshape2::melt(id.vars=c("SSBS","prop"))
levels(o$variable) <- c("Linear","Spline","MedianAverageDay","SeasonalKalmanFilter","CurrentApproach")

# Calculate average per method
oo <- o %>% dplyr::group_by(prop,variable) %>% 
  dplyr::summarise(avg = mean(value,na.rm=T),
                   se.min = mean(value,na.rm=T) - sd(value,na.rm=T),
                   se.max = mean(value,na.rm=T) + sd(value,na.rm=T) ) %>% 
  as.data.frame()

o %>% dplyr::group_by(variable) %>% dplyr::summarise(avg=mean(value,na.rm=T))

ggplot(aes(x=prop,y=avg,ymin=se.min, ymax = se.max),data=oo) + 
  geom_pointrange() + facet_wrap(~variable,nrow = 1) +
  labs(x = "Proportion of data missing", y= expression(R^2)) +
  theme_light() + scale_x_continuous(breaks=pretty_breaks(4)) + scale_y_continuous(breaks=pretty_breaks(5))


#### Maximum length of gap experiment ####

# Experiment gap width
# create gaps at constant positions (.25,50,75) of length of time series
# stepwise increase gap-width from 1-10
# which approach performs best in terms of varying gap-width

# Quick visualization
r <- readRDS("FullTimeSeriesSmoothedCorrected.rds")
missing = unlist( lapply(r,function(x) length(which(is.na(x[["EVI2"]]))) / length(x[["EVI2"]])) )
id <- names(missing)[which(missing==0)]

# Experiment start
seed = 1024 # starting seed for reproduciability
maxGapLength = 16 # Take from 5% to 45% out
res <- data.frame() # output 

for(i in 2:maxGapLength) {
  set.seed(seed+i)
  gapl <- i
  for(tid in id) {
    myLog(gapl," - ",tid)
    x  = r[[tid]]$EVI2
    # Make a proportional sample
    pos <- round(sample(2:(length(x)-gapl-1),1))
    if(pos <= 1 | pos == length(x)) pos <- round(sample(2:(length(x)-gapl-1),1)) # reroll dice
    sam <- seq( pos,pos+gapl) # middle
    
    # Now replace with NA
    x_ori <- x
    x[sam] <- NA
    
    # Pass if time-series is shorter than 10
    if( length(x) < 10) next()
    
    obs <- x_ori[sam] # Get observed values
    x <- na.trim(x,'both')
    
    # And fill
    # Linear
    pred_l <- na.approx(x)[sam]
    # Spline
    pred_sl <- na.spline(x)[sam]
    # Median of annual values and then overall
    pred_ml <- na.approx( na.aggregate(x,by = lubridate::yday,na.rm = T,FUN = median))[sam]
    # Kalman smoother
    fit <- (StructTS(ts(x,frequency = 8),type = "BSM"))
    ks <- KalmanSmooth(coredata(x),fit$model)$smooth[,1]
    pred_ks <- zoo(ks,time(x))[sam]
    
    # Current approach
    xx <- na.approx(na.aggregate(na.approx(x,maxgap=6),by = lubridate::yday,na.rm = T,FUN = median) )
    pred_lml <- xx[sam]
    # Save output if suitable
    if(length(obs)==0 | length(obs) == 1 | all(is.na(pred_l)) ) next()
    res <- rbind(res,data.frame(
      SSBS = tid,
      prop = gapl,
      Run = paste0("Gaplength: ",i),
      TL = length(x),
      rsq_linear = summary(lm(obs~pred_l))$r.squared,
      rsq_spline = summary(lm(obs~pred_sl))$r.squared,
      rsq_medannual = summary(lm(obs~pred_ml))$r.squared,
      rsq_kalman = summary(lm(obs~pred_ks))$r.squared,
      rsq_mixedcur = summary(lm(obs~pred_lml))$r.squared
    ))
  }
}



o <- res %>% dplyr::select(SSBS,prop,rsq_linear:rsq_mixedcur) %>% reshape2::melt(id.vars=c("SSBS","prop"))
levels(o$variable) <- c("Linear","Spline","MedianAverageDay","SeasonalKalmanFilter","CurrentApproach")

# Calculate average per method
oo <- o %>% dplyr::group_by(prop,variable) %>% 
  dplyr::summarise(avg = mean(value,na.rm=T),
                   se.min = mean(value,na.rm=T) - sd(value,na.rm=T),
                   se.max = mean(value,na.rm=T) + sd(value,na.rm=T) ) %>% 
  as.data.frame()

o %>% dplyr::group_by(variable) %>% dplyr::summarise(avg=mean(value,na.rm=T))

ggplot(aes(x=prop,y=avg,ymin=se.min, ymax = se.max),data=oo) + 
  geom_pointrange() + facet_wrap(~variable,nrow = 1) +
  labs(x = "Gap length (randomly inserted sequence of NA)", y= expression(R^2)) +
  theme_light() + scale_x_continuous(breaks=pretty_breaks(4)) + scale_y_continuous(breaks=pretty_breaks(5))

