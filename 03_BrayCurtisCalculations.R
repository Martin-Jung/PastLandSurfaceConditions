# Package loading
par.ori <- par(no.readonly = T)
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
# Custom packages
library(marfunky)
library(ggplot2)
library(scales)
library(TSdist)
library(zoo)
library(xts)
#library(piecewiseSEM)
myLog <- function(...) {
  cat(paste0("[TimeSeries] ", Sys.time(), " | ", ..., "\n"))
}
#source("02c_GenerateAnnualMetrics.R")
sites <- readRDS("AllPredictsSites_full.rds") # All
sites <- subset(sites, startyear > 2001)

#refy <- readRDS("Center_ReferenceYearMonthlyFilled.rds")
specPast <- readRDS("Center_FullTimeSeriesSmooth_Final.rds")

#### Hypothesis - Distance timeseries vs distance biodiversity ####
# -------------------------------------------------#
R.utils::detachPackage(c("greenbrown","plyr"))
source("000_HelperFunction.R")
# Are pairwise distance betwen timeseries related to pairwise difference between biodiversity values?
names(specPast[[1]])
band = "EVI2"
metric = "dBC"
matr = TRUE # Save a raw matrix?
zstand = FALSE # Should absolute z-score standardization be used?
ylist = c(0,1,2,3,4,5)#ylist=c(5)
# START PROCESSING # 
# ---------------------------------------------------- #
# Vector of maximal historic conditions (-1 ... -10)
for(yrs in ylist ){
  # Filter study sites if necessary
  if(yrs > 0){
    focal.period  = ymd("2000.02.18") + (365*(yrs+1)) # Available past
    s <- subset(sites,Sample_start_earliest > focal.period) # Suitable study sites
    
  } else {s <- sites} #Do not need to subset as all time series are greater than a year
  myLog("Assessing period -",yrs, "| NR of studies = ",length(unique(s$SS)))
  
  # subset to suitable studies
  py <- specPast[which(names(specPast) %in% unique(s$SSBS))] 
  py_short <- list() # Shortened time series
  
  # Subset them to target band 
  for(id in unique(s$SSBS)){
    x = py[[id]][[band]]
    if(is.null(x)) next()
    x[which(x < 0 )] <- NA # Overwrite snow with zeros
    py_short[[id]] <- x
  }
  rm(py,id,x,xx) # cleanup
  
  myLog("Calculate pairwise differences - Past period -",yrs)
  ## Calculate pairwise time-series distances per study
  res_yearbefore <- list() # A between sites the year before year of sampling
  for(study in unique(s$SS) ){
    myLog(study)
    sub <- subset(s,SS == study) %>% mutate(SSBS = as.character(SSBS))
    n = length(sub$SSBS)
    if(n==1) next() # Single study site kinda uninteressting
    # Very inefficient code, but it works for now
    mat <- matrix(nrow = n, ncol = n,dimnames = list(sub$SSBS,sub$SSBS))
    for(i in 1:n){
      for(j in 1:n){
        x = py_short[[sub$SSBS[i]]] # Site x
        y = py_short[[sub$SSBS[j]]] # Site y
        
        if(length(x) == 0 | length(y) == 0 ) next() # Skip if one of them has no values (only NA)
        # If both fall into the same cell, set to NA (no distance)
        if(sub$cells[i]==sub$cells[j]) {
          mat[i,j] <- NA
        } else {
          # Construct sampling window 
          # Minimum pairwise consensus date
          sdate <- min(c(sub$Sample_start_earliest[i],sub$Sample_start_earliest[j]))
          # Also skip if they differ further than 3 months
          if( as.numeric(abs(sub$Sample_start_earliest[i] - sub$Sample_start_earliest[j])) > 90) next()
          # Sampling conditions sampled or past period?
          if(yrs == 0){ # yrs0
            earliest <- (sdate) - 365 # Get first year before sampling start
            x <- window(x,start = earliest, end = sdate,extend=F)
            y <- window(y,start = earliest, end = sdate,extend=F)
            
          } else{ # Past conditions
            earliest <- sdate - ((yrs+1) * 365) # get past period start from yr0 off
            x <- window(x,start = earliest, end = sdate - 365,extend=F)# minus yrs 0
            y <- window(y,start = earliest, end = sdate - 365,extend=F)# minus yrs 0
          }
          
          # Merge them and take only the overlapping part
          zz = merge(x,y)
          # Overwrite all with NA to NA (equalize gaps)
          non.perfect.overlap <- which(is.na(apply(zz,1,mean)))
          zz[non.perfect.overlap,1:2] <- NA
          
          # Calculate missing data proportion
          missing = length( which(is.na(zz[,1])) ) / length(zz[,1])
          # If missing is over 50% then skip
          if(missing > .50) next()
          
          if(metric == "correctedBC"){
            o <- (sum( abs(zz[,1] - zz[,2]),na.rm = T) / (sum(zz[,1],na.rm = T) + sum(zz[,2],na.rm=T)) )  * mean(c(coredata(zz[,1]),coredata(zz[,2])),na.rm=T)
          } else{
            # Calculate the bray
            o <- ( sum( abs(zz[,1] - zz[,2]),na.rm = T) / (sum(zz[,1],na.rm = T) + sum(zz[,2],na.rm=T)) ) # bray-curtis
          }
          mat[i,j] <- o
          rm(o,missing,non.perfect.overlap,zz) # Clean up 
        }
        
      }
    }
    if(!matr){
      mat[upper.tri(mat)] <- -1; diag(mat) <- -1  # Remove duplicates and set diagonal to -1
      mat <- reshape2::melt(mat) %>% dplyr::rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
        mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
      mat <- mat[-which(mat$value<0),] # Remove values smaller than 0 - duplicates (not possible with manhattan)
      mat <- dplyr::rename(mat,RS = value)
      mat$SS <- study # Finally append the study ID
    }
    res_yearbefore[[study]] <- mat
    
  }
  saveRDS(res_yearbefore,paste0("res_yearbefored",ifelse(metric=="correctedBC","correctedBC","BC"),"_minStart_",band,paste0(ifelse(yrs == 0,"samplingperiod","pastperiod"),"_",yrs,"y"),ifelse(matr,"-matrix",""),".rds"));rm(res_yearbefore)
  myLog("Done - Past period -",yrs)
  
}
# DONE
myLog("-- !!!DONE!!! --")

# Assess how many are stationary ?
#tseries::adf.test(coredata(na.approx(x)),"stationary")

# ----------# 
# ## Vector of historic years BY Year!
# for(yrs in c(0,1,2,3,4,5) ){
#   # Filter study sites if necessary
#   if(yrs > 0){
#     focal.period  = ymd("2000.02.18") + (365*(yrs+1)) # Available past
#     s <- subset(sites,Sample_start_earliest > focal.period) # Suitable study sites
#     
#   } else {s <- sites} #Do not need to subset as all time series are greater than a year
#   myLog("Assessing period -",yrs, "| NR of studies = ",length(unique(s$SS)))
#   
#   # subset to suitable studies
#   py <- specPast[which(names(specPast) %in% unique(s$SSBS))] 
#   py_short <- list() # Shortened time series
#   
#   # Subset them to target band 
#   for(id in unique(s$SSBS)){
#     x = py[[id]][[band]]
#     if(is.null(x)) next()
#     x[which(x < 0 )] <- NA # Overwrite snow with zeros
#     py_short[[id]] <- x
#   }
#   rm(py,id,x,xx) # cleanup
#   
#   myLog("Calculate pairwise differences - Past period -",yrs)
#   ## Calculate pairwise time-series distances per study
#   res_yearbefore <- list() # A between sites the year before year of sampling
#   for(study in unique(s$SS) ){
#     myLog(study)
#     sub <- subset(s,SS == study) %>% mutate(SSBS = as.character(SSBS))
#     n = length(sub$SSBS)
#     if(n==1) next() # Single study site kinda uninteressting
#     # Very inefficient code, but it works for now
#     mat <- matrix(nrow = n, ncol = n,dimnames = list(sub$SSBS,sub$SSBS))
#     for(i in 1:n){
#       for(j in 1:n){
#         x = py_short[[sub$SSBS[i]]] # Site x
#         y = py_short[[sub$SSBS[j]]] # Site y
#         
#         if(length(x) == 0 | length(y) == 0 ) next() # Skip if one of them has no values (only NA)
#         # If both fall into the same cell, set to NA (no distance)
#         if(sub$cells[i]==sub$cells[j]) {
#           mat[i,j] <- NA
#         } else {
#           # Construct sampling window 
#           # Minimum pairwise consensus date
#           sdate <- min(c(sub$Sample_start_earliest[i],sub$Sample_start_earliest[j]))
#           # Also skip if they differ further than 3 months
#           if( as.numeric(abs(sub$Sample_start_earliest[i] - sub$Sample_start_earliest[j])) > 90) next()
#           # Sampling conditions sampled or past period?
#           if(yrs == 0){ # yrs0
#             earliest <- (sdate) - 365 # Get first year before sampling start
#             x <- window(x,start = earliest, end = sdate,extend=F)
#             y <- window(y,start = earliest, end = sdate,extend=F)
#             
#           } else{ # Past conditions
#             earliest <- sdate - ((yrs+1) * 365) # get past period start from yr0 off
#             x <- window(x,start = earliest, end = sdate - (yrs)*365,extend=F)# minus yrs 0
#             y <- window(y,start = earliest, end = sdate - (yrs)*365,extend=F)# minus yrs 0
#           }
#           
#           # Merge them and take only the overlapping part
#           zz = merge(x,y)
#           # Overwrite all with NA to NA (equalize gaps)
#           non.perfect.overlap <- which(is.na(apply(zz,1,mean)))
#           zz[non.perfect.overlap,1:2] <- NA
#           
#           # Calculate missing data proportion
#           missing = length( which(is.na(zz[,1])) ) / length(zz[,1])
#           # If missing is over 50% then skip
#           if(missing > .50) next()
#           
#           if(metric == "correctedBC"){
#             o <- (sum( abs(zz[,1] - zz[,2]),na.rm = T) / (sum(zz[,1],na.rm = T) + sum(zz[,2],na.rm=T)) )  * mean(c(coredata(zz[,1]),coredata(zz[,2])),na.rm=T)
#           } else{
#             # Calculate the bray
#             o <- ( sum( abs(zz[,1] - zz[,2]),na.rm = T) / (sum(zz[,1],na.rm = T) + sum(zz[,2],na.rm=T)) ) # bray-curtis
#           }
#           mat[i,j] <- o
#           rm(o,missing,non.perfect.overlap,zz) # Clean up 
#         }
#         
#       }
#     }
#     if(!matr){
#       mat[upper.tri(mat)] <- -1; diag(mat) <- -1  # Remove duplicates and set diagonal to -1
#       mat <- reshape2::melt(mat) %>% dplyr::rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
#         mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
#       mat <- mat[-which(mat$value<0),] # Remove values smaller than 0 - duplicates (not possible with manhattan)
#       mat <- dplyr::rename(mat,RS = value)
#       mat$SS <- study # Finally append the study ID
#     }
#     res_yearbefore[[study]] <- mat
#     
#   }
#   saveRDS(res_yearbefore,paste0("res_yearbefored",ifelse(metric=="correctedBC","correctedBC","BC"),"_minStart_",band,paste0(ifelse(yrs == 0,"samplingperiod","pastperiod"),"_INDIVIDUAL_",yrs,"y"),ifelse(matr,"-matrix",""),".rds"));rm(res_yearbefore)
#   myLog("Done - Past period -",yrs)
#   
# }
# 
# # DONE
# myLog("-- !!!DONE!!! --")
# # ----------# 


####  Turnover biodiversity calculations single and permutations ####
# NOTE: This has been outsourced to a seperate script to run on a cluster
source("000_HelperFunction.R")
library(roquefort)
r <- readRDS("../../Data/diversity-2016-02-03-03-37-46.rds") 
r <- DropInvalidMetricsAndMethods(r)
r <- CorrectSamplingEffort(r)
s <- readRDS("AllPredictsSites_full.rds")
s$midyear = year(s$Sample_midpoint)
s <- subset(s,startyear > 2001)

rr <- subset(r,SS%in%s$SS)
rm(r)
rr$SS <- droplevels(rr$SS)
# TEST #
# Normal turnover (adapt code)
#rrr <- subset(rr,SS =="SE1_2006__Armbrecht 1")
rrr <- subset(rr,SS %in% unique(s$SS[which(s$TGrouping=="Reptilia")]))
s.metric1 <- CompDissim2(rrr,"Sor",binary=F) # dBC_abundance computed using vegdist

## Calculate turnover permutations
n = 100
metric = "BCVeg" #"BCVeg"
binary = TRUE
path = "" # Set to "" if no specific folder set
for(i in 1:n){
  myLog("-----",i,"-----")
  # First select one random site per cell within study
  ss <- singleRandomCellPick(s)
  # Subset existing full set to those sites
  rrr <- subset(rr,SSBS%in%ss$SSBS)
  
  # Get studies with only one site
  ag <- aggregate(list(SSBS=rrr$SSBS),by=list(SS = rrr$SS),function(x) length(unique(x)))
  studies.in <- ag$SS[which(ag$SSBS>1)]
  
  rrr <- subset(rrr,SS %in% studies.in)
  
  # Generate matrices
  s.metricm <- CompDissim2(rrr,metric,binary=binary)
  
  # Process list
  s.metric4 <- data.frame()
  for(study in names(s.metricm)){
    myLog(study)
    m <- s.metricm[[study]]
    # Permute / rotate the matrix
    p <- sample.int(dim(m)[1])
    m <- m[p, p]
    # Kick first row out
    m <- m[-1,]
    
    # Take subdiagonal
    if(class(m)=="numeric"){
      val <- m[1]
      s.x <- names(m)[1]
      s.y <- names(m)[2]
    } else {
      val <-   diag(m)
      s.x <- rownames(m)
      s.y <- colnames(m)[-ncol(m)]
    }
    s.metric4 <- rbind(s.metric4,data.frame(
      SS = study,
      SSBS_x = s.x, SSBS_y = s.y,
      value = val
    ))
    
  }
  saveRDS(s.metric4,paste0(path,"s_metric",metric,"_permute_",i,".rds"))
  myLog(" ","DONE"," ")
}
#### Simulation - Testing time-series dissimilarity metrics ####
# Simulate artificial time series
library(greenbrown)
source("000_HelperFunction.R")

# Calculate distance matrix
calcSim <- function(x,y){
  # Calculate distance matrices
  man = sum( abs(x - y),na.rm = T ) # manhattan
  can = sum( abs(x - y) / (x+y) ,na.rm = T) # canberra
  braycurt = 1 - ( sum( abs(x - y)) / (sum(x)+sum(y)) ) # bray-curtis
  DIV =   as.vector( abs( simpson(coredata(x),1,length(index(x))) - simpson(coredata(y),1,length(index(y))) )  )
  Euc = TSDistances(x,y,distance = "euclidean")
  SorAbd = (2 * sum(x) * sum(y)) / (sum(x) + sum(y))

  return( data.frame(Manhattan = man,
                     Euclidean = Euc,
                    Canberra = can,
                    BrayCurtis = braycurt,
                    SorAbundance = SorAbd,
                    AbsoluteIntegral = DIV)
          )
}
result <- data.frame() # result
# Parameters
run = 200 # How many random generated time series?

# Just random generated time-series
for(i in seq(1,run)) {
  myLog("Run : ", i)
  # Higher time-series
  x <- SimTs(M=sample(seq(0.55,0.85,by = 0.05),size = 1),
             Tslope=0, Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years
  # Second
  y <- SimTs(M=sample(seq(0.1,0.5,by = 0.05),size = 1),
             Tslope=0, Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years

  # Save result
  result <- rbind(result,
                  data.frame(run = i,
                             type = "StationarySeperateTimeSeries",
                             calcSim(x,y) # Append results
                             )
                  )
  # -- -- #
  # Overlapping time-series stationary
  # One of them being more variable than the other
  x <- SimTs(M=0.65, Tslope=0, Isd=0.015, Irange=0.03, Srange=0.5, Rsd=0.05, 
             Rrange=0.1, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years
  # Second
  y <- SimTs(M=0.60, Tslope=0, Isd=0.45+(i/1000), Irange=0.05, Srange=0.5, Rsd=0.5, 
             Rrange=0.5+(i/1000), breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years

  # Save result
  result <- rbind(result,
                  data.frame(run = i,
                             type = "StationaryOverlappingIncVariance",
                             calcSim(x,y) # Append results
                  )
  )


  # -- -- #
  # One time series with random slope 
  x <- SimTs(M=0.65, Tslope=sample(seq(-0.01,0.01,by = 0.001),1),
             Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.25, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years
  # Second
  y <- SimTs(M=0.25, Tslope=0, Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.25, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years

  # Save result
  result <- rbind(result,
                  data.frame(run = i,
                             type = "GradualChangesTrend",
                             calcSim(x,y) # Append results
                  )
  )
  
  # -- -- #
  # Seperate time-series stationary
  # increasing variability
  x <- SimTs(M=0.50, Tslope=c(-0.01,-0.05,0.005), Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=c(20,40), abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years
  # Second
  y <- SimTs(M=0.45, Tslope=0, Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years

  # Save result
  result <- rbind(result,
                  data.frame(run = i,
                             type = "OverlapBreakpointRecovery",
                             calcSim(x,y) # Append results
                  )
  )
  
  # -- -- #
  # Introduce breakpoint
  x <- SimTs(M=0.65, Tslope=c(-0.025,0.015), Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=sample(10:40,1), abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years
  # Second
  y <- SimTs(M=0.25, Tslope=0, Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years  

  result <- rbind(result,
                  data.frame(run = i,
                             type = "SeperateSingleBreakpointRecovery",
                             calcSim(x,y) # Append results
                  )
  )
  
  # -- -- #
  # Two breakpoints
  x <- SimTs(M=0.5, Tslope=sample(seq(-0.01,0.01,by=0.0005),7), Isd=0.15, Irange=0.25, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=sample(10:55,5), abrupt=F,
             n=60, start=c(2002, 1), freq=12) # past 5 years
  # Second
  y <- SimTs(M=0.5, Tslope=0, Isd=0.015, Irange=0.025, Srange=0.15, Rsd=0.05, 
             Rrange=0.1, breaks=NULL, abrupt=T,
             n=60, start=c(2002, 1), freq=12) # past 5 years  
  


  zz = na.contiguous(merge(zoo(x[,1]),zoo(y[,1])))

  result <- rbind(result,
                  data.frame(run = i,
                             type = "CompleteChaos",
                             calcSim(zz[,1],zz[,2]) # Append results
                  )
  )
  #plot(x[,1],ylim=c(0,1),col="darkgreen")
  #lines(y[,1],col="red")

}

# Visualize performance
library(GGally)


dir.create("Simulation_Distance")

for(i in 1:6){
  val = unique(result$type)[i]
  png(paste0("Simulation_Distance/",val,".png"),width=1024,height=786)
  g <- result %>% dplyr::filter(type == val) %>% 
    dplyr::select(-run,-type) %>% 
    ggpairs(.,title = val) + theme_bw()
  print(g)
  dev.off()
}

# -------------------------------------------------#