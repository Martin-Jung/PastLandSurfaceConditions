library(dplyr)
library(lubridate)
# Load data
load("HPC_o2000studies.RData")
stopifnot(exists("rr"))
stopifnot(exists("s"))

myLog <- function(...) {
  cat(paste0("[TimeSeries] ", Sys.time(), " | ", ..., "\n"))
}

source("000_HelperFunction.R")


## Calculate turnover permutations
n = 100
metric = "BCVeg"
binary = FALSE # Only for distinction of Bray and Sor
path = "" # Set to "" if no specific folder set


# Bodymass corrections
if(metric =="SorBM"){
  load("AmnioteLBM_PREDICTS_raw.rdata")
  r$Measurement <- r$log10_bodymass # Overwrite Measurement with Bodymass
  rr <- r 
  rm(r)
  # Also subset available studies to only the ones in the amniote database
  s <- subset(s,TGrouping %in% c("Aves","Mammalia","Reptilia"))
}

# -------------------------------------------------#
for(i in 1:n){
  myLog("-----",i,"-----")
  # First select one random site per cell within study
  ss <- singleRandomCellPick(s)
  # Subset existing full set to those sites
  rrr <- subset(rr,SSBS%in%ss$SSBS)
  
  # Get studies with more than only one site
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
    # Kickout empty studies ?
    if( nrow(m) == 0){
      warning("Empty study in ", study)
      next()
    }
    
    # Permute / rotate the matrix
    p <- sample.int(dim(m)[1])
    m <- m[p, p]
    # Kick first row out
    m <- m[-1,]
    
    # Take subdiagonal
    if(class(m)=="numeric"){
      # If only two sites
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
