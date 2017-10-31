library(gdalUtils)
library(raster)
library(stringr)
library(zoo)
#library(rts)

p = "../../MCD43A4_Year0/"
f = sort( list.files(p,"*.tif",full.names = T,ignore.case = T) )
years = c("2007")
s <- stack(f) # Load Stack
# Bray curtis function
bc <- function(x,y){
  zz = merge(x,y)
  if(length(which(is.na(apply(zz,1,mean))))>23) return(NA) # Harmonize and return NA if too few overlapping
  return( # bray-curtis
    sum( abs(zz[,1] - zz[,2]),na.rm = T) / (sum(zz[,1],na.rm = T) + sum(zz[,2],na.rm=T))
    ) 
  }
emptyraster <- function(x, ...) { # add name, filename, 
  
  emptyraster <- raster(nrows=nrow(x), ncols=ncol(x),
                        crs=x@crs, 
                        ext=extent(x), ...)
  
  return(emptyraster)
}

# ---------------------- #
# Process row by row
pb <- txtProgressBar(min = 1,max = nrow(s),initial = 1,style = 3)
template <- emptyraster(s[[1]])
out <- writeStart(template, filename='dBC.tif', format='GTiff', overwrite=TRUE)
ngb = 3 # nbg

for(rr in 1:nrow(s)){
  mat <- raster::getValuesFocal(s,row = rr,nrows = 1,ngb = ngb,names=TRUE,padValue=NA,array=T)
  mat[mat<0] <- NA # Reset values smaller than zero to NA
  
  # Now for every cell (also 2400)
  rowVal = vector()
  for(cc in 1:nrow(mat)){
    # Calculate reference time series
    ref = zoo::na.approx( zooreg(mat[cc,5,]),maxgap=5 )
    # If more than half is missing skip    
    if(length(which(is.na(ref)))>23){
      rowVal <- c(rowVal,NA) # Append empty value
      next()
    } 
    if(length(ref)==0){ # If empty (nothing to approx)
      rowVal <- c(rowVal,NA) # Append empty value
      next()
    } 
    ref.bc = vector()
    for(o in c(1,2,3,4,6,7,8,9)){
      test = zoo::na.approx(zooreg(mat[cc,o,]),maxgap=5)
      if(length(test)==0) next()
      ref.bc <- c(ref.bc,
                  bc(ref,test)
      )
    }
    # Average those over zero. Zero occurs if test is empty
    val <- mean(ref.bc[which(ref.bc!=0)],na.rm=T)
    rowVal <- c(rowVal,val)
  }
  stopifnot(length(rowVal)==nrow(mat)) # ThiS should be as long as the matrix
  
  out <- writeValues(out,rowVal,rr)
  setTxtProgressBar(pb,rr)# Update progressbar
  rm(rowVal) # Delete computd values to be sure
  
}
close(pb)
out <- writeStop(out)
stop("DONE!!")
