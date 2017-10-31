### Helper functions ####
library(dplyr)

myLog <- function(...) {
  cat(paste0("[Processing] ", Sys.time(), " | ", ..., "\n"))
}


new.freq <- function(d, freq = 46) {
  y <- as.Date(cut(range(d), "years")) + c(0, 366)
  yd <- seq(y[1], y[2], "year")
  yy <- as.numeric(format(yd, "%Y"))
  floor(freq * approx(yd, yy, xout = d)$y) / freq
}

# Load in al MODIS BRDF Bands
readInFormat <- function(x,idv = "SSBS"){
  require(data.table)
  val = as.data.frame(data.table::fread(x,showProgress = T))[,-1] %>% dplyr::select(-.geo,-latitude, -longitude) %>% 
    reshape2::melt(., id.vars = idv) %>%
    distinct() %>% mutate(variable = as.character(variable)) %>% 
    mutate(year = as.numeric( str_sub(variable,13,16)), # Get year out of file name
           month = as.numeric( str_sub(variable,18,19) ), # get month
           day = as.numeric( str_sub(variable,21,22) )) %>%  # get day
    # Make a date column
    mutate(date = ymd(paste(year,month,day,sep="-")))
  return(val)
}

# Function to convert from MODIS QA scores - BRDF Albedo band scores
# Inspired by MODISTools -> Tuck et al.
BRDFconvertQAScores <- function(qsc,band=NA,QualityThreshold=3){
  #MCD43A4 = c(0,4294967294,4294967295), # BRDF albedo band quality, taken from MCD43A2, for reflectance data
  QualityScores <- qsc
  if(max(QualityScores,na.rm = T) == 0) num.binary.digits <- 1
  if(max(QualityScores,na.rm = T) != 0) num.binary.digits <- floor(log(max(QualityScores,na.rm = T), base = 2)) + 1
  
  binary.set<- matrix(nrow = length(QualityScores), ncol = num.binary.digits)
  for(n in 1:num.binary.digits){
    binary.set[ ,(num.binary.digits - n) + 1] <- QualityScores %% 2
    QualityScores <- QualityScores %/% 2
  }
  # Construct quality binary score
  quality.binary <- apply(binary.set, 1, function(x) paste(x, collapse = ""))
  rm(binary.set)
  
  if(!is.na(band)){
    band.num <- as.numeric(substr(band, nchar(band), nchar(band)))
    
    # Select respective subset of quality score
    qa.binary <- substr(quality.binary, (nchar(quality.binary) - (((band.num - 1) * 2) + 2)),
                        (nchar(quality.binary) - ((band.num - 1) * 2)))
    
    # Make a result vector
    qa.int <- numeric(length(qa.binary))
    qa.int[qa.binary == "000"] <- 0 # best quality, full inversion (WoDs, RMSE majority good)
    qa.int[qa.binary == "001"] <- 1 # good quality, full inversion
    qa.int[qa.binary == "010"] <- 2 # Magnitude inversion (numobs >=7)
    qa.int[qa.binary == "011"] <- 3 # Magnitude inversion (numobs >=3&<7)
    qa.int[qa.binary == "100"] <- 4 # Fill value
    qa.int[qa.binary == "ANA" | is.na(qa.binary) ] <- NA # NA
    
    # Finally replace everything above the treshold with NA
    qa.int[qa.int > QualityThreshold] <- NA
    # And return
    return(qa.int)
  } else {
   # If band score is NA, it is assumed that we intend to calcualte the solar zenith angle
    
    # Solar zenith information is stored between 8-14
    qa.binary <- substr(quality.binary,8,14)
    qa.binary[qa.binary == "ANANANA"] <- NA
    
    x = strtoi(qa.binary, base = 2)
    # Which values have to high zenith values -> Encode with 1 otherwise 0
    x = ifelse(x>QualityThreshold,1,0)
  }
}

# Backtransformation of asin(sqrt)
bt.asin <- function(x) {
  z <- sin(x)^2
  return(z)
}

#### Max CCF ####
# Returns the maximal value of the crosscorrelation function,
# thus indicating the time when there is a certain lag
aMaxCCF <- function(a,b)
{
  d <- ccf(a, b, plot = FALSE, lag.max = length(a)-5)
  cor = d$acf[,,1]
  abscor = abs(d$acf[,,1])
  lag = d$lag[,,1]
  res = data.frame(cor,lag)
  absres = data.frame(abscor,lag)
  absres_max = res[which.max(absres$abscor),]
  return(absres_max)
}


#### Compdis update #####

# Requires a full predicts subset
CompDissim2 <- function (data, metric,binary=F)
{
  require(vegan)
  require(data.table)
  data <- as.data.table(data) # For speed improvement
  if (metric == "SorAbd") {
    data <- data[data$Diversity_metric_type == "Abundance", 
                 ]
  }
  if (metric == "SorCorr") {
    data <- data[((data$Diversity_metric == "abundance") & 
                    (data$Diversity_metric_unit == "individuals")), ]
  }
  if (metric == "BCVeg" & binary == F) {
    data <- data[data$Diversity_metric_type == "Abundance",]
  }
  if (metric == "BrayCurt" & binary == F) {
    data <- data[data$Diversity_metric_type == "Abundance",]
  }
  
  data <- subset(data, select = c("SS", "SSBS", "Measurement", 
                                  "Taxon_name_entered"))
  data <- na.omit(data)
  if (metric == "SorCorr") {
    study.all.int.meas <- tapply(data$Measurement, data$SS, 
                                 function(m) all(floor(m) == m))
    int.meas <- study.all.int.meas[match(data$SS, names(study.all.int.meas))]
    data <- data[int.meas, ]
  }
  
  # Results
  result <- list()
  for (st in unique(data$SS)) {
    cat(paste("Processing ", st, "\n", sep = ""))
    sub.data <- data[data$SS == st, ]
    if (metric != "SorCorr") {
      sub.data <- sub.data[sub.data$Measurement > 0, ]
    }
    if(metric == "SorVeg"){
      if(length(unique(sub.data$SSBS)) < 2) next()
      m <- reshape2::acast(data=sub.data,SSBS~Taxon_name_entered,value.var = "Measurement",fill = 0) # Fill with zero assuming absence
      sites.matrix <- as.matrix( suppressWarnings( vegan::betadiver(m,method = "sor",binary=binary) ) )

      #sites.matrix[upper.tri(sites.matrix)] <- -9999
      diag(sites.matrix) <- NA  # Set diagonal to NA
      #mat_a <- reshape2::melt(sites.matrix) %>% rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
      #  mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
      #mat_a <- mat_a[-which(mat_a$value==-9999),] # remove triangle and missing data
      #mat_a$SS <- st # Finally append the study ID
      
      result[[st]] <- sites.matrix#mat_a
      
    } else if(metric == "BCVeg"){
      if(binary==T) stop("Binary should be FALSE for BcVEG")
      if(length(unique(sub.data$SSBS)) < 2) next()
      m <- reshape2::acast(data=sub.data,SSBS~Taxon_name_entered,value.var = "Measurement",fill = 0) # Fill with zero assuming absence
      sites.matrix <- as.matrix( suppressWarnings( vegdist(m,method="bray",binary=binary,na.rm = T) ) )
      
      #sites.matrix[upper.tri(sites.matrix)] <- -9999
      diag(sites.matrix) <- NA  # Set diagonal to NA
      #mat_a <- reshape2::melt(sites.matrix) %>% rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
      #  mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
      #mat_a <- mat_a[-which(mat_a$value==-9999),] # remove triangle and missing data
      #mat_a$SS <- st # Finally append the study ID
      
      result[[st]] <- sites.matrix#mat_a
      
    } else {
      sites.matrix <- matrix(nrow = length(unique(sub.data$SSBS)), 
                             ncol = length(unique(sub.data$SSBS)))
      i1 <- 1
      for (s1 in unique(sub.data$SSBS)) {
        i2 <- 1
        for (s2 in unique(sub.data$SSBS)) {
          if (metric == "Sor") {
            # sorrensen
            u <- length(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                            s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                               s2]))
            i <- length(intersect(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                   s2]))
            sor <- (2 * i)/((2 * i) + (u - i))
          }
          else if (metric == "SorAbd") {
            # abundance corrected sorrensen
            u <- sum(sub.data$Measurement[(sub.data$SSBS == 
                                             s1) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))])
            v <- sum(sub.data$Measurement[(sub.data$SSBS == 
                                             s2) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))])
            sor <- (2 * u * v)/(u + v)
          }
          else if (metric == "SorBM") {
            # abundance corrected sorrensen
            # Sum of estimates of species recorded at both sites
            A <- sum(sub.data$Measurement[union(sub.data$SSBS == s1,sub.data$SSBS == s2) &
                                            (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))])
            
            B <- sum(sub.data$Measurement[(sub.data$SSBS == s1) & 
                                            (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                  s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                     s2]))])
            C <- sum(sub.data$Measurement[(sub.data$SSBS == s2) & 
                                            (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                  s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                     s2]))])
            sor <- (2 * A)/(2*A + B + C)
          }
          else if (metric == "BrayCurt"){
           # Bray curtis similarity 
            u <- sub.data$Measurement[(sub.data$SSBS == 
                                             s1) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))]
            v <- sub.data$Measurement[(sub.data$SSBS == 
                                             s2) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                         s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                            s2]))]
            sor <- ( sum( abs(u - v) )  / ( sum(u) + sum(v) ) ) # Similarity 1- following Faith
            # Tested and completly identical to vegan vegdist
            
          }
          else if (metric == "SorCorr") {
            # Sampling corrected sorrensen
            n <- sum(sub.data$Measurement[sub.data$SSBS == 
                                            s1])
            m <- sum(sub.data$Measurement[sub.data$SSBS == 
                                            s2])
            if ((n > 0) & (m > 0)) {
              xi <- sub.data$Measurement[sub.data$SSBS == 
                                           s1][(match(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                          s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                             s2]), sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                                 s1]))]
              yi <- sub.data$Measurement[sub.data$SSBS == 
                                           s2][(match(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                          s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                             s2]), sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                                 s2]))]
              xi[is.na(xi)] <- 0
              yi[is.na(yi)] <- 0
              f1. <- length(which((xi == 1) & (yi > 0)))
              f2. <- max(1, length(which((xi == 2) & (yi > 
                                                        0))))
              f.1 <- length(which((xi > 0) & (yi == 1)))
              f.2 <- max(1, length(which((xi > 0) & (yi == 
                                                       2))))
              p1 <- sum(xi[yi > 0]/n)
              p2 <- ((m - 1)/m) * (f.1/(2 * f.2))
              p3 <- sum(xi[yi == 1]/n)
              u <- min(1, p1 + p2 * p3)
              q1 <- sum(yi[xi > 0]/m)
              q2 <- ((n - 1)/n) * (f1./(2 * f2.))
              q3 <- sum(yi[xi == 1]/m)
              v <- min(1, q1 + q2 * q3)
              if ((u > 0) & (v > 0)) {
                sor <- (2 * u * v)/(u + v)
              }
              else {
                sor <- 0
              }
            }
            else {
              sor <- 0
            }
          }
          else {
            stop("Error: specfied dissimilarity metric is not supported")
          }
          if (s1 != s2) 
            sites.matrix[i1, i2] <- sor
          i2 <- i2 + 1
        }
        i1 <- i1 + 1
      }
      #sites.matrix[upper.tri(sites.matrix)] <- -9999
      #diag(sites.matrix) <- -9999  # Remove duplicates and set diagonal to -9999
      rownames(sites.matrix) <- unique(sub.data$SSBS); colnames(sites.matrix) <- unique(sub.data$SSBS)
      #mat_a <- reshape2::melt(sites.matrix) %>% rename(SSBS_x = Var1, SSBS_y = Var2 ) %>% #Melt and rename
      #  mutate(SSBS_x = as.character(SSBS_x), SSBS_y = as.character(SSBS_y)) # format again
      #mat_a <- mat_a[-which(mat_a$value==-9999),] # remove triangle and missing data
      #mat_a$SS <- st # Finally append the study ID
      
      
      result[[st]] <- sites.matrix
    } 
  }
  return(result)
}



simpson <- function(y, a, b, n = 100) {
  # numerical integral of y from a to b
  # using Simpson's rule with n subdivisions
  #
  # y is a function of a single variable
  # we assume a < b and n is a positive even integer
  
  n <- max(c(2*(n %/% 2), 4),na.rm=T)
  h <- (b-a)/n
  x.vec1 <- seq(a+h, b-h, by = 2*h)
  x.vec2 <- seq(a+2*h, b-2*h, by = 2*h)
  f.vec1 <- y[x.vec1]
  f.vec2 <- y[x.vec2]
  S <- h/3*(y[a] + y[b] + 4*sum(f.vec1,na.rm = T) + 2*sum(f.vec2,na.rm = T))
  return(S)
}
eudis = function(x, y) { sqrt( sum( (x-y)^2 ) ) } # define Euclidean distance

trapezoid = function(x, y) 
{ # computes the integral of y with respect to x using trapezoidal integration. 
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}


#### Subsetting functions - pairwise ####
# Subsetting function
pairwiseRandomCellPick <- function(sout,cluster=T, pr=T){
  if(cluster) {
    library(doParallel)
    cl <- makeCluster(getOption("cl.cores", 2))
    registerDoParallel(cl)
  }
  ### Picks a random combination of cells for biodiversity values
  res <- data.frame()
  # Do a progress bar
  if(pr) pb <- txtProgressBar(min = 0, max = length(unique(unique(sout$SS))), style = 3)
  
  for(study in unique(sout$SS)){
    sub <- subset(sout,SS==study)
    if(nrow(sub) < 2) next() # To the rare case of only a single study site
    
    # Kickout those where cells.x and cell.y are identical
    sub <- sub[which(sub$cells.x!=sub$cells.y),]
    if(nrow(sub)==0) next()
    # Get unique pairwise combinations
    sub <- transform(sub, cells.u = paste(pmin(cells.x,cells.y), pmax(cells.x,cells.y), sep="_")) 
    
    if(cluster){
      df2 <- parLapply(cl,split(sub, sub$cells.u),
                       function(subdf) subdf[sample(1:nrow(subdf), 1),]
      )
    } else {
      df2 <- lapply(split(sub, sub$cells.u),
                       function(subdf) subdf[sample(1:nrow(subdf), 1),]
      )
    }
    
    res <- rbind(res, do.call('rbind', df2) )
    if(pr) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
  }
  if(pr) close(pb) # close progressbar
  if(cluster) stopCluster(cl)
  return(res)
}


# Get a Single permutation round of the subdiagonal 
pairwiseSubDiagonalPermutation <- function(sout,pr=T) {
  res <- data.frame()
  # Do a progress bar
  if(pr) pb <- txtProgressBar(min = 0, max = length(unique(unique(sout$SS))), style = 3)
  
  # Take one value per cell
  for(study in unique(sout$SS)){
    sub <- subset(sout,SS==study)

    # Kickout those where cells.x and cell.y are identical
    sub <- sub[which(sub$cells.x!=sub$cells.y),]
    if(nrow(sub)==0) next()
    # Get unique pairwise combinations
    sub <- transform(sub, cells.u = paste(pmin(cells.x,cells.y), pmax(cells.x,cells.y), sep="_")) 
    
    # Sample from combinations
    df2 <- do.call('rbind', lapply(split(sub, sub$cells.u),
                  function(subdf) subdf[sample(1:nrow(subdf), 1),]) )
    
    # if(nrow(df2)>2) {
    #   df2$SSBS_c <- paste0(df2$SSBS_x,"_",df2$SSBS_y) # variable to merge on
    #   # Build pairwise matrix
    #   m <- reshape2::acast(data=df2,SSBS_y~SSBS_x,value.vdar = "value",fill = NA) # Fill with zero assuming absence
    #   # Permute / rotate the matrix
    #   #p <- sample.int(dim(m)[1])
    #   #m <- m[p, p]
    #   # take diagonal. A bit fishy to get the row.names
    #   diag(m) <- -9999 # Take diagonal
    #   m <- reshape2::melt(m)
    #   m <- m[which(m$value == -9999),];m$SSBS_c <- paste0(m$Var1,"_",m$Var2);m$Var1 <- NULL;m$Var2 <- NULL;m$value <- NULL
    #   # Merge back and remove merging columns
    #   df2 <- merge.data.frame(m,df2,by="SSBS_c")
    #   df2$SSBS_c <- NULL; df2$cells.u <- NULL
    # }
    # df2$cells.u <- NULL
    res <- rbind(res,df2)
    if(pr) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    
  }
  if(pr) close(pb) # close progressbar
  return(res)
}

singleRandomCellPick <- function(sout){
  ### Picks a random combination of cells for biodiversity values
  res <- data.frame()
  # Do a progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(unique(sout$SS))), style = 3)
  
  for(study in unique(sout$SS)){
    sub <- subset(sout,SS==study)
    # For each cell individual cell pick a random entry
    df2 <- lapply(split(sub, sub$cells),
                  function(subdf) subdf[sample(1:nrow(subdf), 1),]
    )
    res <- rbind(res, do.call('rbind', df2) )
    setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
  }
  close(pb) # close progressbar
  return(res)
}


#### Filter out permanent water ####
# Returns site ID and if site falls into water (MODE) or not
# Test if this has an effect
waterFilter <- function(){
  library(jsonlite);require(dplyr)
  wa <- fromJSON("../P5_MagnitudeBreakpoints/extracts/PREDICTS_Permanentwater_mode.geojson.json",flatten=T)$features
  wa <- wa %>% dplyr::select(properties.SS,properties.SSBS,properties.mode) %>% 
    dplyr::rename(SS = properties.SS, SSBS = properties.SSBS, WaterMode = properties.mode)
  return(wa)
}