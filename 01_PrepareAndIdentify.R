# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
# Custom packages
library(roquefort)
library(marfunky)
library(xts)
myLog <- function(...) {
  cat(paste0("[Extract] ", Sys.time(), " | ", ..., "\n"))
}
# Exclude single cell studies?
exclude = T

#### Prepare PREDICTS data to be uploaded as G-Fusion table ####
# Load
r <- readRDS("../../Data/diversity-2016-02-03-03-37-46.rds") 

# Convert all SSBS names to correct encoding
r$SSBS <- iconv(r$SSBS,"UTF-8","latin1")
r <- DropInvalidMetricsAndMethods(r)
r <- CorrectSamplingEffort(r)
sites <- SiteMetrics(diversity=r,
                    extra.cols=c("SSB","SSBS","Longitude","Latitude","Sample_start_earliest","Sample_end_latest","Sample_midpoint","Sample_date_resolution",
                                 "Ecoregion","Biome","Country","UN_subregion","Site_name","Order","Family",
                                 "Sampling_method","Study_common_taxon","Max_linear_extent","Coordinates_precision_metres"
                    ),srEstimators = F)


# Proportion of studies with abundance data
length(which(sites$Diversity_metric_type == "Abundance")) / nrow(sites)

#studies = unique(mlist$Reptilia$evi$pa@frame$SS)
#sub <- subset(sites,SS %in% studies)
#subr <- subset(r, SS  %in% studies)

#com <- CompDissim2(subr,"SorVeg",binary = T)
#com2 <- CompDissim2(subr,"Sor",binary = T)
#all_equal(com2$`DL1_2009__Woinarski 2`,com$`DL1_2009__Woinarski 2`) # TRUE

#fit <- glmer(Species_richness ~ Predominant_habitat + (1|SS),data=sub,family=poisson)
#sjPlot::sjp.glmer(fit,type="fe")

# Load specPast
#specPast <- readRDS("Center_FullTimeSeriesSmooth_Final.rds")
# Subset
#specPast <- specPast[which(names(specPast) %in%sub$SSBS)]
#specPast <- lapply(specPast, function(x) return(x$EVI2))
# Get last year
#sub$Earliest <- (sub$Sample_start_earliest - 365)
#test = sub %>% dplyr::select(SSBS,Earliest,Sample_start_earliest)
#for(n in names(specPast)){
#  earliest = test$Earliest[which(test$SSBS==n)]
#  sdate = test$Sample_start_earliest[which(test$SSBS==n)]-1
#  x = specPast[n]
#  specPast[n] <- window(x,start = earliest, end = sdate,extend=F)
#}


# Calculate average and merge back
#x <- as.data.frame(unlist(lapply(specPast, function(x) mean(x,na.rm=T))))
#names(x) <- "meanEVI"
#x$SSBS <- row.names(x)
#sub <- left_join(sub,x)

#fit <- glmer(Species_richness ~ meanEVI + (1|SS/cells),data=sub,family=poisson)
#sjPlot::sjp.glmer(fit,type="fe")

# Alternative with similarity
#sor = com %>% 
 # melt(.) %>% rename(SS = L1, SSBS_x = Var1, SSBS_y = Var2,sor=value)
# Subset from below
#sor <- left_join(per.df,sor)

#sor$meanEVI_x = sub$meanEVI[match(sor$SSBS_x,sub$SSBS)]
#sor$meanEVI_y = sub$meanEVI[match(sor$SSBS_y,sub$SSBS)]
# Difference
#sor$diffMeanEVI = abs(sor$meanEVI_x - sor$meanEVI_y)
#fit <- glmer(value ~ EVIpast + (1|SS),data=sor,family=gaussian)
#sjPlot::sjp.glmer(fit,type="fe")
#infl <- influence.ME::influence(fit,group = "SS")
#cooks.distance(infl)

#res_yearsampling <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"samplingperiod_y0-matrix.rds")) # The year of before start
#res_yearbefore <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"pastperiod_",y,"y-matrix.rds")) # The years before
# Load permutation
#s.metric <- readRDS(f.permut[p])

# Available past data
#focal.period  = ymd("2000.02.18") + (365*(y+1)) # Available past in total
# Add in startdate
#s.metric$startdate_x <- sites$Sample_start_earliest[match(s.metric$SSBS_x,sites$SSBS)]
#s.metric$startdate_y <- sites$Sample_start_earliest[match(s.metric$SSBS_y,sites$SSBS)]
# Calculate distance in days and remove sites with greater dist of 3 months (~90days)
#s.metric <- mutate(s.metric,startd = as.numeric( abs(startdate_x - startdate_y))) %>% 
#  dplyr::filter(startd <= 90) %>% # 3 months (~ 90 days)
#  dplyr::filter(startdate_x > focal.period & startdate_y > focal.period) %>%  # Must have sufficient historical coverage
#  dplyr::select(-startdate_x,-startdate_y,-startd) # remove again

# # Subset and assemble the dataset of the permutations and full matrix
# per.df <- data.frame()
# for(study in unique( s.metric$SS) ){
#     sub <- subset(s.metric,SS == study)
#     if(!(study %in% names(res_yearbefore))) next()
#     # Melt past and sampling cond
#     
#     m <- reshape2::melt( res_yearsampling[[as.character(study)]] ) %>% rename(SSBS_x = Var1,SSBS_y = Var2,EVIpres = value)
#     m2 <- reshape2::melt( res_yearbefore[[as.character(study)]] ) %>% rename(SSBS_x = Var1,SSBS_y = Var2,EVIpast = value)
#     mm <- merge.data.frame(m2,m,by =c("SSBS_x","SSBS_y"),all.x =T) # Merge present with past
#     rm(m,m2)
#     # Merge with permutation to get only the permuteted pairs
#     per.df <- rbind(per.df,
#                     merge.data.frame(sub,mm,by =c("SSBS_x","SSBS_y"),all.x = T)
#     )
#   }
# stopifnot(nrow(s.metric)==nrow(per.df)) # Check. Should be equal to the permutation dataset
# 
# per.df$py <- y # Add year
# # -------------------------------------------- #
# myLog("Running permutation ",p," --- ","Built model sets")
# # Formatting
# per.df$TGrouping <- sites$TGrouping[match(per.df$SS,sites$SS)] # Get taxonomic grouping
# 
# 
# # Subset to reptilia
# per.df <- subset(per.df,TGrouping == "Reptilia")
# fit1 <- glmer(value ~ EVIpres + (1|SS), data=per.df, family = gaussian )
# fit2 <- glmer(value ~ EVIpast + (1|SS), data=per.df, family = gaussian )
# # Use a kenwards roger approximation for testing the p-value
# ( kr <- KRmodcomp.mer(fit2,glmer(value~1 + (1|SS),data=fit2@frame)) )
# 
# sjPlot::sjp.glmer(fit2,type="ma")
# 
# # Calc inv-simpson
# # sites$InvSimpson <- with(r,{
# #   tapply(Measurement, SSBS, function(x) vegan::diversity(x,index="invsimpson"))
# # })
# # 
# # # Calculate exp. shannon index per site
# # sites$Shannon_diversity <- with(r,{
# #   tapply(Measurement, SSBS, function(x) vegan::diversity(x,index="shannon"))
# # })

rm(r)
# Kickout 
sites$startyear <- year(ymd(sites$Sample_start_earliest))
sites$endyear <- year(ymd(sites$Sample_end_latest))
sites <- subset(sites,startyear > 2000) #Kickout everything that started before 2000

# Get distinct 500m grid
library(rgdal)
library(sp)
library(raster)
library("snow")
library("parallel")
sp = subset(sites,select = c("SSBS","Longitude","Latitude"))
# Make spatial file
coordinates(sp) <- ~Longitude+Latitude
proj4string(sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# read in MCD12Q1 for fishnet proxy
ref <- raster("R:/ecocon_d/shared/GIS_LandCover/MODIS_MCD12Q1/MCD12Q1_2001.tif")
# Harmonize projections
if( proj4string(sp) != proj4string(ref) ){
  warning("Projection of input files is not equal. Will try to reproject point file")
  if(is.na(proj4string(ref))) stop("CRS of RasterLayer must be set!")
  sp <- spTransform(sp,CRSobj = CRS(proj4string(ref)))
}

### Make cluster object
beginCluster( detectCores()-2 ) # leave two core for background processes
#extraxt point with df
df <- raster::extract(ref,sp,method="simple",cellnumbers=T,df=T)
endCluster()
R.utils::detachPackage("raster")
# Append to data
sites$cells <- df$cells 
sites$SameCellokay <- FALSE

if(exclude){
  # Now loop through all studies 
  for( study in unique(sites$SS)) {
    print(study)
    sub <- subset(sites,SS==study)
    # Filter those studies out which fall into only one MODIS 500m cell
    check = colSums(with(sub,table(SSBS,cells)))
    if(length(check)>1){
      sites$SameCellokay[which(sites$SS==study & sites$cells %in% names(check))] <- TRUE
    } else {
      print(paste("All sites of",study,"fall within one 500m cell -> Exclude"))
    }
  }
  
  # Then filter out those which only fall in one cell
  sites <- sites %>% dplyr::filter(SameCellokay ==TRUE)
}

# Make an alternative grouping
sites$TGrouping <- as.character(sites$Study_common_taxon)
sites$TGrouping[grep("Hymenoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Insecta",x = sites$TGrouping,ignore.case = T)]<- "Invertebrates"
sites$TGrouping[grep("Chordata",x = sites$TGrouping,ignore.case = T)] <- "Other"
sites$TGrouping[grep("Animalia",x = sites$TGrouping,ignore.case = T)] <- "Other"
sites$TGrouping[grep("Formicidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Scarabaeidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Tracheophyta",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Strigiformes",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Isoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Coleoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Anogeissus",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Poaceae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Colubridae",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Chiroptera",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Ascomycota",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Lepidoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Bryophyta",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Sarcoptiformes",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[which(str_length(sites$TGrouping)==0)] <- "Other"
sites$TGrouping[grep("Bombus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Apidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Arthropoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Drosophilidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Colletes floralis",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Phasianidae",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Lophophorus impejanus",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Gastropoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Araneae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Arachnida",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Clitellata",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Anura",x = sites$TGrouping,ignore.case = T)] <- "Amphibia"
sites$TGrouping[grep("Carabidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Hemiptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Isopoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Collembola",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Agaricomycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Curculionidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Pongo pygmaeus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Squamata",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Culicidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Phyllostomidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Maerua subcordata",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Oryctolagus cuniculus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Arecaceae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Pteropus tonganus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Nymphalidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Diptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Staphylinidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Opiliones",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Orthoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Swietenia macrophylla",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Aenictus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Dorylus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Vespertilionidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Primates",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Panthera pardus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Odocoileus virginianus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Cephalophus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Geometridae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Rodentia",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Magnoliopsida",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Sciomyzidae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Liolaemus",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Dolichopus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Muridae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Soricidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Lumbricidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Lecanoromycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Clethrionomys gapperi",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Passeriformes",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Dipteryx oleifera",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Nematoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Diprotodontia",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Glomeromycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Strabomantidae",x = sites$TGrouping,ignore.case = T)] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2008__Schon 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2008__Schon 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DG1_2012__Ge 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2009__Woinarski 2")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2012__Dominguez 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HP1_2010__Bicknell 1")] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HP1_2010__Bicknell 2")] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MJ1_2009__Lehouck 2")] <- "Sturnidae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2014__Kurz 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2010__Gaigher 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2012__Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2014a_Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2014b_Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SH1_2011__Todd 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SH1_2013__Peri 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SH1_2014__Walker 3")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__Carpenter 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__Carpenter 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__Carpenter 6")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2015__Mumme 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DI1_2012__Muchane 1")] <- "Invertebrates"
sites$TGrouping[grep("Sturnidae",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "AR1_2008__Basset 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2005__Barratt 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2005__Barratt 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "YP1_2012__Sung 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2006__Norton 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VK1_2007__StLaurent 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2013__Burton 3")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2013__Burton 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012a_Carpenter 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2009__Boutin 3")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2009__Boutin 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2009__Boutin 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2012__LeightonGoodall 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2008__Smith 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2008a_Smith 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2006__Smith 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2006__Smith 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2005__Eggleton 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "VB1_2005__Eggleton 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2010__McCarthy 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2009__Christensen 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "TN1_2008__Ngai 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "TN1_2007__Gardner 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2006__UrbinaCardona 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SC1_2005__Richardson 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MJ1_2009__Lehouck 1")] <- "Aves"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2010__Schon 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2010__Schon 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE2_2009__Craig 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "SE1_2012__Lopez 1")] <- "Fungi"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 5")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 4")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 3")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 2")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MH1_2010__CATIE 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MG1_2011__Schon 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "MG1_2008__Buscardo 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "LK1_2010__Endo 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "LK1_2009__Hayward 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "JD1_2004__Alcala 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HZ1_2012__Kutt 1")] <- "Amphibia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HW1_2011__Robinson 1")] <- "Fungi"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "GP1_2009__Vasconcelos 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "HB1_2009__Parry 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "GP1_2007__Kutt 1")] <- "Reptilia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2013__deThoisy 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2008__MacipRios 1")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "CC1_2010__Schon 2")] <- "Invertebrates"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DB1_2010__Garden 1")] <- "Mammalia"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DI1_2008__Noeske 1")] <- "Plantae"
sites$TGrouping[which(sites$TGrouping=="Other" & sites$SS == "DL1_2008__MacipRios 1")] <- "Invertebrates"

#sites$Family[which(sites$SS == "DL1_2008__MacipRios 1")]
stopifnot( length(sites$SS[which(sites$TGrouping=="Other")]) == 0 )

# Check number of studies
rowSums(table(sites$TGrouping,sites$Source_ID)>0)

# Assess if there are multiple Tgrouping within a study
which( rowSums( table(sites$SS,sites$TGrouping)>1 ) >1  ) # NO

## Reclass Max linear extent ##
sites.ori <- sites
d <- sites
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  dplyr::group_by(Study_common_taxon,Sampling_method) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- dplyr::left_join(d,temp) # Join back
# Insert where empty
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Then check again this time with Higher Taxa at the remaining empty values
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,TGrouping) %>% 
  group_by(TGrouping,Sampling_method) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Last try
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  group_by(Study_common_taxon) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
rm(temp)

# How many left
print(paste("How many remaining without MaxLinExtent:", length(which(is.na(d$Max_linear_extent)))))
# What proportion is bigger than 500 meters
length(which(d$Max_linear_extent>500)) / nrow(d)
sites <- d
rm(d)
# Lookup skewed Max_Linear_Extent
unique( sites$Source_ID[which(sites$Max_linear_extent>5000)] )

# Do some sensiblity testing using jackniving 
d <- sites.ori
# Replace 25
d$Max_linear_extent[sample( which(!is.na(d$Max_linear_extent)),size = 0.25 * length( which(!is.na(d$Max_linear_extent)) ) )] <- NA
d$EmptyMLE <- ifelse(is.na(d$Max_linear_extent),1,0)

temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  dplyr::group_by(Study_common_taxon,Sampling_method) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- dplyr::left_join(d,temp) # Join back
# Insert where empty
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Then check again this time with Higher Taxa at the remaining empty values
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,TGrouping) %>% 
  group_by(TGrouping,Sampling_method) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Last try
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  group_by(Study_common_taxon) %>% 
  dplyr::summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
rm(temp)

# Compare new with old estimates
d <- d %>% dplyr::filter(EmptyMLE==1)
s2 <- sites[which(sites$SSBS %in% d$SSBS),] # Subset to those same sites with previous gaps
plot(log10(s2$Max_linear_extent),log10(d$Max_linear_extent))
cor.test( log10(s2$Max_linear_extent),log10(d$Max_linear_extent) )
rm(d)


# New log-transformed columns
sites$logabund<- log10(sites$Total_abundance+1)
sites$logsimp <- log10(sites$Simpson_diversity)


#Save
saveRDS(sites,"sites_center_all.rds")
#saveRDS(sites,"sites_center.rds")

# Export center coordinates as KML
library(rgdal)
sp = sites %>% dplyr::select(SSBS,Longitude,Latitude,startyear,endyear)
coordinates(sp) <- ~Longitude+Latitude
proj4string(sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
writeOGR(sp,"PREDICTS_centerCoordinates_all.kml","PREDICTS_center",driver = "KML",overwrite_layer = T)

# --------------------------------------------------- # 
stop("Preprocessing finished")
# --------------------------------------------------- #


#### Test sampling duration with studies ####
library(lubridate)
library(ggplot2)
library(marfunky)

# Load in
r <- readRDS("../../Data/diversity-2016-02-03-03-37-46.rds") 
r <- DropInvalidMetricsAndMethods(r)
r <- CorrectSamplingEffort(r)
sites <- SiteMetrics(diversity=r,
                     extra.cols=c("SSB","SSBS","Longitude","Latitude","Sample_start_earliest","Sample_end_latest","Sample_midpoint","Sample_date_resolution",
                                  "Ecoregion","Biome","Country","UN_subregion","Site_name","Order","Family",
                                  "Sampling_method","Study_common_taxon","Max_linear_extent","Coordinates_precision_metres"
                     ))
# Get study duration in days
sites$StudyLength <- (sites$Sample_end_latest - sites$Sample_start_earliest) # Study duration

# Load my subset
sites <- readRDS("AllPredictsSites_full.rds") # All
sites <- subset(sites, startyear > 2000)

# ---- #
# Sampling start

# How many studies with differing sampling start of sites ?
x <- sites %>% dplyr::group_by(SS) %>% 
  dplyr::summarise(Nr.Sites = length(unique(SSBS)), # Number of sites
                   Nr.UniqueSamplingStart = length(unique(Sample_start_earliest)), #Overall
                   Nr.UniqueSamplingStart.year = length(unique(year(Sample_start_earliest))), #Year
                   Nr.UniqueSamplingStart.month = length(unique(month(Sample_start_earliest))) #Month
                   ) %>% 
  dplyr::filter(Nr.UniqueSamplingStart > 1) %>% # Remove studies where N = 1
  arrange(desc(Nr.UniqueSamplingStart))
save.xlsx("Assessment_SamplingStartUnique.xlsx",as.data.frame(x))

# -
# What is range of sampling duration (Sampling end date – Sampling start date) of sites within studies ?
s <- sites %>% group_by(SS) %>% 
  dplyr::summarise(Study.dur.min = min(StudyLength),
                   Study.dur.max = max(StudyLength)
  ) %>% 
  dplyr::filter(Study.dur.min != Study.dur.max ) %>% # Kick out those where study duration min is equal to max
  mutate(Abs.diff = abs(Study.dur.max - Study.dur.min) ) %>% 
  arrange(desc(Abs.diff)) %>% 
  mutate(Abs.diff = as.numeric(Abs.diff))
s

g <- ggplot(s,aes(Abs.diff)) + theme_few() +
  geom_histogram() +
  scale_x_continuous(breaks=pretty_breaks(10))  +
  scale_y_continuous(breaks=pretty_breaks(10),expand=c(0,0)) +
  labs(y = "Number studies",x = "Absolute difference of sites sampling duration in days \n (Longest sampling duration - shortest sampling duration)",
       title = paste0("Number of studies = ",nrow(s)))
g
ggsave("Assessment_SamplingDuration.png",plot=g)


#-
# What is the range of delay (min – max) in sampling start dates between sites within a study?

# Create subset based on those studies with known different sampling starts
s <- sites %>% dplyr::group_by(SS) %>% 
  dplyr::summarise(Nr.UniqueSamplingStart = length(unique(Sample_start_earliest))) %>% 
  dplyr::filter(Nr.UniqueSamplingStart > 1) %>% dplyr::select(SS)
s <- left_join(s,sites) %>% dplyr::select(SS,SSB,SSBS,Sample_start_earliest)

# Now loop through and calculate the delay in start date
res <- data.frame()
for(study in unique(s$SS)){
  print(study)
  sub <-  subset(s,SS == study)   
  mat <- matrix(nrow = length(sub$SSBS), ncol = length(sub$SSBS),dimnames = list(sub$SSBS,sub$SSBS))
  for(i in 1:length(sub$SSBS)){
    for(j in 1:length(sub$SSBS)){
      # Calc absolute difference
      mat[i,j] <- abs( sub$Sample_start_earliest[i] - sub$Sample_start_earliest[j]  )
    }
  }
  diag(mat) <- NA
  mat[upper.tri(mat)] <- NA
  res <- rbind(res,
               data.frame(SS= study, min.delay = min(mat,na.rm=T), max.delay = max(mat,na.rm=T),avg.delay = mean(mat,na.rm=T) )
               )
}

# Make a range plot
library(ggplot2)
library(scales)
library(ggExtra)
library(ggthemes)
g1 <- ggplot(res,aes(x = max.delay, y = avg.delay )) + theme_few() +
  geom_point() +
  scale_x_continuous(breaks=pretty_breaks(10)) + scale_y_continuous(breaks=pretty_breaks(10)) + 
  labs(x="Maximal delay \n(in days)",y="Average delay \n(in days)")
g1
gg <- ggMarginal(g1,size=2,type="histogram")
ggsave("Assessment_SamplingstartRange.png",plot=gg)

# Does the Sample_date_resolution (accucary of sampling date entries) can help explain some of those differences?
# Calculate range above
s <- sites %>% dplyr::select(SS,Sample_date_resolution) %>% distinct() %>% 
  right_join(.,res)

g <- ggplot(s,aes(x = Sample_date_resolution, y =  avg.delay)) + theme_few() +
  geom_boxplot() +
  geom_jitter(size = 2,alpha=.5) +
  scale_y_continuous(breaks=pretty_breaks(10))  +
  labs(y="Average delay \n(in days)")
g
ggsave("Assessment_SamplingStartDateResolution.png",plot=g)
