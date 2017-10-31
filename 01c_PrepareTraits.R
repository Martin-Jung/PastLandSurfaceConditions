#### Load packages and classify into taxonomic groups ####
library(data.table)
library(lubridate)
library(stringr)
library(dplyr)
library(roquefort)
library(yarg)

r <- readRDS("../../Data/diversity-2016-02-03-03-37-46.rds")
r <- DropInvalidMetricsAndMethods(r)
r <- CorrectSamplingEffort(r)

# Format to relevant studies
r$startyear <- year(ymd(r$Sample_start_earliest))
r$endyear <- year(ymd(r$Sample_end_latest))
r <- subset(r,startyear >= 2000) #Kickout everything that started before 2000

### Format taxonomic grouping ###
r$TGrouping <- as.character(r$Study_common_taxon)
r$TGrouping[grep("Hymenoptera",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Insecta",x = r$TGrouping,ignore.case = T)]<- "Invertebrates"
r$TGrouping[grep("Chordata",x = r$TGrouping,ignore.case = T)] <- "Other"
r$TGrouping[grep("Animalia",x = r$TGrouping,ignore.case = T)] <- "Other"
r$TGrouping[grep("Formicidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Scarabaeidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Tracheophyta",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Strigiformes",x = r$TGrouping,ignore.case = T)] <- "Aves"
r$TGrouping[grep("Isoptera",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Coleoptera",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Anogeissus",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Poaceae",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Colubridae",x = r$TGrouping,ignore.case = T)] <- "Reptilia"
r$TGrouping[grep("Chiroptera",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Ascomycota",x = r$TGrouping,ignore.case = T)] <- "Fungi"
r$TGrouping[grep("Lepidoptera",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Bryophyta",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Sarcoptiformes",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[which(str_length(r$TGrouping)==0)] <- "Other"
r$TGrouping[grep("Bombus",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Apidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Arthropoda",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Drosophilidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Colletes floralis",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Phasianidae",x = r$TGrouping,ignore.case = T)] <- "Aves"
r$TGrouping[grep("Lophophorus impejanus",x = r$TGrouping,ignore.case = T)] <- "Aves"
r$TGrouping[grep("Gastropoda",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Araneae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Arachnida",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Clitellata",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Anura",x = r$TGrouping,ignore.case = T)] <- "Amphibia"
r$TGrouping[grep("Carabidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Hemiptera",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Isopoda",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Collembola",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Agaricomycetes",x = r$TGrouping,ignore.case = T)] <- "Fungi"
r$TGrouping[grep("Curculionidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Pongo pygmaeus",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Squamata",x = r$TGrouping,ignore.case = T)] <- "Reptilia"
r$TGrouping[grep("Culicidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Phyllostomidae",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Maerua subcordata",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Oryctolagus cuniculus",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Arecaceae",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Pteropus tonganus",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Nymphalidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Diptera",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Staphylinidae",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Opiliones",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Orthoptera",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Swietenia macrophylla",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Aenictus",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Dorylus",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Vespertilionidae",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Primates",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Panthera pardus",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Odocoileus virginianus",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Cephalophus",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Geometridae",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Rodentia",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Magnoliopsida",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Sciomyzidae",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Liolaemus",x = r$TGrouping,ignore.case = T)] <- "Reptilia"
r$TGrouping[grep("Dolichopus",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Muridae",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Soricidae",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Lumbricidae",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Lecanoromycetes",x = r$TGrouping,ignore.case = T)] <- "Fungi"
r$TGrouping[grep("Clethrionomys gapperi",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Passeriformes",x = r$TGrouping,ignore.case = T)] <- "Aves"
r$TGrouping[grep("Dipteryx oleifera",x = r$TGrouping,ignore.case = T)] <- "Plantae"
r$TGrouping[grep("Nematoda",x = r$TGrouping,ignore.case = T)] <- "Invertebrates"
r$TGrouping[grep("Diprotodontia",x = r$TGrouping,ignore.case = T)] <- "Mammalia"
r$TGrouping[grep("Glomeromycetes",x = r$TGrouping,ignore.case = T)] <- "Fungi"
r$TGrouping[grep("Strabomantidae",x = r$TGrouping,ignore.case = T)] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2008__Schon 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2008__Schon 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DG1_2012__Ge 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DL1_2009__Woinarski 2")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DL1_2012__Dominguez 1")] <- "Plantae"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "HP1_2010__Bicknell 1")] <- "Aves"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "HP1_2010__Bicknell 2")] <- "Aves"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MJ1_2009__Lehouck 2")] <- "Sturnidae"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SC1_2014__Kurz 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SE2_2010__Gaigher 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SE2_2012__Craig 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SE2_2014a_Craig 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SE2_2014b_Craig 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SH1_2011__Todd 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SH1_2013__Peri 1")] <- "Plantae"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SH1_2014__Walker 3")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2012__Carpenter 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2012__Carpenter 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2012__Carpenter 6")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2015__Mumme 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DI1_2012__Muchane 1")] <- "Invertebrates"
r$TGrouping[grep("Sturnidae",x = r$TGrouping,ignore.case = T)] <- "Aves"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "AR1_2008__Basset 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2005__Barratt 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2005__Barratt 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "YP1_2012__Sung 1")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2006__Norton 1")] <- "Plantae"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VK1_2007__StLaurent 1")] <- "Mammalia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2013__Burton 3")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2013__Burton 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2012a_Carpenter 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2009__Boutin 3")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2009__Boutin 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2009__Boutin 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2012__LeightonGoodall 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2008__Smith 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2008a_Smith 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2006__Smith 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2006__Smith 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2005__Eggleton 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "VB1_2005__Eggleton 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SE2_2010__McCarthy 1")] <- "Mammalia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SC1_2009__Christensen 1")] <- "Plantae"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "TN1_2008__Ngai 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "TN1_2007__Gardner 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SC1_2006__UrbinaCardona 1")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SC1_2005__Richardson 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MJ1_2009__Lehouck 1")] <- "Aves"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2010__Schon 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2010__Schon 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SE2_2009__Craig 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "SE1_2012__Lopez 1")] <- "Fungi"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MH1_2010__CATIE 5")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MH1_2010__CATIE 4")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MH1_2010__CATIE 3")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MH1_2010__CATIE 2")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MH1_2010__CATIE 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MG1_2011__Schon 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "MG1_2008__Buscardo 1")] <- "Plantae"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "LK1_2010__Endo 1")] <- "Mammalia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "LK1_2009__Hayward 1")] <- "Mammalia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "JD1_2004__Alcala 1")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "HZ1_2012__Kutt 1")] <- "Amphibia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "HW1_2011__Robinson 1")] <- "Fungi"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "GP1_2009__Vasconcelos 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "HB1_2009__Parry 1")] <- "Mammalia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "GP1_2007__Kutt 1")] <- "Reptilia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DL1_2013__deThoisy 1")] <- "Mammalia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DL1_2008__MacipRios 1")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "CC1_2010__Schon 2")] <- "Invertebrates"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DB1_2010__Garden 1")] <- "Mammalia"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DI1_2008__Noeske 1")] <- "Plantae"
r$TGrouping[which(r$TGrouping=="Other" & r$SS == "DL1_2008__MacipRios 1")] <- "Invertebrates"

#r$Family[which(r$SS == "DL1_2008__MacipRios 1")]
stopifnot( length(r$SS[which(r$TGrouping=="Other")]) == 0 )
assertthat::assert_that(!("Other"%in% unique(r$TGrouping) ))

# Now get only studies with full dataset included
s <- readRDS("AllPredictsSites_full.rds") %>% dplyr::select(SS,SSB,SSBS,TGrouping)
r <- left_join(s,r)
rm(s)
# --- Done --- #
stop("Stop loading here")
#### Load and aggregate bodymass from Amniote database ####
# MYHRVOLD et al. 2011
# db.amni <- read.csv("../../Data/Traits/Amniote_Database_Aug_2015.csv",header=T,na.strings = "-999")
# db.amni$Best_guess_binomial <- paste(db.amni$genus,db.amni$species) # Assemble species names
# db.amni <- db.amni %>% dplyr::select(class,family,genus,Best_guess_binomial,adult_body_mass_g) %>% 
#   rename(TGrouping = class,Family = family, Genus = genus)
# # Subset PREDICTS sites to covered taxonomic classes and distinct species names
# r.sub <- r %>% dplyr::filter(TGrouping %in% c("Aves","Mammalia","Reptilia")) %>% 
#   dplyr::select(TGrouping,Family,Genus,Best_guess_binomial) %>% distinct() %>% dplyr::filter(str_length(Best_guess_binomial)>1)
# 
# # First match with binomial species name
# r.sub <- merge.data.frame(r.sub, db.amni[,c("Best_guess_binomial","adult_body_mass_g")],by = c("Best_guess_binomial"),all.x = T)
# length(which(is.na(r.sub$adult_body_mass_g))) / nrow(r.sub) # How many missing entries
# 
# # Interpolate genus
# temp <-  r.sub %>% dplyr::filter(!is.na(Genus) | str_length(Best_guess_binomial)!=0 ) %>% 
#   dplyr::group_by(TGrouping,Genus) %>% 
#   dplyr::summarise(bm_g_avg = mean(adult_body_mass_g,na.rm=T))
# r.sub <- dplyr::left_join(r.sub,temp) # Join back
# 
# # interpolate among family and grouping to fill missing data
# temp <-  r.sub %>% dplyr::filter(!is.na(Family)) %>% 
#   dplyr::group_by(TGrouping,Family) %>% 
#   dplyr::summarise(bm_f_avg = mean(adult_body_mass_g,na.rm=T))
# r.sub <- dplyr::left_join(r.sub,temp) # Join back
# 
# # Do the same per taxonomic group only
# temp <-  r.sub %>%
#   dplyr::group_by(TGrouping) %>% 
#   dplyr::summarise(bm_tg_avg = mean(adult_body_mass_g,na.rm=T))
# r.sub <- dplyr::left_join(r.sub,temp) # Join back
# # Insert stepwise where empty
# r.sub$adult_body_mass_g <- ifelse(is.na(r.sub$adult_body_mass_g),r.sub$bm_g_avg,r.sub$adult_body_mass_g) # Per Genus
# r.sub$adult_body_mass_g <- ifelse(is.na(r.sub$adult_body_mass_g),r.sub$bm_f_avg,r.sub$adult_body_mass_g) # Per family
# r.sub$adult_body_mass_g <- ifelse(is.na(r.sub$adult_body_mass_g),r.sub$bm_tg_avg,r.sub$adult_body_mass_g) # Per TGroup
# r.sub$bm_g_avg <- NULL;r.sub$bm_f_avg <- NULL;r.sub$bm_tg_avg <- NULL # Kickout the previously calc. median
# # What proportion is missing ?
# length(which(is.na(r.sub$adult_body_mass_g))) / nrow(r.sub)
# 
# # Calculate log10 biomass
# r.sub <- mutate(r.sub,log10_bodymass =  log10(adult_body_mass_g+1) )
# 
# # Now on the global dataset of bodymass calculate 1/3 quantile 
# TG.quantiles <- db.amni %>% group_by(TGrouping) %>% 
#   summarise(lq.bm33 = quantile(adult_body_mass_g,probs=.33,na.rm = T),
#             mq.bm66 = quantile(adult_body_mass_g,probs=.66,na.rm = T),
#             hq.bm1 = quantile(adult_body_mass_g,probs=1,na.rm = T))
# # Merge with the filled bodymass estimates
# r.sub <- left_join(r.sub,TG.quantiles)
# # Bodymass grouping
# r.sub$BMgroup <- ifelse(r.sub$adult_body_mass_g <= r.sub$lq.bm33,"small",
#                         ifelse(r.sub$adult_body_mass_g > r.sub$lq.bm33 & r.sub$adult_body_mass_g <= r.sub$mq.bm66, "middle",
#                                ifelse(r.sub$adult_body_mass_g > r.sub$mq.bm66,"large","")
#                                )
#                         )
# # Remove quantiles
# r.sub$lq.bm33 <- NULL; r.sub$mq.bm66 <- NULL; r.sub$hq.bm1 <- NULL
# 
# # Merge back to full vertebrate
# r2 <- r %>% dplyr::filter(TGrouping %in% c("Aves","Mammalia","Reptilia")) %>% 
#   left_join(.,r.sub,by=c("TGrouping","Family","Genus","Best_guess_binomial")) %>% 
#   # Filter out those with no species guess
#   dplyr::filter(str_length(Best_guess_binomial) >0) %>% 
#   dplyr::select(SS,SSBS,TGrouping,Best_guess_binomial,adult_body_mass_g,BMgroup,Measurement) %>% 
#   dplyr::filter(Measurement > 0) %>% dplyr::select( -Measurement)# Only species presence
# 
# # NOWWW... Classify studies according to their predominant species bodymass group
# r2_SS <- r2 %>% group_by(SS) %>% 
#   summarise(BMgroup_SS = names(which.max( table(BMgroup) ))
#             )
# # Merge back
# r2 <- left_join(r2,r2_SS)
# # The same for sites (is a certain size class dominant in a class ?)
# ag <- aggregate(list(BMgroup_SSBS = r2$BMgroup),list(SSBS = r2$SSBS),function(x) names(which.max( table(x) )))
# # Merge back
# r2 <- left_join(r2,ag)
# 
# # Filter out empty species (does nothing as all species are filled)
# r2 <- r2 %>% dplyr::select(-Best_guess_binomial,-adult_body_mass_g,-BMgroup) %>% distinct()
# 
# saveRDS(r2,file = "AmnioteLBM_PREDICTS_raw.rds")
# 
# #### TRY - Plants #####
# # Get Plant species in PREDICTS
# r.sub2 <- r %>% dplyr::filter(TGrouping %in% c("Plantae")) %>% 
#   dplyr::select(TGrouping,Family,Genus,Best_guess_binomial) %>% distinct() 
# # First match with genus (species basically not existing consistenly)
# # Need to construct the genus name wherever missing 
# r.sub2$Best_guess_binomial <- as.character(r.sub2$Best_guess_binomial)
# r.sub2$Genus2 <- do.call(rbind, lapply(str_split(r.sub2$Best_guess_binomial," "), function(x) return(x[1])))
# r.sub2$Genus <- ifelse(str_length(r.sub2$Genus)==0,r.sub2$Genus2,as.character(r.sub2$Genus))
# r.sub2$Genus2 <- NULL
# 
# # Try and Load TRY (nice pun)
# trydb <- fread("../../Data/Traits/PlantsTRY_05042016123306/1886.txt",showProgress = T)
# unique(trydb$TraitName)
# # Dispersal syndromes
# #trydb %>% dplyr::filter(TraitName == "Dispersal syndrome") %>% dplyr::select(OrigValueStr) %>% distinct()
# # Plant lifeform
# #trydb %>% dplyr::filter(TraitName == "Plant life form (Raunkiaer life form)") %>% dplyr::select(OrigValueStr) %>% distinct()
# 
# # Height only
# trydb2 <- trydb %>% dplyr::filter(TraitName == "Plant height") %>% # Only height
#          dplyr::select(AccSpeciesName,UnitName,StdValue) %>% # Relevant columns
#         dplyr::rename(Plant.height = StdValue)
# # Using the standardized values
# # Thus already converted into numeric and error corrected (not the case with original values)
# 
# # Now convert and filter and aggregate per species
# plh <- trydb2 %>% 
#   dplyr::group_by(AccSpeciesName) %>%  # First aggregate by species name
#   dplyr::summarise(Plant_height_avg = mean(Plant.height,na.rm=T)) %>% 
#   ungroup() %>% distinct() %>% 
#   rename(Best_guess_binomial = AccSpeciesName) %>% 
#   # Filter out missing entries
#   subset(., complete.cases(.))
# 
# # Assemble artifical genus name
# plh <- plh %>% mutate(Genus = do.call(rbind, lapply(str_split(plh$Best_guess_binomial," "), function(x) return(x[1])))) %>% 
#   mutate(Genus = str_to_lower(Genus)) %>% mutate(Genus = Hmisc::capitalize(Genus))
# 
# # Check how many in PREDICTS can be matched to this Genus name
# plh$Level <- NA
# plh$Level[union(which(plh$Genus %in% r.sub2$Genus),which(is.na(plh$Level)))] <- "Genus"
# # Do an iniational check which species occur in PREDICTS
# plh$Level[which(plh$Best_guess_binomial %in% r.sub2$Best_guess_binomial)] <- "Species"
# 1 - (length(which(!is.na(plh$Level))) / nrow(plh)) # How many of try genus are missing in PREDICTS ?
# # NONE!
# #How many are species vs genus reso?
# length(which(plh$Level=="Genus"));length(which(plh$Level=="Species"))
# length(which(plh$Level=="Species")) / nrow(plh) # Only 14.83 % can be matched at species level
# # --------------------------------------------- #
# # Calculate average plant height for Genus of all species present within the genus
# plh2 <- plh %>% 
#   dplyr::select(Genus,Plant_height_avg) %>% 
#   group_by(Genus) %>%  # Calculate median plant height by genus
#   summarise(Plant_height_avg = mean(Plant_height_avg,na.rm=T)) %>% 
#   rename(PH_genus = Plant_height_avg)
# # Left_join with plh
# plh <- left_join(plh,plh2)
# 
# # Generall TH-quantiles
# TH.quantiles <- plh %>%  
#   summarise(lq.bm33 = quantile(Plant_height_avg,probs=.33,na.rm = T),
#             mq.bm66 = quantile(Plant_height_avg,probs=.66,na.rm = T),
#             hq.bm1 = quantile(Plant_height_avg,probs=1,na.rm = T))
# 
# # Now join with PREDICTS data
# r.subp <- r %>% dplyr::filter(TGrouping %in% c("Plantae")) %>% 
#   dplyr::select(SS,SSBS,Genus,Best_guess_binomial,Measurement) %>% 
#   # Filter out absent species
#   dplyr::filter(Measurement > 0) %>% dplyr::select(-Measurement)
# # Construct Genus wherever missing
# r.subp$Best_guess_binomial <- as.character(r.subp$Best_guess_binomial)
# r.subp$Genus2 <- do.call(rbind, lapply(str_split(r.subp$Best_guess_binomial," "), function(x) return(x[1])))
# r.subp$Genus <- ifelse(str_length(r.subp$Genus)==0,r.subp$Genus2,as.character(r.subp$Genus))
# r.subp$Genus2 <- NULL
# # Then Filter out the few where both species and genus is NA (based on assumption that average per family is to inaccurate)
# r.subp <- r.subp[which(r.subp$Genus!=""),] # Now everything has at least a genus or species entry
# 
# # Merge quantiles with the PREDICTS
# r.subp <- cbind(r.subp,TH.quantiles)
# 
# # Insert tree height
# r.subp$PlantHeight <- NA
# # Then overwrite with species depending on availabilty (match)
# r.subp$PlantHeight <- plh$Plant_height_avg[match(r.subp$Best_guess_binomial,plh$Best_guess_binomial)]
# # Fill the missing with genus averaged height
# r.subp2 <- r.subp[which(is.na(r.subp$PlantHeight)),]
# # Fill missing
# r.subp2$PlantHeight <- plh$Plant_height_avg[match(r.subp2$Genus,plh$Genus)]
# # Combine both to get a height estimate for each site
# r.subp <- rbind( subset(r.subp,complete.cases(r.subp)),subset(r.subp2,complete.cases(r.subp2) ))
# r.subp2 <- NULL
# # All entries there...
# summary( r.subp$PlantHeight )
# 
# # Plant height grouping
# r.subp$THgroup <- ifelse(r.subp$PlantHeight <= r.subp$lq.bm33,"small",
#                         ifelse(r.subp$PlantHeight > r.subp$lq.bm33 & r.subp$PlantHeight <= r.subp$mq.bm66, "middle",
#                                ifelse(r.subp$PlantHeight > r.subp$mq.bm66,"large","")
#                         )
# )
# 
# r.subp$lq.bm33 <- NULL; r.subp$mq.bm66 <- NULL; r.subp$hq.bm1 <- NULL
# 
# # NOWWW... Classify studies according to their predominant species bodymass group
# ag <- aggregate(list(BMgroup_SS = r.subp$THgroup),list(SS = r.subp$SS),function(x) names(which.max( table(x) )))
# 
# # Merge back
# r.subp <- left_join(r.subp,ag)
# # The same for sites (is a certain size class dominant in a class ?)
# ag <- aggregate(list(BMgroup_SSBS = r.subp$THgroup),list(SSBS = r.subp$SSBS),function(x) names(which.max( table(x) )))
# # Merge back
# r.subp <- left_join(r.subp,ag)
# 
# # Filter out empty species (does nothing as all species are filled)
# r.subp <- r.subp %>% dplyr::select(-Best_guess_binomial,-Genus,-PlantHeight,-THgroup) %>% distinct() %>% 
#   mutate( TGrouping = "Plantae")
# 
# # Save
# saveRDS(r.subp,file = "TryDBHeight.rds")
# 
# # Quick summ
# Hmisc::describe(r.subp$BMgroup_SS)
# psych::describeBy(r.subp$BMgroup_SSBS,r.subp$TGrouping)
# 
# 
# ### ---- ###
# # Plant lifespan
# r.sub2 <- r %>% dplyr::filter(TGrouping %in% c("Plantae")) %>% 
#   # For lifespan filter per only scientific name and longer than 0
#   dplyr::filter(Resolution_entered == "Scientific",str_length(Best_guess_binomial)>0) %>% 
#   dplyr::select(TGrouping,Family,Genus,Best_guess_binomial) %>% distinct() 
# # First match with genus (species basically not existing consistenly)
# # Need to construct the genus name wherever missing 
# r.sub2$Best_guess_binomial <- as.character(r.sub2$Best_guess_binomial)
# r.sub2$Genus2 <- do.call(rbind, lapply(str_split(r.sub2$Best_guess_binomial," "), function(x) return(x[1])))
# r.sub2$Genus <- ifelse(str_length(r.sub2$Genus)==0,r.sub2$Genus2,as.character(r.sub2$Genus))
# r.sub2$Genus2 <- NULL
# 
# trydb.lo <- trydb %>% dplyr::filter(TraitName == "Plant lifespan (longevity)") %>%
#   dplyr::select(AccSpeciesName,UnitName,StdValue,OrigValueStr) %>% # Relevant columns
#   dplyr::rename(Plant.lifespan = OrigValueStr)
# # Assemble Genus
# trydb.lo$Genus <- as.character(do.call(rbind, lapply(str_split(trydb.lo$AccSpeciesName," "), function(x) return(x[1]))))
# trydb.lo <- trydb.lo %>% 
#   mutate(Genus = str_to_lower(Genus)) %>% mutate(Genus = Hmisc::capitalize(Genus)) %>% 
#   dplyr::select(Genus,AccSpeciesName,StdValue,Plant.lifespan) %>% distinct() %>% 
#   rename(Best_guess_binomial = AccSpeciesName)
# 
# # Filter that
# trydb.lo$Level <- NA
# trydb.lo$Level[match(trydb.lo$Genus,r.sub2$Genus)] <- "Genus"
# # Species]
# trydb.lo$Level[match(trydb.lo$Best_guess_binomial,r.sub2$Best_guess_binomial)] <- "Species"
# # Do an iniational check which species occur in PREDICTS
# 1 - (length(which(!is.na(trydb.lo$Level))) / nrow(trydb.lo)) # How many of plants do I have
# trydb.lo <- subset(trydb.lo,!is.na(Level))
# 
# # Reclassify into new coarse categories 
# # First on standardized age column
# trydb.lo$ShortLong <- NA
# trydb.lo$ShortLong <- ifelse(trydb.lo$StdValue>5,"Long-lived","Short-lived")
# t1 <- subset(trydb.lo,is.na(ShortLong))
# t2 <- subset(trydb.lo,!is.na(ShortLong))
# # Then manually on original str
# unique(t1$Plant.lifespan[which(is.na(t1$ShortLong))])
# t1$ShortLong[which(t1$Plant.lifespan=="Annual")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="Short")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="Annual, Biennial")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always annual")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always pluriennial-pollakanthic")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always annual, always biennial, always pluriennial-hapaxanthic, always pluriennial-pollakanthic")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always annual, always biennial")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always annual, always biennial, always pluriennial-pollakanthic")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always annual, always pluriennial-pollakanthic")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always annual, sometimes biennial")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always biennial")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="always biennial, sometimes pluriennial-hapaxanthic")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="sometimes biennial, always pluriennial-pollakanthic")] <- "Short-lived"
# t1$ShortLong[which(t1$Plant.lifespan==">5")] <- "Short-lived"
# 
# # Long-lived
# t1$ShortLong[which(as.numeric( gsub("\\D","",t1$Plant.lifespan) )>5)] <- "Long-lived"
# t1$ShortLong[which(t1$Plant.lifespan=="Long")] <- "Long-lived" 
# t1$ShortLong[which(t1$Plant.lifespan=="Moderate")] <- "Long-lived" # Assumption
# t1$ShortLong[which(t1$Plant.lifespan=="Perennial")] <- "Long-lived" 
# t1$ShortLong[which(t1$Plant.lifespan=="always pluriennial-hapaxanthic")] <- "Long-lived" 
# t1$ShortLong[which(t1$Plant.lifespan=="always pluriennial-hapaxanthic, always pluriennial-pollakanthic")] <- "Long-lived" 
# 
# # How many missing still
# length(which(is.na(t1$ShortLong)==T)) / nrow(t1)
# 
# tfull <- rbind(t1,t2)
# 
# 
# # --------------------------------------------- #
# # Calculate average plant height for Genus of all species present within the genus
# 
# # Then Filter out the few where both species and genus is NA (based on assumption that average per family is to inaccurate)
# r.sub2 <- r.sub2[which(r.sub2$Genus!=""),] # Now everything has at least a genus or species entry
# 
# 
# ag <- aggregate(list(BMgroup_SS = r.subp$THgroup),list(SS = r.subp$SS),function(x) names(which.max( table(x) )))
# 
# # Merge back
# r.subp <- left_join(r.subp,ag)
# # The same for sites (is a certain size class dominant in a class ?)
# ag <- aggregate(list(BMgroup_SSBS = r.subp$THgroup),list(SSBS = r.subp$SSBS),function(x) names(which.max( table(x) )))
# # Merge back
# r.subp <- left_join(r.subp,ag)
# 
# 
# 
# #### Proportion of carnivores - EltonianTraits ####
# # Schippers et al 2016
# # The abundance of carnivorous species (i.e., species with a diet consisting of at least 60% meat, fish, and/or carrion)
# # relative to the total abundance of all species (%).
# tr.elton.bird <- fread("../../Data/Traits/Eltonian_Foraging/BirdFuncDat.txt",header=T)
# names(tr.elton.bird) <- str_replace_all(names(tr.elton.bird),"-",".")
# # Calculate sum of meat-based nutrition
# tr.elton.bird <- tr.elton.bird %>% mutate(MeatDiet = Diet.Inv + Diet.Vend + Diet.Vect + Diet.Vfish  + Diet.Vunk + Diet.Scav)
# tr.elton.bird$Carnivore <- ifelse(tr.elton.bird$MeatDiet >= 60,"Carnivore","Non_Carnivore")
# 
# tr.elton.mamm <- fread("../../Data/Traits/Eltonian_Foraging/MamFuncDat.txt",header=T)
# names(tr.elton.mamm) <- str_replace_all(names(tr.elton.mamm),"-",".")
# tr.elton.mamm <- tr.elton.mamm %>% mutate(MeatDiet = Diet.Inv + Diet.Vend + Diet.Vect + Diet.Vfish  + Diet.Vunk + Diet.Scav)
# tr.elton.mamm$Carnivore <- ifelse(tr.elton.mamm$MeatDiet >= 60,"Carnivore","Non_Carnivore")
# 
# # Merge both
# a = tr.elton.bird %>% dplyr::select(SpecID,Scientific,Carnivore) %>% mutate(TGrouping = "Aves")
# b = tr.elton.mamm %>% dplyr::select(MSW3_ID ,Scientific,Carnivore) %>% mutate(TGrouping = "Mammalia") %>% rename(SpecID = MSW3_ID)
# tr.elton <- rbind(a,b) %>% rename(Best_guess_binomial = Scientific)
# 
# # Now merge with PREDICTS mammal and bird species names
# r.sub <- r %>% dplyr::filter(TGrouping %in% c("Aves","Mammalia")) %>% 
#   dplyr::select(TGrouping,Best_guess_binomial) %>% distinct() %>% dplyr::filter(str_length(Best_guess_binomial)>1)
# 
# # Merge with species names
# r.sub <- merge.data.frame(r.sub, tr.elton,by = c("TGrouping","Best_guess_binomial"),all.x = T)
# length(which(is.na(r.sub$Carnivore))) / nrow(r.sub) # How many missing entries
# 
# ## --- ##
# # Merge back with r and calculate proportion relative to total abundance
# r.full <- r %>% dplyr::filter(TGrouping %in% c("Aves","Mammalia")) %>% 
#   left_join(.,r.sub,by=c("TGrouping","Best_guess_binomial"))
# 
# sites.ca <- r.full %>% dplyr::filter(Diversity_metric == "abundance") %>% group_by(SS,SSBS) %>% 
#   summarise(AbundanceCarnivores = sum(Measurement[which(Carnivore == "Carnivore")]))
# 
# saveRDS(sites.ca,"Eltonian_Carnivores.rds")
# 
# #### Prepare my trait files for PREDICTS match ####
# # Trait data dir
EltonianTraitsBirds <- function(path, sep = ",",...) {
  raw.data <- read.csv(path, sep = sep, ...)

  names(raw.data) <- str_replace_all(names(raw.data),"-",".")
  # Assemble trophic level
  yarg:::.Log("Loading Elton bird trophic levels")
  raw.data$Trophic_level <- NA
  raw.data$Trophic_level[which(raw.data$Diet.5Cat == "Omnivore")] <- "Omnivore"
  raw.data$Trophic_level[which(raw.data$Diet.5Cat == "Invertebrate")] <- "Carnivore"
  raw.data$Trophic_level[which(raw.data$Diet.5Cat == "PlantSeed")] <- "Herbivore"
  raw.data$Trophic_level[which(raw.data$Diet.5Cat == "VertFishScav")] <- "Carnivore"
  # Placing fruit nectar eaters under herbivores
  raw.data$Trophic_level[which(raw.data$Diet.5Cat == "FruiNect")] <- "Herbivore"
  raw.data$Trophic_level[which(raw.data$Diet.Fruit < 50)] <- "Omnivore" # If fruit lower than 50, omnivore
  # Subset Trophic levels
  raw.data <- raw.data %>% dplyr::select(Scientific,Trophic_level,Diet.Certainty) %>%
    rename(TrophicLevel = Trophic_level, Taxon = Scientific, Rank = Diet.Certainty) %>%
    mutate(Kingdom = "Animalia", Rank = as.character(Rank))
  raw.data$Rank[which(raw.data$Rank %in% c("A","B"))] <- "Species"
  raw.data$Rank[which(raw.data$Rank %in% c("C","D1"))] <- "Genus"
  # Wrong entrance of uncertainty
  raw.data$Rank[which(raw.data$Diet.Certainty==" Ref_19")][1] <- "Family"
  raw.data$Rank[which(raw.data$Diet.Certainty==" Ref_19")][2] <- "Genus"
  raw.data$Rank[which(raw.data$Diet.Certainty==" Ref_24")] <- "Family"
  raw.data$Rank[which(raw.data$Diet.Certainty==" Ref_22")] <- "Family"
  raw.data$Rank[which(raw.data$Rank %in% c("D2"))] <- "Family"
  # Assemble Genus
  raw.data$Genus <- as.character(do.call(rbind, lapply(str_split(raw.data$Taxon," "), function(x) return(x[1]))))
  raw.data <- raw.data %>% mutate(Genus = str_to_lower(Genus)) %>% mutate(Genus = Hmisc::capitalize(Genus))

  levels <- yarg:::.TrophicLevels()
  raw.data$TrophicLevel <- EnsureFactorLevels(raw.data$TrophicLevel,
                                              levels = levels)

  # The same for bodymass
  cat("\n")
  yarg:::.Log("Loading Elton bird bodymass data")
  eltonbm <- read.csv(path, sep = sep, ...)
  names(eltonbm) <- str_replace_all(names(eltonbm),"-",".")

  trb <- eltonbm  %>% dplyr::select(Scientific,BodyMass.SpecLevel,BodyMass.Value) %>%
    rename(Taxon = Scientific,Rank = BodyMass.SpecLevel,Value = BodyMass.Value) %>%
    mutate(Kingdom = "Animalia",Rank = as.character(Rank))
  trb$Rank[which(trb$Rank==1)] <- "Species"
  trb$Rank[which(trb$Rank==0)] <- "Genus"
  trb$Rank[which(trb$Rank=="GenAvg")] <- "Genus"
  trb$Rank[which(trb$Rank=="Dunning08")] <- "Species"
  trb$Rank[which(trb$Rank=="PrimScale")] <- "Species"

  # Set rank to species as trait data was imputed for species
  return(list(
    Adult_wet_mass = TraitDataset(reference = "EltonianBird",
                                            trait = "Adult wet mass", kingdom = trb$Kingdom, unit = "log10 (g)",
                                            taxa = trb$Taxon, ranks = "Species", values = log10(trb$Value)),
    Trophic_level = TraitDataset(reference = "EltonianBird",
                                           trait = "Trophic level", kingdom =raw.data$Kingdom, unit = "Category",
                                           taxa = raw.data$Taxon, ranks = "Species",
                                           values = raw.data$TrophicLevel)
      )
    )
}
EltonianTraitsMammals <- function(path, sep = ",",...) {
  raw.data <- read.csv(path, sep = sep, ...)

  names(raw.data) <- str_replace_all(names(raw.data),"-",".")
  # Assemble trophic level
  yarg:::.Log("Loading Elton mammal trophic levels")
  raw.data$Trophic_level <- NA
  raw.data$Trophic_level[which(raw.data$Diet.Inv > 50)] <- "Carnivore"
  raw.data$Trophic_level[which(raw.data$Diet.Vend > 50)] <- "Carnivore"
  raw.data$Trophic_level[which(raw.data$Diet.Vfish > 50)] <- "Carnivore"
  raw.data$Trophic_level[which(raw.data$Diet.Vunk > 50)] <- "Carnivore"
  raw.data$Trophic_level[which(raw.data$Diet.Fruit >50)] <- "Herbivore"
  raw.data$Trophic_level[which(raw.data$Diet.Nect >50)] <- "Herbivore"
  raw.data$Trophic_level[which(raw.data$Diet.Seed >50)] <- "Herbivore"
  raw.data$Trophic_level[which(raw.data$Diet.Nect >50)] <- "Herbivore"
  raw.data$Trophic_level[which(raw.data$Diet.PlantO >50)] <- "Herbivore"
  # The rest are supposedly omnivores
  raw.data$Trophic_level[which(is.na(raw.data$Trophic_level))] <- "Omnivore"
  # Subset Trophic levels
  raw.data <- raw.data %>% dplyr::select(Scientific,Trophic_level,Diet.Certainty) %>%
    rename(TrophicLevel = Trophic_level, Taxon = Scientific, Rank = Diet.Certainty) %>%
    mutate(Kingdom = "Animalia", Rank = as.character(Rank))
  raw.data$Rank[which(raw.data$Rank %in% c("ABC"))] <- "Species"
  raw.data$Rank[which(raw.data$Rank %in% c("D1"))] <- "Genus"
  raw.data$Rank[which(raw.data$Rank %in% c("D2"))] <- "Family"
  # Assemble Genus
  raw.data$Genus <- as.character(do.call(rbind, lapply(str_split(raw.data$Taxon," "), function(x) return(x[1]))))
  raw.data <- raw.data %>% mutate(Genus = str_to_lower(Genus)) %>% mutate(Genus = Hmisc::capitalize(Genus))

  levels <- yarg:::.TrophicLevels()
  raw.data$TrophicLevel <- EnsureFactorLevels(raw.data$TrophicLevel,
                                              levels = levels)

  # The same for bodymass
  cat("\n")
  yarg:::.Log("Loading Elton mammal bodymass data")
  eltonbm <- read.csv(path, sep = sep, ...)
  names(eltonbm) <- str_replace_all(names(eltonbm),"-",".")

  trb <- eltonbm  %>% dplyr::select(Scientific,BodyMass.SpecLevel,BodyMass.Value) %>%
    rename(Taxon = Scientific,Rank = BodyMass.SpecLevel,Value = BodyMass.Value) %>%
    mutate(Kingdom = "Animalia",Rank = as.character(Rank))
  trb$Rank[which(trb$Rank==1)] <- "Species"
  trb$Rank[which(trb$Rank==0)] <- "Genus"
  trb$Rank[which(trb$Rank==2)] <- "Family" # Assuming phylogenetically imputed values reflect family level cogenics

  # Rank to species as data were imputed for species
  return(list(
    Adult_wet_mass = TraitDataset(reference = "EltonianMammal",
                                  trait = "Adult wet mass", kingdom = trb$Kingdom, unit = "log10 (g)",
                                  taxa = trb$Taxon, ranks = "Species", values = log10(trb$Value)),
    Trophic_level = TraitDataset(reference = "EltonianMammal",
                                 trait = "Trophic level", kingdom =raw.data$Kingdom, unit = "Category",
                                 taxa = raw.data$Taxon, ranks = "Species",
                                 values = raw.data$TrophicLevel)
  )
  )
}
TRY2 <- function(path, sep = "\t", fileEncoding = "Latin-1", ...) {
  yarg:::.Log("Reading TRY raw data\n")
  #raw.data <- read.csv(path, sep = sep, fileEncoding = fileEncoding,...)
  raw.data <- data.table::fread(path, sep = sep, encoding = fileEncoding, ...)

  traits <- c("Plant height")
  raw.data <- raw.data[raw.data$TraitName %in% traits, ]
  bad <- duplicated(raw.data)
  bad <- bad | raw.data$OrigObsDataID %in% raw.data$ObsDataID
  bad <- bad | is.na(raw.data$StdValue) | is.infinite(raw.data$StdValue) |
    0 == raw.data$StdValue
  if (any(bad)) {
    yarg:::.Log("Removing", sum(bad), "duplicated, NA, Inf or 0 TRY measurements\n")
    raw.data <- raw.data[!bad, ]
  }
  raw.data <- dplyr::select(raw.data,TraitName,StdValue,AccSpeciesName)
  yarg:::.Log("Dropping levels\n")
  raw.data <- droplevels(raw.data)
  raw.data$Log10StdValue <- log10(raw.data$StdValue)
  stopifnot(!any(is.na(raw.data$Log10StdValue) | is.infinite(raw.data$Log10StdValue)))
  # Calculate mean trait value over all species names
  F <- function(trait.name) {
    trait <- droplevels(raw.data[trait.name == raw.data$TraitName,
                                 ])
    return(tapply(trait$Log10StdValue, trait$AccSpeciesName,
                  mean))
  }
  yarg:::.Log("Computing species-level means\n")
  #seed.mass <- F("Seed mass")
  veg.height <- F("Plant height")
  #gen.height <- F("Plant height generative")
  return(list(
    # Seed_mass = TraitDataset(reference = "TRY", trait = "Seed mass",
    #                                    kingdom = "Plantae", unit = "log10 (g)", taxa = names(seed.mass),
    #                                    ranks = "Species", values = as.vector(seed.mass)), Vegetative_height = TraitDataset(reference = "TRY",
    #                                                                                                                        trait = "Vegetative height", kingdom = "Plantae", unit = "log10 (m)",
    #                                                                                                                        taxa = names(veg.height), ranks = "Species", values = as.vector(veg.height)),
              Plant_height = TraitDataset(reference = "TRY", trait = "Plant height",
                                               kingdom = "Plantae", unit = "log10 (m)", taxa = names(veg.height),
                                               ranks = "Species", values = as.vector(veg.height))))
}

#### Offical PREDICTS trait assembleance ####
library(yarg)
# Trait data dir
dataDir <- "../../Data/Traits/"

cat('Adding trait data\n') #need to get files and correct pathways
bentley<-BentleyArthropods(paste(dataDir,"Ryan Arth  4.8.15.csv",sep=""))
edgar<-EdgarBeetleLength(paste(dataDir,"EdgarBeetleLength-2014-01-09.csv",sep=""))
horton <- HortonArachnid(paste(dataDir,"HortonArachnids-2014-05-22.csv",sep=""))
su<-SuFormicidae(paste(dataDir,"SuFormicidaeBodyLength-2014-05-23.csv",sep=""))
pantheria <- PantheriaMammals(path = paste(dataDir,"PanTHERIA_1-0_WR05_Aug2008.txt",sep=""))
cooper <- CooperHerptileSVL(paste(dataDir,"Cooper_svl_data.csv",sep=""))
senior <- SeniorHerptileSVL(paste(dataDir,"Senior_svl_data.csv",sep=""))
myhrvold <- MyhrvoldAmniotes(paste(dataDir,"Amniote_Database_Aug_2015.csv",sep=""))
kissling<-KisslingMammals(paste(dataDir,"MammalDIET_v1.0.txt",sep=""))
#sekerciog<-SekerciogluAvianDiet(paste(dataDir,"SekercoigluDiet.csv",sep=""))
BirdDietBM = EltonianTraitsBirds(paste0(dataDir,"Eltonian_Foraging/BirdFuncDat.txt"))
MammalDietBM = EltonianTraitsMammals(paste0(dataDir,"Eltonian_Foraging/MamFuncDat.txt"))
PlantHeightTry = TRY2(paste0(dataDir,"PlantsTRY_05042016123306/1886.txt"))
bentleymissing <-BentleyMissingTaxa(paste(dataDir,"Vertebrate Families.csv",sep=""))

## Add trophic level ##
contrb = data.frame(SS = unique(r$SS))
#r2 = r # Security copy
r = r2 # Copy
r$Adult_Trophic_level_Category <- bentley$Adult_Trophic_Level$Values$Value[
  match(r$Family,bentley$Adult_Trophic_Level$Values$Taxon)]
n.bentley <- rowSums( table(r$SS,r$Adult_Trophic_level_Category) >0)
length( which(n.bentley>0) )
contrb$bentley <- NA
contrb$bentley[ as.character(contrb$SS) %in% names(which(n.bentley>0)) ] <- 1
# For how many studies and which did that add
r = r2 # Copy
r <- AddTraits(r,edgar$Trophic_level)
r <- AddGenusInterpolatedCategoricalTraits(r,0.95,edgar$Trophic_level)
n.edgar <- rowSums( table(r$SS,r$Adult_Trophic_level_Category) >0)
length( which(n.edgar>0) )
contrb$edgar <- NA
contrb$edgar[ as.character(contrb$SS) %in% names(which(n.edgar>0)) ] <- 1
# - #
r = r2 # Copy
r <- AddTraits(r,su$Trophic_level)
r <- AddGenusInterpolatedCategoricalTraits(r,0.95,su$Trophic_level)
n.su <- rowSums( table(r$SS,r$Trophic_level_Category) >0)

contrb$su <- NA
contrb$su[ as.character(contrb$SS) %in% names(which(n.su > 0)) ] <- 1
saveRDS(contrb,"999_CoAuthorShipMatch.rds")

r <- AddTraits(r,kissling$Trophic_level) 
r <- AddGenusInterpolatedCategoricalTraits(r,0.95,kissling$Trophic_level)

#r <- AddTraits(r,sekerciog$Trophic_level) 
#r <- AddGenusInterpolatedCategoricalTraits(r,0.95,sekerciog$Trophic_level)
# Replace sekercioglu with EltonianTraits
r <- AddTraits(r,BirdDietBM$Trophic_level) 
r <- AddGenusInterpolatedCategoricalTraits(r,0.95,BirdDietBM$Trophic_level)

r$Adult_Trophic_level_Category[is.na(r$Adult_Trophic_level_Category)] <-
  r$Trophic_level_Category[is.na(r$Adult_Trophic_level_Category)] 

r$Missing_Trophic_level <- bentleymissing$Values$Value[match(r$Family,bentleymissing$Values$Taxon)]

r$Adult_Trophic_level_Category[is.na(r$Adult_Trophic_level_Category)]<-
  r$Missing_Trophic_level[is.na(r$Adult_Trophic_level_Category)] 

# Assign autotroph to plants
r$Adult_Trophic_level_Category <- factor(r$Adult_Trophic_level_Category,levels=c(
  levels(r$Adult_Trophic_level_Category),"Autotroph"))
r$Adult_Trophic_level_Category[(r$Kingdom=="Plantae")]<-"Autotroph"

# Fungi being fungi
r$Adult_Trophic_level_Category<-factor(r$Adult_Trophic_level_Category,levels=c(
  levels(r$Adult_Trophic_level_Category),"Fungi"))
r$Adult_Trophic_level_Category[(r$Kingdom=="Fungi")]<-"Fungi"

# Still Missing
yarg:::.Log("Still Missing", length(which(is.na(r$Adult_Trophic_level_Category))) / nrow(r)) # 0.1073
# Cleanup
r$Missing_Trophic_level <- NULL 
r$Trophic_level_Category <- NULL

r$ThermStrat <- factor(NA,levels=c("Endothermy","Ectothermy"))
r$ThermStrat[(r$Phylum!="") & (r$Phylum!="Chordata")]<-"Ectothermy"
r$ThermStrat[(r$Class=="Amphibia") | (r$Class=="Reptilia")]<-"Ectothermy"
r$ThermStrat[(r$Class=="Aves") | (r$Class=="Mammalia")]<-"Endothermy"

### Now do bodymass ###
r$Adult_Minimum_Mass_g <- bentley$Adult_Min_Mass$Values$Value[
  match(r$Family,bentley$Adult_Min_Mass$Values$Taxon)]
r$Adult_Maximum_Mass_g <- bentley$Adult_Max_Mass$Values$Value[
  match(r$Family,bentley$Adult_Max_Mass$Values$Taxon)]

# Calculate mean adult mass
r$Adult_Mean_Mass_g <- 10^((log10(r$Adult_Minimum_Mass_g)+
                                      log10(r$Adult_Maximum_Mass_g))/2)

r <- AddTraits(r,edgar$Length_derived_mass,
                       horton$Length_derived_mass,
                       su$Length_derived_mass,pantheria$Adult_wet_mass,
                       cooper$Mass,senior$Mass,myhrvold$Adult_wet_mass,
                       BirdDietBM$Adult_wet_mass,MammalDietBM$Adult_wet_mass,
                       PlantHeightTry$Plant_height)

r <- AddGenusAveragedTraits(r,edgar$Length_derived_mass,
                                    horton$Length_derived_mass,
                                    su$Length_derived_mass,
                                    pantheria$Adult_wet_mass,
                                    cooper$Mass,senior$Mass,
                                    myhrvold$Adult_wet_mass,
                                    BirdDietBM$Adult_wet_mass,MammalDietBM$Adult_wet_mass,
                                    PlantHeightTry$Plant_height
                            )


# Combine estimates into one column
r$mass <- 10^(r$Adult_wet_mass_log10_g)
r$mass[is.na(r$mass)] <- r$Length_derived_mass_g[is.na(r$mass)]
r$mass[is.na(r$mass)] <- r$Adult_Mean_Mass_g[
  is.na(r$mass)]
r$mass[is.na(r$mass)] <- 10^(r$Adult_wet_mass_log10_g[
  is.na(r$mass)])
r$mass[is.na(r$mass)] <- 10^(r$Plant_height_log10_m[is.na(r$mass)]) * 100 # to cm for better equivalenz to g

r$mass_log10 <- log10(r$mass+1)
# How many missing
yarg:::.Log("Still Missing", length(which(is.na(r$mass))) / nrow(r)) # 0.284
#clean up
r$Plant_height_log10_m <- NULL
r$Adult_wet_mass_log10_g <- NULL
r$Length_derived_mass_g <- NULL
r$Adult_Mean_Mass_g <- NULL
r$Adult_Maximum_Mass_g <- NULL
r$Adult_Minimum_Mass_g <- NULL


# ------------------------------------- #
# Now for each study assemble bodymass and trophic groups
# But first remove absence rows from the raw data
r <- subset(r,Measurement>0)
r <- subset(r, TGrouping != "Fungi")
# Descriptive stats
psych::describeBy(r$mass,r$TGrouping)
psych::describeBy(r$mass_log10,r$TGrouping)


# Tim's massbins
r$MassBin.Tim <- cut(r$mass,breaks=c(0,2,20,100,5000000),labels=c(1,2,3,4))
# My Massbins
hist(r$mass_log10)
r$MassBin.Own <- cut(r$mass_log10,breaks=c(0,1,2,7),labels=c(1,2,3))

# Overall terciles 
rt <- r %>% dplyr::summarise(low.33 = quantile(mass,.33,na.rm=T),
                   mid.66 = quantile(mass,.66,na.rm=T),
                   high.1 = quantile(mass,1,na.rm=T))
r$MassBin.Tercile <- cut(r$mass,breaks=c(0,rt$low.33,rt$mid.66,rt$high.1),labels=c(1,2,3))

# Terciles per taxonomic group
r <- r %>% dplyr::select(TGrouping,mass) %>% 
  dplyr::group_by(TGrouping) %>% 
  dplyr::summarise(low.33 = quantile(mass,.33,na.rm=T),
                   mid.66 = quantile(mass,.66,na.rm=T),
                   high.1 = quantile(mass,1,na.rm=T)) %>% ungroup() %>% 
  # Join back with Tgrouping
  full_join(r,.,by="TGrouping")
# Cut manually
r$MassBin.TGTercile <- ifelse(r$mass > 0 & r$mass <= r$low.33,1,
                        ifelse(r$mass > r$low.33 & r$mass <= r$mid.66, 2,
                               ifelse(r$mass > r$mid.66,3,NA)
                        )
                    )
r$MassBin.TGTercile <- factor(r$MassBin.TGTercile)
table(r$MassBin.TGTercile)
r$low.33 <- NULL; r$mid.66 <- NULL; r$high.1 <- NULL


### Now aggregate per study ###
res <- data.frame()

## Overall
ag <- aggregate(list(MassBin.Tim = r$MassBin.Tim),list(SS = r$SS),function(x) names(which.max( table(x) )))
res <- ag

# MassBin.Tim * Measurement
rt <- r %>% dplyr::select(SS,Measurement,MassBin.Tim) %>% 
  dplyr::group_by(SS) %>% 
  dplyr::summarise(nb.1 = sum(Measurement[which(MassBin.Tim=="1")]),
                   nb.2 = sum(Measurement[which(MassBin.Tim=="2")]),
                   nb.3 = sum(Measurement[which(MassBin.Tim=="3")]),
                   nb.4 = sum(Measurement[which(MassBin.Tim=="4")]),
                   rNA = sum(Measurement[which(is.na(MassBin.Tim))])
  )
# Removal criterium. 
# Kick out trophic studies where more missing recordings are present than non missing
rt <- rt[which(apply(rt[,c("nb.1","nb.2","nb.3","nb.4")],1,sum) > rt$rNA),]
rt$rNA <- NULL
res <- full_join(res,data.frame(SS = rt$SS,MassBin.Tim.abd = apply(rt[,-1],1,which.max) ),by="SS")

# For my own groups
ag <- aggregate(list(MassBin.Own = r$MassBin.Own),list(SS = r$SS),function(x) names(which.max( table(x) )))
res <- full_join(res,ag,by="SS")

# MassBin.Own * Measurement
rt <- r %>% dplyr::select(SS,Measurement,MassBin.Own) %>% 
  dplyr::group_by(SS) %>% 
  dplyr::summarise(nb.1 = sum(Measurement[which(MassBin.Own=="1")]),
                   nb.2 = sum(Measurement[which(MassBin.Own=="2")]),
                   nb.3 = sum(Measurement[which(MassBin.Own=="3")]),
                   rNA = sum(Measurement[which(is.na(MassBin.Own))])
  )
# Removal criterium. 
# Kick out trophic studies where more missing recordings are present than non missing
rt <- rt[which(apply(rt[,c("nb.1","nb.2","nb.3")],1,sum) > rt$rNA),]
rt$rNA <- NULL
res <- full_join(res,data.frame(SS = rt$SS,MassBin.Own.abd = apply(rt[,-1],1,which.max) ),by="SS")


## The same for terciles
ag <- aggregate(list(MassBin.Tercile = r$MassBin.Tercile),list(SS = r$SS),function(x) names(which.max( table(x) )))
res <- full_join(res,ag,by="SS")
# MassBin * Measurement
rt <- r %>% dplyr::select(SS,Measurement,MassBin.Tercile) %>% 
  dplyr::group_by(SS) %>% 
  dplyr::summarise(nb.1 = sum(Measurement[which(MassBin.Tercile=="1")]),
                   nb.2 = sum(Measurement[which(MassBin.Tercile=="2")]),
                   nb.3 = sum(Measurement[which(MassBin.Tercile=="3")]),
                   rNA = sum(Measurement[which(is.na(MassBin.Tercile))])
  )
# Kick out trophic studies where more missing recordings are present than non missing
rt <- rt[which(apply(rt[,c("nb.1","nb.2","nb.3")],1,sum) > rt$rNA),]
rt$rNA <- NULL

res <- full_join(res,data.frame(SS = rt$SS,MassBin.Tercile.abd = apply(rt[,-1],1,which.max) ),by="SS")

# Terciles within group
ag <- aggregate(list(MassBin.TGTercile = r$MassBin.TGTercile),list(SS = r$SS),function(x) names(which.max( table(x) )))
res <- full_join(res,ag,by="SS")
# MassBin * Measurement
rt <- r %>% dplyr::select(SS,Measurement,MassBin.TGTercile) %>% 
  dplyr::group_by(SS) %>% 
  dplyr::summarise(nb.1 = sum(Measurement[which(MassBin.TGTercile=="1")]),
                   nb.2 = sum(Measurement[which(MassBin.TGTercile=="2")]),
                   nb.3 = sum(Measurement[which(MassBin.TGTercile=="3")]),
                   rNA = sum(Measurement[which(is.na(MassBin.TGTercile))])
  )
# Kick out trophic studies where more missing recordings are present than non missing
rt <- rt[which(apply(rt[,c("nb.1","nb.2","nb.3")],1,sum) > rt$rNA),]
rt$rNA <- NULL

res <- full_join(res,data.frame(SS = rt$SS,MassBin.TGTercile.abd = apply(rt[,-1],1,which.max) ),by="SS")

# Finally calculate mean and sd overall mass per study
rt <- r %>% dplyr::group_by(SS) %>% 
  dplyr::summarise(Mass.SSavg = mean(mass,na.rm=T),
                   Mass.SSSD = sd(mass,na.rm = T))
res <- full_join(res,rt,by="SS")

# And also for dominant trophic structure
ag <- aggregate(list(Trophic_SS = r$Adult_Trophic_level_Category),list(SS = r$SS),function(x) names(which.max( table(x) )))
res <- full_join(res,ag,by="SS")

# Trophic * Measurement
rt <- r %>% dplyr::select(SS,Measurement,Adult_Trophic_level_Category) %>% 
  dplyr::group_by(SS) %>% 
  dplyr::summarise(Herbivore = sum(Measurement[which(Adult_Trophic_level_Category=="Herbivore")]),
                   Omnivore = sum(Measurement[which(Adult_Trophic_level_Category=="Omnivore")]),
                   Carnivore = sum(Measurement[which(Adult_Trophic_level_Category=="Carnivore")]),
                   Fungivore = sum(Measurement[which(Adult_Trophic_level_Category=="Fungivore")]),
                   Detritivore = sum(Measurement[which(Adult_Trophic_level_Category=="Detritivore")]),
                   Autotroph = sum(Measurement[which(Adult_Trophic_level_Category=="Autotroph")]),
                   Non.feeding = sum(Measurement[which(Adult_Trophic_level_Category=="Non-feeding")]),
                   Fungi = sum(Measurement[which(Adult_Trophic_level_Category=="Fungi")]),
                   rNA = sum(Measurement[which(is.na(Adult_Trophic_level_Category))])
                   )

rt <- rt[which(rt %>% dplyr::select(Herbivore:Fungi) %>% apply(.,1,sum) > rt$rNA),]
rt$rNA <- NULL

res <- full_join(res,data.frame(SS = rt$SS,Trophic_SS.abd = apply(rt[,-1],1,function(x)names(which.max(x)) ) ),by="SS")

# Finally also get biomass (abundance weighted bodymass)
rt <- r %>% dplyr::group_by(SS) %>% 
  dplyr::filter(Diversity_metric_type == "Abundance") %>% 
  dplyr::summarise(Biomass = weighted.mean(mass,w = Measurement,na.rm=T))
res <- full_join(res,rt,by="SS")

# Filter all those out with missing mean bodymass (grouping associated with non-data == 0)
#res <- res[which(!is.nan(res$Mass.SSavg)),]
# Then save
saveRDS(res,"resSaves/MassBinStudy.rds")



