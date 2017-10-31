# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
# Custom packages
library(marfunky)
library(ggplot2)
library(scales)
myLog <- function(...) {
  cat(paste0("[Extract] ", Sys.time(), " | ", ..., "\n"))
}

sites <- readRDS("PlantStudiesOnFocus.rds")
sites$TGrouping <- factor(sites$TGrouping)
sites$startyear <- factor(sites$startyear)

# What is the average spread of sampling startyear for all taxonomic groups
library("gplots")
dt <- table(sites$startyear,sites$TGrouping)
balloonplot(t(dt), main ="", xlab ="", ylab="Starting year\n of sampling",
            label = T,label.lines=T,
            show.margins = T,show.zeros=F,cum.margins=F)

g <- sites %>% filter(TGrouping == "Plantae") %>% 
  ggplot(aes(startyear,group=TGrouping)) + theme_minimal() + 
    geom_density(alpha=.5,stat = "count") + 
    labs(x="",y="Number of Plant sites") +
    scale_y_continuous(breaks = pretty_breaks(7))
ggsave("Plots/NumberofPlantsSites.png",plot=g)

# Calculate average time of sampling
temp <- sites %>% 
  mutate(Sample_start_earliest = as.Date(ymd(Sample_start_earliest)),
         Sample_end_latest = as.Date(ymd(Sample_end_latest))) %>% # correct sampling starts
  select(SSBS,Sample_start_earliest,Sample_end_latest) %>% 
  group_by(SSBS) %>% 
  summarize(TL = time_length(interval(Sample_start_earliest, Sample_end_latest),unit = "day") )
# Join back
sites <- left_join(sites,temp)

# Make a plot showing the average sampling timespan per taxonomic group and startyear
g <- sites %>% mutate(TL = TL / 30) %>% filter(TGrouping == "Plantae") %>% 
  #filter(TL <= 12) %>% # Get only sites within a year
  ggplot(aes(TL)) + theme_light() +
    geom_density(stat="count") +
    scale_y_continuous(breaks = pretty_breaks(5)) +
    scale_x_continuous(breaks = pretty_breaks(6),minor_breaks = pretty_breaks(3)) +
    facet_wrap(~UN_subregion,scales = "free") +
    labs(x = "Months",y="Total number of plant sites sampled that long", title = "PREDICTS plant sites - Total sampling duration \ (in months)")
g
ggsave("Plots/PlantSamplingDurationUNRegion.png",plot=g)

# What is the time of sampling spread over regions
sites$Lat2 <- round(sites$Latitude,0)
sites$TL2 <- sites$TL / 30
s <- sites %>% group_by(TGrouping,Lat2) %>% summarize(N = mean(TL2))

g <- ggplot(s,aes(x=Lat2,y=N,fill=TGrouping)) + theme_light()+
  geom_bar(stat="identity",position = position_dodge()) +
  scale_x_continuous(breaks=pretty_breaks(10)) +
  facet_wrap(~TGrouping) +
  coord_flip() +
  labs(x = "Latitude (rounded)", y = "Average sampling duration (months)")
ggsave("Plots/LatitudeAverageMonth.png",plot=g)


#### Table - How many sites per cell ####
library(gridExtra)
plants <- readRDS("PlantStudiesOnFocus.rds") %>% 
  dplyr::select(SS,SSBS,cells) %>% 
  group_by(SS) %>% 
  summarise(N_Site = length(unique(SSBS)),
            N_Cells = length(unique(cells))) %>% 
  mutate(Dif = N_Site / N_Cells) %>% 
  arrange(desc(Dif))
# Get only suitable studies
a = plants %>% filter(Dif == 1)

names(plants) <- c("Plant.Study","Number.unique.sites","Number.unique.cells",expression(Sites/Cells))
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(head(plants,33), rows=NULL, theme=tt)
grid.arrange(tbl,as.table=T)


#### What is the disribution off inversed values ####
# Make time series for all plots
library(data.table)
library(zoo)
sites <- readRDS("PlantStudiesOnFocus.rds")
qs <- as.data.table(readRDS("Extracts/PREDICTS_AlbedoOverallQUality.rds")) # Quality scores

s <- sites %>% dplyr::filter(SSBS == unique(sites$SSBS)[1])
# Construct the interval of interest
i <- interval("2000-01-01",unique(s$Sample_start_earliest))
sub_qs <- subset(qs,SSBS %in% s$SSBS)
sub_qs <- sub_qs[which(ymd(sub_qs$date) %within% i),]
sub_qs$date <- as.Date(ymd(sub_qs$date),origin="2000-01-01")

sub1_qs <- zooreg(sub_qs$value,order.by = sub_qs$date,frequency = 16)
out <- sub1_qs

for(id in unique(sites$SSBS)[-1]){
  myLog(id)
  # Create a subset
  s <- sites %>% dplyr::filter(SSBS == id)
  # Construct the interval of interest
  i <- interval("2000-01-01",unique(s$Sample_start_earliest))
  sub_qs <- subset(qs,SSBS %in% s$SSBS)
  sub_qs <- sub_qs[which(ymd(sub_qs$date) %within% i),]
  sub_qs$date <- as.Date(ymd(sub_qs$date),origin="2000-01-01")

  sub1_qs <- zooreg(sub_qs$value,order.by = sub_qs$date,frequency = 16)
  out <- merge(out,sub1_qs,suffixes=id)
}  


library(ggfortify)
library(scales)
a <- zoo(x = rowSums(coredata(out)==1,na.rm = T) / apply(coredata(out),1,function(x) length(which(!is.na(x)))) ,order.by = index(out))
autoplot(a) + theme_light() + labs(x ="", y = "Proportion of bad observations (QA Flag = 1)\n of all observations",title="MODIS BRDF Overall quality") +scale_y_continuous(breaks = pretty_breaks())
