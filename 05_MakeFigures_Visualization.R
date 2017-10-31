#### Package and data loading ####
library(ggplot2)
library(ggthemes)
library(ggsci)
library(scales)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(marfunky)
library(lubridate)
library(stringr)
library(reshape2)
library(sp)
library(maps)
library(tmap)
library(rgeos)
library(broom) # New. Force using it
par.ori <- par(no.readonly = T)
OutputPath = "PubFigures/"
PermPath ="permutations/"
sites <- readRDS("AllPredictsSites_full.rds") # All
sites <- subset(sites, startyear > 2000)

# # Get the ID of those studies with clear site / pixel ratio
# ID <- sites %>% 
#   dplyr::select(SS,SSBS,cells) %>% 
#   group_by(SS) %>% 
#   summarise(N_Site = length(unique(SSBS)),
#             N_Cells = length(unique(cells))) %>% 
#   mutate(Dif = N_Site / N_Cells) %>% 
#   arrange(desc(Dif)) %>% 
#   dplyr::select(SS,Dif)
# sites <- left_join(sites,ID)
# 
# # Exclude Fungi from the analysis due to only two studies in total
# sites <- dplyr::filter(sites,TGrouping != "Fungi")

# Sites with zero estimates
stop("Don't run this")
spast <- readRDS("Center_FullTimeSeriesSmooth_Final.rds")
ss <- lapply(spast,FUN = function(x) which( coredata(x[["EVI2"]])==0) )
unlist(ss)
ss[["SH1_2013__Peri 1  691"]]

# ------------------------------------------- #

#### Figure 1 - Composite map with barcharts ####

# Subset to only studies in analysis
# f.permut <- readRDS( list.files(PermPath,pattern = "s_metricSorVeg_permute_1.rds")[1] ) %>% 
#   dplyr::select(SS) %>% distinct()
# New approach. Based on fitted studies
f.permut <- readRDS("resSaves/outBCpermute_minStart5yEVI2_BC_BCVeg.rds")$StudyCoeff %>% 
  dplyr::select(SS) %>% distinct()

s <- left_join(f.permut,sites,by="SS") %>% 
  dplyr::select(SS,SSBS) %>% group_by(SS) %>% summarise(N = n())
s <- left_join(s, sites %>% dplyr::select(SS,Longitude,Latitude) ) %>% distinct %>% as.data.frame %>% 
  # Calculate mean for simple centroid
  aggregate(.~SS, data=., mean) %>% left_join(., sites %>% dplyr::select(SS,TGrouping)) %>% distinct %>% as.data.frame
coordinates(s) = ~Longitude+Latitude
proj4string(s) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

s$TGrouping <- factor(s$TGrouping,
                  levels = c("Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))

data(World)
w <- subset(World, continent != "Antarctica")

# color vector
levels(s$TGrouping)
cv <- c("#3B5C00","tan1","lightskyblue","plum","#5D959A","#D4605A")


mp <- tm_shape(w, projection="eck4") + 
  tm_grid(projection="longlat",n.x = 10, n.y = 10, labels.size = 0,lwd = .5) +
  tm_polygons(col="#FFF8DC",alpha = 1,border.alpha = .25,border.col = "black") + 
  tm_shape(s) + tm_bubbles(col="TGrouping", title.col="",palette = cv,scale = 1.5,
                           border.col = "black", border.alpha = .5,
                           legend.col.is.portrait=F,size=.3,
                           legend.hist = F,legend.hist.title = "") +
  tm_format_World_wide(inner.margins=c(.04,.03, .02, .01),outer.margins=c(.01,.01,.01,.01), 
                       title.size=2,title = paste0(""),title.bg.color="white",
                       title.position=c("right","top"),title.bg.alpha=0.75,frame=F,
                       bg.color="white",earth.boundary = T,earth.boundary.lwd =.5,space.color="transparent",
                       legend.show=F)
mp
save_tmap(mp, paste0(OutputPath,"GeographicCoverage.png"), width=3840, height=2160,dpi = 400)

#### Figure 1 - Chart - Taxonomic coverage ####
# Taxonomic coverage 
data(World)
sp <- sites %>% dplyr::select(SS,Longitude,Latitude) %>% as.data.frame
coordinates(sp) = ~Longitude+Latitude
proj4string(sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
sp <- spTransform(sp,CRS(proj4string(World)))
sites$continent <- over(sp,World)$continent
# Manually set small islands
sites$continent[which(sites$Country=="Solomon Islands")] <- "Oceania"
sites$continent[which(sites$Country=="New Zealand")] <- "Oceania"
sites$continent[which(sites$Country=="South Africa")] <- "Africa"
sites$continent[which(sites$Country=="Mozambique")] <- "Africa"
sites$continent[which(sites$Country=="Japan")] <- "Asia"
sites$continent[which(sites$Country=="Sweden")] <- "Europe"
sites$continent[which(sites$Country=="Comoros")] <- "Africa"
sites$continent[which(sites$Country=="Portugal")] <- "Europe"
sites$continent[which(sites$Country=="Equatorial Guinea")] <- "Africa"
sites$continent[which(sites$Country=="Argentina")] <- "South America"
sites$continent[which(sp$SS =="AD1_2010__Davis 1")] <- "Europe"
sites$continent[which(sp$SS =="AD1_2010__Redpath 1")] <- "Europe"
sites$continent[which(sp$SS =="AD1_2011__Cameron 1")] <- "North America"
sites$continent[which(sp$SS =="AD1_2011__Cameron 2")] <- "North America"
sites$continent[which(sp$SS =="AD1_2011__Tonietto 1")] <- "North America"
sites$continent[which(sp$SS =="AD1_2012__Osgathorpe 2")] <- "Europe"
sites$continent[which(sp$SS =="AD1_2012__Yoon 1")] <- "Asia"
sites$continent[which(sp$SS =="DL1_2010__Luskin 1")] <- "Oceania"
sites$continent[which(sp$SS =="DL1_2012__Dallimer 1")] <- "Africa"
sites$continent[which(sp$SS =="GP1_2011__Andersen 1")] <- "Oceania"
sites$continent[which(sp$SS =="GP1_2011__Andersen 2")] <- "Oceania"
sites$continent[which(sp$SS =="HP1_2010__Krauss 1")] <- "Europe"
sites$continent[which(sp$SS =="JD1_2010__Sodhi 1")] <- "Asia"
sites$continent[which(sp$SS =="MH1_2010__CATIE 4")] <- "South America"
sites$continent[which(sp$SS =="SC2_2012__Giordani 1")] <- "Europe"
sites$continent[which(sp$SS =="VB1_2013a_Jones 1")] <- "Europe"
sites$continent[which(sp$SS =="VB1_2013__Jones 1")] <- "Europe"
sites$continent[which(sp$SS =="VB1_2013__Jones 2")] <- "Europe"
sites$continent[which(sp$SS =="AD1_2011__Nielsen 1")] <- "Europe"
sites$continent[which(sp$SS =="AD1_2011__Nielsen 2")] <- "Europe"
sites$continent[which(sp$SS =="AD1_2011__Nielsen 3")] <- "Europe"
sites$continent[which(sp$SS =="DI1_2005__Gove 1")] <- "North America"
sites$continent[which(sp$SS =="DI1_2005__Gove 2")] <- "North America"
sites$continent[which(sp$SS =="HP1_2008__Ranganathan 1")] <- "Asia"
sites$continent[which(sp$SS =="HP1_2010__Krauss 4")] <- "Europe"
sites$continent[which(sp$SS =="HP1_2010__Lasky 1")] <- "North America"
sites$continent[which(sp$SS =="HP1_2010__Lasky 2")] <- "North America"
sites$continent[which(sp$SS =="AD1_2011__Nielsen 4")] <- "Europe"
sites$continent[which(sp$SS =="SC1_2006__UrbinaCardona 1")] <- "North America"

#unique(sp$SS[which(is.na(sites$continent))])
#sites$Country[which(sites$SS == "SC1_2006__UrbinaCardona 1")]

# Calculate groupings
s <- left_join(f.permut,sites,by="SS") %>% dplyr::select(SS,continent,TGrouping) %>% 
  group_by(continent,TGrouping) %>% summarise(N = n()) %>% mutate(freq = N / sum(N))
# Correct ordering of levels to reverse
s$TGrouping <- factor(s$TGrouping,
                      levels = c("Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))
# Adapt continent labels for americas
s$continent <- as.character(s$continent)
s$continent[which(s$continent=="North America")] <- "North\n America"
s$continent[which(s$continent=="South America")] <- "South\n America"
s$continent <- factor(s$continent, levels = c("North\n America","South\n America","Europe","Africa","Asia","Oceania"))

# Calculate number of studies per continent
N <- left_join(f.permut,sites,by="SS") %>% dplyr::select(SS,continent) %>% 
  distinct() %>% group_by(continent) %>% summarise(Ns = n()) %>% mutate(continent = as.character(continent))
N$continent[which(N$continent=="North America")] <- "North\n America"
N$continent[which(N$continent=="South America")] <- "South\n America"
N$continent <- factor(N$continent, levels = c("North\n America","South\n America","Europe","Africa","Asia","Oceania"))
s <- left_join(s,N)

# Calculate total per continent
sN = s %>% group_by(continent) %>% summarise(fs = sum(unique(Ns)))
sN$lab <- paste(sN$continent,"\n(",sN$fs,")")
s <- left_join(s,sN,by="continent")
s$lab <- factor(s$lab,levels = c("North\n America \n( 42 )","South\n America \n( 21 )","Europe \n( 63 )",
                                 "Africa \n( 29 )","Asia \n( 20 )","Oceania \n( 23 )"))
gt <- ggplot(s,aes(x = lab, y = Ns, fill=TGrouping) ) + theme_classic(base_size = 18) +
  geom_bar(position = "fill",stat = "identity") + 
  scale_fill_manual(values = cv, # Color scale used above
                    guide_legend(title="Taxonomic groups")) +
  scale_y_continuous(expand = c(0,0),labels = percent_format(),position = "right") +
#  annotate("text",aes(x = continent,y = 1.05,values=Ns)) + 
  labs(x="",y="Proportion of studies") +
  theme(axis.text.x = element_text(vjust = .5), axis.text.y = element_text(color="#696969"),
        legend.position = "left")
gt
ggsave( paste0(OutputPath,"TaxonomicCoverage.png"),plot=gt)

## --------------- Plot Figure 2-4 --------------- ##
indic = "EVI2"#
dism = "BCVeg"# #Dissimilarity metric
metric =  "BC" #"_meanSPD"# "correctedBC"

res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds"))
#res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5yINDIVIDUAL",indic,"_",metric,"_",dism,".rds"))

out <- res.pred$PermRuns

# ---- # 
# Summary stats to be included
out$py <- factor(out$py)
results <- out %>% dplyr::filter(Time != "PresentPast") %>% 
  group_by(Model,Type,Time,py,xx) %>% 
  dplyr::summarise(avg.x = mean(x),
                   avg.fit = mean(fit),
                   avg.se = mean(se),
                   fit.low = (mean(fit) - mean(se)),
                   fit.high = (mean(fit) + mean(se)),
                   min.N = max(N.studies))

#mres$Type <- factor(mres$Type,levels = rev(levels(mres$Type)))
results$Model <- factor(results$Model,
                        levels = c("All","Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))
results$Time <- droplevels(results$Time)
levels(results$Time) <- c("Sampling conditions","Past conditions")

#results$N <- N$N[match(results$Model,N$TGrouping)]
lab = out %>% group_by(Model) %>% summarise(N = min(N.studies))
lab$label <- factor(paste0(lab$Model,"\n","(N = ",lab$N,")"))
results$label <- lab$label[match(results$Model,lab$Model)]
results$label <- factor(results$label,levels = c("All\n(N = 187)","Plantae\n(N = 37)","Invertebrates\n(N = 98)","Amphibia\n(N = 8)","Reptilia\n(N = 4)",
                                                 "Aves\n(N = 21)","Mammalia\n(N = 15)"))

# Add marginal plots 
# Need to sample the original distribution of values within permutations
f.permut <- sort(list.files(PermPath,pattern = paste0("s_metric",dism,"_permute"),full.names = T))
res.pred <- list() # Output for all predictions
# Load matrices and permSam based on 10% random of the original
permSamp <- data.frame()
pp = 1
y = 5
res_yearsampling <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"samplingperiod_0y-matrix.rds")) # The year of before start
res_yearbefore <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"pastperiod_5y-matrix.rds")) # The years before
for(p in sample(1:100,pp)){
    myLog("Running permutation ",p)
    # Load permutation
    s.metric <- readRDS(f.permut[p])

    # Available past data
    focal.period  = ymd("2000.02.18") + (365*( y +1)) # Available past in total
    # Add in startdate
    s.metric$startdate_x <- sites$Sample_start_earliest[match(s.metric$SSBS_x,sites$SSBS)]
    s.metric$startdate_y <- sites$Sample_start_earliest[match(s.metric$SSBS_y,sites$SSBS)]
    # Calculate distance in days and remove sites with greater dist of 3 months (~90days)
    s.metric <- mutate(s.metric,startd = as.numeric( abs(startdate_x - startdate_y))) %>%
      dplyr::filter(startd <= 90) %>% # 3 months (~ 90 days)
      dplyr::filter(startdate_x > focal.period & startdate_y > focal.period) %>%  # Must have sufficient historical coverage
      # Filter to only those studies that are within the maximum of ylist, thus equalizing nr of studies for all years
      dplyr::filter(startdate_x > (ymd("2000.02.18") + (365*(max(ylist)+1))) & startdate_y > (ymd("2000.02.18") + (365*(max(ylist)+1)))) %>%
      dplyr::select(-startdate_x,-startdate_y,-startd) # remove again

    # Subset and assemble the dataset of the permutations and full matrix
    per.df <- data.frame()
    for(study in unique( s.metric$SS) ){
      sub <- subset(s.metric,SS == study)
      if(!(study %in% names(res_yearbefore))) next()
      # Melt past and sampling cond

      m <- reshape2::melt( res_yearsampling[[as.character(study)]] ) %>% dplyr::rename(SSBS_x = Var1,SSBS_y = Var2,EVIpres = value)
      m2 <- reshape2::melt( res_yearbefore[[as.character(study)]] ) %>% dplyr::rename(SSBS_x = Var1,SSBS_y = Var2,EVIpast = value)
      mm <- merge.data.frame(m2,m,by =c("SSBS_x","SSBS_y"),all.x =T) # Merge present with past
      rm(m,m2)
      # Merge with permutation to get only the permuteted pairs
      per.df <- rbind(per.df,
                      merge.data.frame(sub,mm,by =c("SSBS_x","SSBS_y"),all.x = T)
      )
    }
    permSamp <- rbind(permSamp,data.frame(p=p,value=per.df$value,EVIpres = per.df$EVIpres,EVIpast = per.df$EVIpast) )
  }

a =  results %>% dplyr::filter(Model == "All",py %in% c(0,5)) %>% ungroup()
a$Time <- factor(a$Time,labels= c("Sampling conditions","Past conditions (5 years ago)") )

# Figure 2 with rugs
g1 <- a %>% 
  # Get sampling conditions only once!
  ggplot(.,aes(x = avg.x, y = avg.fit, ymin = avg.fit - avg.se, ymax =  avg.fit + avg.se,
               group = Time ,color=Time,fill=Time)) + theme_few(base_size = 18) +
  # Add rugs as points
  #stat_density2d(data = permSamp, aes(x=EVIpres,y=value),inherit.aes = F, contour = T,color=ggthemes_data$colorblind[1]) +
  #stat_density2d(data = permSamp, aes(x=EVIpast,y=value),inherit.aes = F, contour = T,color=ggthemes_data$colorblind[2]) +
  #geom_rug(data = permSamp,aes(y=value),inherit.aes = F,sides = "l",show.legend = F) +
  # Rugs for LULC
  geom_rug(data = permSamp,aes(x=EVIpres),inherit.aes = F,color = ggthemes::ggthemes_data$colorblind[1],sides = "b",show.legend = F) +
  geom_rug(data = permSamp,aes(x=EVIpast),inherit.aes = F,color = ggthemes::ggthemes_data$colorblind[2],sides = "t",show.legend = F) +
  geom_line(lwd=1.05) + geom_ribbon(alpha=.1,show.legend = F) + 
  scale_x_continuous(breaks = pretty_breaks(7),limits = c(0, 0.9) ) +
  scale_y_continuous(breaks = pretty_breaks(7),expand = c(0,0),limits = c(0.6,0.98) ) +
  scale_fill_colorblind() +  scale_color_colorblind() +
  #labs(x = expression(BC[EVI]), y = expression(BC[Biodiversity]),title=paste0("All (N=", unique(a$min.N),")") ) +
  labs(x = "Difference between remote-sensing time series", y = "Difference between species assemblages") +
  # Make an arrow to indicate trend
  geom_segment(x=0.83, xend=0.83, y=0.89, yend=0.92,color="black",size=1.25,
               arrow=arrow(length=unit(2, "mm"))) +
  theme(legend.position = "bottom") + guides(fill="none",color="none") +
  theme(plot.title = element_text(hjust = 0.5)) #+ theme(plot.margin = unit(c(-0.5,-0.74,0,0), "lines")) 
g1
ggsave("PubFigures/Figure2.png",plot=g1,dpi = 400)

# Build vertical and horizontal density curves
library(gridExtra)
blankPlot <- ggplot() +
  geom_blank(aes(1,1)) +
  theme(line = element_blank(),
        text  = element_blank(),
        title = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
  )

xDense <- ggplot(permSamp %>% dplyr::select(EVIpres:EVIpast) %>% reshape2::melt(), aes(x=value,fill=variable)) + 
  geom_density(aes(y= ..count..),trim=F,alpha=.5) + 
  xlab("") + ylab("") + xlim(range(a$avg.x)) +
  theme_void() + scale_fill_colorblind()+
  #scale_x_continuous(limits = range(a$avg.x)) +
  theme(legend.position = "none",plot.margin = unit(c(0,-0.5,0,3.5), "lines"))

# Y margin density 
yDense <- ggplot(permSamp, aes(x=value)) + 
  geom_density(aes(y= ..count..),trim=F,alpha=.5) + 
  xlab("") + ylab("") + xlim(range(a$avg.fit)) + 
  coord_flip() + 
  theme_void() + scale_fill_manual(values="black") +
  #scale_x_continuous(limits = range(a$avg.fit)) +
  theme(legend.position = "none", plot.margin = unit(c(-1.25,4,4.5,0.268), "lines")) 


gg <- grid.arrange(xDense, blankPlot, g1, yDense, ncol=2, nrow=2, widths=c(3.5, 1.4), heights=c(1.4, 3))
ggsave("PubFigures/Figure2.png",plot=gg)

# Seperate one for the Taxonomic groups
g2 <- results %>% dplyr::filter(Model != "All",Time == "Sampling conditions") %>% 
  ggplot(.,aes(x = avg.x, y = avg.fit, ymin = avg.fit - avg.se, ymax = avg.fit + avg.se, group = Time,fill=Time,linetype = Time)) + theme_few() +
  geom_line() + geom_ribbon(stat="identity",alpha=.15) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  facet_wrap(~label,nrow=2,scales = "free_x") +
  scale_fill_colorblind() + 
  labs(x = expression(BC[Indicator]), y = expression(dBC[Biodiversity]),
       title = "") +
  theme(strip.text = element_text(size=12))
g2

# Combine them both
library(grid)
library(gridExtra)
library(cowplot)
grobs <- ggplotGrob(g1 + theme(legend.position="bottom",legend.direction = "vertical") )$grobs

pgrid <- plot_grid(g1 + theme(legend.position = "none"),
                   g2 + theme(legend.position = "bottom"),
                   ncol = 2,rel_widths = c(.4, .6),
                   labels = "AUTO")
# add legend
#p <- plot_grid(pgrid, legend_b)#, nrow = 2, rel_widths = c(1, .1))
#p

cowplot::ggsave(paste0(OutputPath,"Hyp3_ModelPastOnly.png"),plot=pgrid)
cowplot::ggsave(paste0(OutputPath,"Hyp3_ModelPastOnly.pdf"),plot=pgrid)

#### Figure 2 - Average relative coefficient + CI/SE overall per year ####
res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds"))
#res.pred <- readRDS(paste0("resSaves/outBCpermute_minStartALL__meanSPD_SorVeg.rds"))
#res.pred <- readRDS(paste0("resSaves/outBCpermute_minStartEVI2_BC_SorVeg.rds"))

results <- res.pred$Coeff

# Average
results2 <- results %>% dplyr::filter(Time != "PresentPast",py<=6,Model!="Fungi") %>% 
  group_by(Model,Time,py) %>% 
  dplyr::summarise(avg.eff = mean(eff),
                   avg.se = mean(se))

results2$xx <- as.factor(gsub("\\D","",results2$py))
results2$Model <- factor(results2$Model,
                        levels = c("All","Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))
results2$Time <- droplevels(results2$Time)
levels(results2$Time) <- c("Sampling conditions","Past conditions")

lab = results %>% group_by(Model) %>% summarise(N = min(N))
lab$label <- factor(paste0(lab$Model,"\n","(N = ",lab$N,")"))
results2$label <- lab$label[match(results2$Model,lab$Model)]
results2$label <- factor(results2$label,levels = c("All\n(N = 187)","Plantae\n(N = 37)","Invertebrates\n(N = 98)","Amphibia\n(N = 8)","Reptilia\n(N = 4)",
                                                 "Aves\n(N = 21)","Mammalia\n(N = 15)"))

# Plot
# First all
g1 <- results2 %>% dplyr::filter(Model == "All",Time == "Past conditions") %>% 
  ggplot(.,aes(x = xx, y = avg.eff, ymin = avg.eff - avg.se,ymax = avg.eff + avg.se,
                           color = Time)) +
    theme_few(base_size = 18) +
    scale_y_continuous(breaks = pretty_breaks(7)) +
    # Present
    geom_hline(data = dplyr::filter(results2,Model == "All",Time == "Sampling conditions"),aes(yintercept=avg.eff)) +
    # Upper and lower line
    geom_hline(data = dplyr::filter(results2,Model == "All",Time == "Sampling conditions"),aes( yintercept=avg.eff - avg.se),linetype="dotted") +
    geom_hline(data = dplyr::filter(results2,Model == "All",Time == "Sampling conditions"),aes( yintercept=avg.eff + avg.se),linetype="dotted") +
  
    geom_linerange(position = position_dodge(.75)) + geom_point(size = 3,position = position_dodge(.75)) +
    scale_color_manual(values = "#E69F00") +
#    geom_text(aes(x = xx,y = y,label = paste0(min.N ),fontface=2),data = res2label[which(res2label$Model=="All"),], inherit.aes = F,size=3) +
    theme(axis.text = element_text(color = "#808080"),axis.ticks.x = element_blank() ) +
    labs(x = "", y = substitute(paste("Average effect on ",dBC[Biodiversity])) ,title=lab$label[1]) +
    guides(color = "none")
g1
# Now the rest with wraps
g2 <- results2 %>% dplyr::filter(Model != "All",Time == "Past conditions") %>% 
  ggplot(.,aes(x = xx, y = avg.eff, ymin = avg.eff - avg.se,ymax = avg.eff + avg.se,
               color= Time)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  geom_hline(data = dplyr::filter(results2,Model != "All",Time == "Sampling conditions"),aes(yintercept=avg.eff)) +
  # Upper and lower line
  geom_hline(data = dplyr::filter(results2,Model != "All",Time == "Sampling conditions"),aes( yintercept=avg.eff - avg.se),linetype="dotted") +
  geom_hline(data = dplyr::filter(results2,Model != "All",Time == "Sampling conditions"),aes( yintercept=avg.eff + avg.se),linetype="dotted") +
  
  geom_linerange(position = position_dodge(.75)) + geom_point(size = 3,position = position_dodge(.75)) +
  scale_color_manual(values = "#E69F00") +
  #  geom_text(aes(x = xx,y = y,label = paste0(min.N ),fontface=2),data = res2label[which(res2label$Model!="All"),], inherit.aes = F,size=3) +
  facet_wrap(~label,nrow=2,scales="free_y") +
  theme(axis.text = element_text(color = "#808080"),axis.ticks.x = element_blank(),strip.text = element_text(size=12) ) +
  labs(x = "", y = substitute(paste("Average effect on ",dBC[Biodiversity]))) +
  guides(shape = "none",color="none")
g2

# Assemble
library(grid)
library(gridExtra)
library(cowplot)
grobs <- ggplotGrob(g1 + theme(legend.position="bottom",legend.direction = "vertical") )$grobs

pgrid <- plot_grid(g1 + theme(legend.position = "none"),
                   g2 + theme(legend.position = "bottom") + labs(y=""),ncol = 2,rel_widths = c(.4, .6),
                   labels = "AUTO")
pgrid = add_sub(pgrid,label = "Past period (years)",y  = 1.5, vjust = 0)


cowplot::ggsave(paste0(OutputPath,"Hyp3_CoeffCIPast.png"),plot=ggdraw(pgrid) )
cowplot::ggsave("Manuscript/Figure3.png",plot = ggdraw(pgrid))


#### SI - Figure 00X Current conditions line per group  ####
res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds"))
results <- res.pred$PermRuns
results$py[which(results$Time == "Present")] <- 0

# Average
results2 <- results %>% dplyr::filter(Time == "Present",Model!="Fungi") %>% 
  group_by(Model,xx,py) %>% 
  dplyr::summarise(avg.x = mean(x),
                   avg.fit = mean(fit),
                   avg.se = mean(se),
                   min.N = max(N.studies))

results2$Model <- factor(results2$Model,
                  levels = c("All","Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"),
                  labels = c("All\n(195|4044)","Plantae\n(41|775)","Invertebrates\n(102|1676)","Amphibia\n(9|142)","Reptilia\n(5|162)",
                             "Aves\n(23|969)","Mammalia\n(16|400)") # Copied from below
                  )

g1 <- results2 %>% dplyr::filter(Model == "All\n(195|4044)") %>% 
  ggplot(.,aes(x = avg.x, y = avg.fit, ymin= avg.fit - avg.se, ymax = avg.fit + avg.se)) + theme_few() +
  geom_line(lwd=.85) + geom_ribbon(stat="identity",alpha=.15) +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  scale_y_continuous(breaks = pretty_breaks(5)) + # Limit to not overpredict
  scale_color_colorblind() +
  facet_wrap(~Model,scales = "fixed") + 
  labs(x = expression(BC[EVI]), y = expression(BC[Biodiversity]),title = "") +
  theme(legend.position = "bottom") +
  theme(strip.text = element_text(size=12))
g2 <- results2 %>% dplyr::filter(Model != "All\n(195|4044)") %>% 
  ggplot(.,aes(x = avg.x, y = avg.fit, ymin= avg.fit - avg.se, ymax = avg.fit + avg.se)) + theme_few() +
  geom_line(lwd=.85) + geom_ribbon(stat="identity",alpha=.15) +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  scale_y_continuous(breaks = pretty_breaks(5)) + # Limit to not overpredict
  scale_color_colorblind() +
  facet_wrap(~Model,scales = "fixed") + 
  labs(x = expression(BC[EVI]), y = expression(BC[Biodiversity]),title = "") +
  theme(legend.position = "bottom") +
  theme(strip.text = element_text(size=12))
library(grid)
library(gridExtra)
library(cowplot)
pgrid <- plot_grid(g1 + theme(legend.position = "none"),
                   g2 + theme(legend.position = "none"),
                   ncol = 2,rel_widths = c(.4, .6),labels = NULL)
pgrid
cowplot::ggsave("Manuscript/SIFigure5.png",plot=pgrid)



# ----------------------------------------------------- #
#### Figure 3-4 - Relative Difference to Current as baseline ####
indic = "EVI2"#
dism = "BCVeg"#"SorVeg"  #Dissimilarity metric
metric =  "BC" #"_meanSPD"# "correctedBC"
res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds"))
out <- res.pred$Coeff

#  Histogram of predicted standard errors
qplot(se,eff,color=Time,data=out) + facet_wrap(~Model,scales = "free") + labs(x="Standard error",y="Coefficients",title="100 Models") 

# Average0
out$py[which(out$Time=="Present")] <- 0
out$py <- factor(out$py)
results <- out %>% mutate(py) %>% 
  dplyr::filter(Time != "PresentPast",py!=6) %>% 
  group_by(Model,py) %>%  # Group by year
  dplyr::summarise(avg.eff = mean(eff),
                   avg.se = mean(se),
                   avg.pval = mean(pval,na.rm = T),
                   N = max(N),
                   NS = max(N.sites))

# Z score loop - Very nasty programming
results2 <- data.frame()
for(i in unique(results$Model)){
  sub <- as.data.frame( subset(results, Model == i) )
  # Transform and calculate standard errors
  ref.eff = sub$avg.eff[which(sub$py == 0)]; ref.se = sub$avg.se[which(sub$py == 0)]
  for(k in 1:nrow(sub)){
    sub[k,"Z"] = (ref.eff - sub$avg.eff[k]) / sqrt( (ref.se^2) + (sub$avg.se[k]^2)  )
    sub[k,"Z.pval"] = 1 - (pnorm(sub[k,"Z"]))
  }
  results2 <- rbind(results2,sub)
}
results <- results2

#1 - pnorm( ( (0.32787 - 0.40491) / sqrt( (0.04612^2) + (0.05493^2)  ) ) )
#(0.50 / 0.28)

makeRatios <- function(results,pred="py",variable=NULL,exp=T,ci=F,seComp = F){
    df <- data.frame(Model=results$Model,py=results[[pred]],
                     avg.eff = results[["avg.eff"]],
                     avg.se = results[["avg.se"]],
                     avg.eff.low = (results[["avg.eff"]] - results[["avg.se"]]* ifelse(ci,1.96,1)),
                     avg.eff.high = (results[["avg.eff"]] + results[["avg.se"]]* ifelse(ci,1.96,1)),
                     N = results$N,
                     fix=NA,fix.low=NA,fix.high=NA)
    if(!is.null(variable)) { df$variable <- results[[variable]] }
    # Get Standard error
    if(exp==T) {
      df$avg.eff <- exp(df$avg.eff)
      df$avg.eff.low <- exp(df$avg.eff.low)
      df$avg.eff.high <- exp(df$avg.eff.high)
    } 
    # Correct errors for baseline percentage
    for(taxa in unique(df$Model) ){
      if(is.null(variable)){
        for(i in unique(results[[pred]]) ){    
          df$fix[which(df$Model==taxa & df$py==i)] <- ( df$avg.eff[which(df$Model==taxa & df$py==i)]  / df$avg.eff[which(df$Model==taxa & df$py=="0")] ) * 100
          if(seComp){
            df$fix.low[which(df$Model==taxa & df$py==i)] <- (df$avg.eff.low[which(df$Model==taxa & df$py==i)] / df$avg.eff.low[which(df$Model==taxa & df$py=="0")]) * 100
            df$fix.high[which(df$Model==taxa & df$py==i)] <- (df$avg.eff.high[which(df$Model==taxa & df$py==i)]  / df$avg.eff.high[which(df$Model==taxa & df$py=="0")]) * 100
          } else {
            df$fix.low[which(df$Model==taxa & df$py==i)] <- (df$avg.eff.low[which(df$Model==taxa & df$py==i)] / df$avg.eff[which(df$Model==taxa & df$py=="0")]) * 100
            df$fix.high[which(df$Model==taxa & df$py==i)] <- (df$avg.eff.high[which(df$Model==taxa & df$py==i)]  / df$avg.eff[which(df$Model==taxa & df$py=="0")]) * 100
          }
          } 
      } else {
         # Is a an extra lower level seperator (defined variable)
         for(l in unique(results[[variable]])){
           for(i in unique(results[[pred]]) ){    
             df$fix[which(df$Model==taxa & df$variable == l & df$py==i)] <- ( df$avg.eff[which(df$Model==taxa & df$variable == l & df$py==i)]  / df$avg.eff[which(df$Model==taxa & df$variable == l & df$py=="0")] ) * 100
             if(seComp){
               df$fix.low[which(df$Model==taxa & df$variable == l & df$py==i)] <- (df$avg.eff.low[which(df$Model==taxa & df$variable == l & df$py==i)] / df$avg.eff.low[which(df$Model==taxa & df$variable == l & df$py=="0")]) * 100
               df$fix.high[which(df$Model==taxa & df$variable == l & df$py==i)] <- (df$avg.eff.high[which(df$Model==taxa & df$variable == l & df$py==i)]  / df$avg.eff.high[which(df$Model==taxa & df$variable == l & df$py=="0")]) * 100
             } else {
               df$fix.low[which(df$Model==taxa & df$variable == l & df$py==i)] <- (df$avg.eff.low[which(df$Model==taxa & df$variable == l & df$py==i)] / df$avg.eff[which(df$Model==taxa & df$variable == l & df$py=="0")]) * 100
               df$fix.high[which(df$Model==taxa & df$variable == l & df$py==i)] <- (df$avg.eff.high[which(df$Model==taxa & df$variable == l & df$py==i)]  / df$avg.eff[which(df$Model==taxa & df$variable == l & df$py=="0")]) * 100
               }
           } 
         }
      }
    }
  return(df)
}

o <- makeRatios(results,exp=T,ci=F,seComp = T)
o$Time <- F
o$Time[which(o$py=="0")] <- T
o$Model <- factor(o$Model,
                         levels = c("All","Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))

o$py <- factor(o$py,labels = c(expression(yr[0]),expression(yr[1]),expression(yr[1:2]),expression(yr[1:3]),expression(yr[1:4]),expression(yr[1:5])))
g1 <- o  %>% dplyr::filter(Model == "All") %>% 
  ggplot(.,aes(x = py, y = fix-100, ymin = fix.low-100, ymax = fix.high-100)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  scale_x_discrete(labels = c(expression(yr[0]),expression(yr[1]),expression(yr[1:2]),expression(yr[1:3]),expression(yr[1:4]),expression(yr[1:5])) ) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_linerange(position = position_dodge(.75),size=1.5) + geom_point(size = 4,position = position_dodge(.75)) +
  scale_shape_circlefill() +
  #facet_wrap(~Model,scales="free") +
  theme(axis.ticks.x = element_blank() ) +
  labs(x = "Past period (years)", y = "Effect relative to current conditions (%)",title="") +
  guides(color = "none",shape="none")
g1
ggsave(filename = "Manuscript/Figure3.png",plot=g1)

# Transparent and increased font size for poster
# t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "white",colour = NA),
#            axis.text = element_text(colour="black",size = 22),axis.text.x =element_text(colour="black"),axis.text.y =element_text(colour="black"),axis.ticks = element_line(colour = "black"), axis.line = element_line(colour="black"),
#            axis.title = element_text(colour="black",size = 26),
#            strip.text.x = element_text(colour="black"),
#            legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
#            panel.grid.minor = element_blank(), title = element_text(colour="black")
# )
# ggsave(filename = "../../Conferences/110416_BES2016/OverallEffect.png",plot=g1+t,scale=1.1,dpi=400,bg = "transparent")


# Labels
lab = results %>% group_by(Model) %>% summarise(N = min(N),NS = min(NS))
lab$label <- factor(paste0(lab$Model,"\n","(",lab$N,"|",lab$NS,")"))
o$label <- lab$label[match(o$Model,lab$Model)]
o$label <- factor(o$label,levels = c("All\n(195|4044)","Plantae\n(41|775)","Invertebrates\n(102|1676)","Amphibia\n(9|142)","Reptilia\n(5|162)",
                                                   "Aves\n(23|969)","Mammalia\n(16|400)"))

g2 <- o  %>% dplyr::filter(Model != "All") %>% 
  ggplot(.,aes(x = py, y = fix-100, ymin = fix.low-100, ymax = fix.high-100)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = extended_breaks(5)) +
  scale_x_discrete(labels = c(expression(yr[0]),expression(yr[1]),expression(yr[1:2]),expression(yr[1:3]),expression(yr[1:4]),expression(yr[1:5])) ) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_linerange(position = position_dodge(.75)) + geom_point(size = 3,position = position_dodge(.75)) +
  facet_wrap(~label,scales="free_y") +
  theme(axis.text = element_text(color = "black",size=12),axis.ticks.x = element_blank() ) +
  labs(x = "Past period (years)", y = "Effect relative to current conditions (%)") +
  guides(shape="none",color = "none" )
g2
ggsave(filename = "Manuscript/Figure4.png",plot=g2)

# Transparent and increased font size for poster
# t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "white",colour = NA),
#            axis.text = element_text(colour="black",size = 22),axis.text.x =element_text(colour="black"),axis.text.y =element_text(colour="black"),axis.ticks = element_line(colour = "black"), axis.line = element_line(colour="black"),
#            axis.title = element_text(colour="black",size = 26),
#            strip.text.x = element_text(colour="black"),strip.background = element_rect(fill = "transparent",colour = NA),
#            legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
#            panel.grid.minor = element_blank(), title = element_text(colour="black")
# )
# ggsave(filename = "../../Conferences/110416_BES2016/IndividualEffect.png",plot=g2+t,scale=1.1,dpi=400,bg = "transparent")

# ------------------ #

## -- ##
#### Figure 5 - relative effect by Bodymass --- ####
# For Figure 5 we combine the effects overall with those of trophic level per group
res.pred <- readRDS("resSaves/outBCpermute_minStartEVI2_BC_BodyMassINT_OWN.rds")
out <- res.pred$Coeff
out$py[which(out$Time=="Sampling conditions")] <- 0
out$py <- factor(out$py)


results <- out %>% dplyr::filter(focus == "BM") %>% 
  #dplyr::filter(N>2) %>%  # FILTER out the marginal
  group_by(Model,py,variable) %>% 
  # First average per PermRun the similar ones
  dplyr::summarise(avg.eff = mean(estimate),
                   avg.se = mean(std.error),
                   avg.R2m = mean(R2m),
                   avg.R2c = mean(r2c),
                   N = max(N),
                   N.sites = max(N.sites)) %>% 
  ungroup()

# Z score loop
results2 <- data.frame()
for(i in unique(results$Model)){
  sub <- as.data.frame( subset(results, Model == i) )
  # Transform and calculate standard errors
  ref.eff = sub$avg.eff[which(sub$py == 0)]; ref.se = sub$avg.se[which(sub$py == 0)]
  sub$Z = (ref.eff - sub$avg.eff) / sqrt( (ref.se^2) + sub$avg.se^2  )
  sub$Z.pval = 1- (pnorm(sub$Z))
  results2 <- rbind(results2,sub)
}
results <- results2

results %>% dplyr::filter(variable == 1)

o <- makeRatios(results,variable = "variable",exp=T,ci=F,seComp = T)
o$Time <- F
o$Time[which(o$py=="0")] <- T
#o$variable <- factor(o$variable, levels = c(1,2,3,4), 
#                     labels = paste(c("0-2","2-20","20-100","100-5000000"),"(cm/g)"),
#                     ordered = T)
o$variable <- factor(o$variable, levels = c(1,2,3), 
                                          labels = c( (paste("S\n(0 - 9)")),
                                                      (paste("M\n(10 - 99)")),
                                                      (paste("L\n(>= 100)")) ),
                                          ordered = T)
#c( (10^(1)-1),(10^(2)-1),(10^(7)-1) )
o$Model <- as.character(o$Model)
o$Model[which(o$Model == "AllExceptPlantae")] <- c("Animalia (g)")
o$Model[which(o$Model == "Plantae")] <- c("Plants (cm)")
o$Model[which(o$Model == "Aves")] <- c("Birds (g)")
o$Model[which(o$Model == "Mammalia")] <- c("Mammals (g)")
o$Model[which(o$Model == "Reptilia")] <- c("Reptiles (g)")
o$Model[which(o$Model == "Aves")] <- c("Birds (g)")
o$Model[which(o$Model == "Amphibia")] <- c("Amphibians (g)")
o$Model[which(o$Model == "Invertebrates")] <- c("Invertebrates (g)")
o$Model <- factor(o$Model,levels = c("All","Animalia (g)","Plants (cm)","Invertebrates (g)","Amphibians (g)","Reptiles (g)","Birds (g)","Mammals (g)"))

## Accross all taxonomic groups, but removing the reptiles
g2 <- o %>% dplyr::filter(Model %in% c("Plants (cm)","Birds (g)","Mammals (g)") ) %>% 
  ggplot(.,aes(x = py, y = fix-100, ymin = fix.low-100, ymax = fix.high-100,
               group=Model,color=Model)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  scale_x_discrete(labels = c(expression(yr[0]),expression(yr[1]),expression(yr[1:2]),expression(yr[1:3]),expression(yr[1:4]),expression(yr[1:5])) ) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_linerange(position = position_dodge(.75)) + geom_point(size = 3,position = position_dodge(.75)) +
  scale_shape_circlefill() + 
  scale_color_manual(values = c("#3B5C00","#5D959A","#D4605A")) +
  facet_wrap(~variable,scales="free_y")+
  theme(axis.text = element_text(size = 10),axis.ticks.x = element_blank() ) +
  labs(x = "Past period (years)", y = "Effect relative to current conditions (%)") +
  guides(shape="none",color = guide_legend(title="",nrow = 1) ) +
  theme(legend.position = "bottom")
g2
# The number of studies
o$lab <- NA
o$lab[which(o$variable=="S\n(0 - 9)")] <- "S"
o$lab[which(o$variable=="M\n(10 - 99)")] <- "M"
o$lab[which(o$variable=="L\n(>= 100)")] <- "L"
o$lab <- factor(o$lab,levels=c("S","M","L"))
g3 <- o %>% dplyr::filter(Model %in% c("Plants (cm)","Birds (g)","Mammals (g)") ) %>% 
  dplyr::filter(py == 5) %>% 
  ggplot(.,aes(x = lab,y=N, Group=Model,fill=Model)) + theme_few(base_size = 18) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks = pretty_breaks(10)) + 
  scale_fill_manual(values = c("#3B5C00","#5D959A","#D4605A")) +
  labs(x = "") +
  guides(shape="none",fill = guide_legend(title="Body size\n (cm|g)",nrow = 2) ) +
  theme(legend.position = "none",axis.text.x = element_text(angle=0,vjust = .5))
g3

library(grid)
library(gridExtra)
library(cowplot)
grobs <- ggplotGrob(g2 + theme(legend.position="none",legend.direction = "vertical") )$grobs

pgrid <- plot_grid(g2 + theme(legend.position = "none"),
                   g3 + theme(legend.position = "none") + labs(y="N",x=""),
                   ncol = 2,rel_widths = c(.8, .2),labels = NULL)
pgrid
cowplot::ggsave("Manuscript/Figure6.png",plot=pgrid,dpi = 400)


# ---- #
# Individual figure just for the poster #
lab = results %>% dplyr::filter(Model %in% c("Mammalia") ) %>%
  group_by(variable) %>% summarise(N = min(N),NS = min(N.sites))

lab$label <- factor(paste0(lab$Model,"\n","(",lab$N,"|",lab$NS,")"))
o$label <- lab$label[match(o$Model,lab$Model)]
o$label <- factor(o$label,levels = c("All\n(187|3952)","Plantae\n(37|730)","Invertebrates\n(98|1605)","Amphibia\n(8|127)","Reptilia\n(4|138)",
                                     "Aves\n(21|911)","Mammalia\n(15|372)"))

g <- o %>% dplyr::filter(Model %in% c("Mammalia (g)") ) %>% 
  ggplot(.,aes(x = py, y = fix-100, ymin = fix.low-100, ymax = fix.high-100,
               group=Model,color=Model)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_linerange(position = position_dodge(.75)) + geom_point(size = 3,position = position_dodge(.75)) +
  scale_shape_circlefill() + 
  #scale_color_manual(values = c("#3B5C00","#5D959A","#D4605A")) +
  facet_wrap(~variable,scales="free_y")+
  theme(axis.text = element_text(color = "#808080"),axis.ticks.x = element_blank() ) +
  labs(x = "Past period (years)", y = "Effect relative to current conditions\n (%)") +
  guides(shape="none",color = guide_legend(title="",nrow = 1) ) +
  theme(legend.position = "none")
g

# t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
#            axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
#            strip.text.x = element_text(colour="white"),         
#            legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
#            panel.grid.minor = element_blank(), title = element_text(colour="white")
# )
# 
# ggsave(filename="../../Conferences/110416_BES2016/MammalBodmass_trans.png",plot=g+t,scale=1.1,dpi=400,bg = "transparent")
# 
# -------------------------------------------------- #
# -------------------------------------------------- #

# Overall
results <- out %>% dplyr::filter(focus == "BM") %>% 
  #dplyr::filter(Model %in% c("AllExceptPlantae","Plantae")) %>%  # FILTER out the marginal
  group_by(Model,py,variable) %>% 
  # First average per PermRun the similar ones
  dplyr::summarise(avg.eff = mean(estimate),
                   avg.se = mean(std.error),
                   N = max(N),
                   N.sites = max(N.sites)) %>% 
  ungroup()
o <- makeRatios(results,variable = "variable",exp=T,ci=F)
o$Time <- F
o$Time[which(o$py=="0")] <- T
o$Model <- factor(o$Model,
                  levels = c("All","AllExceptPlantae","Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))
o$variable <- droplevels(o$variable)
o$variable <- factor(o$variable, levels = c(1,2,3), 
                     labels = c( (paste("S\n(0 - 9)")),
                                 (paste("M\n(10 - 99)")),
                                 (paste("L\n(100 - 9999999)")) ),
                     ordered = T)

lab = results %>% dplyr::filter(Model %in% c("AllExceptPlantae","Plantae")) %>%
  group_by(variable) %>% summarise(N = min(N),NS = min(N.sites))

lab$label <- factor(paste0(lab$Model,"\n","(",lab$N,"|",lab$NS,")"))
o$label <- lab$label[match(o$Model,lab$Model)]
o$label <- factor(o$label,levels = c("All\n(187|3952)","Plantae\n(37|730)","Invertebrates\n(98|1605)","Amphibia\n(8|127)","Reptilia\n(4|138)",
                                     "Aves\n(21|911)","Mammalia\n(15|372)"))


g <- o %>% 
  ggplot(.,aes(x = py, y = fix-100, ymin = fix.low-100, ymax = fix.high-100,
               group=Model,color=Model)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_linerange(position = position_dodge(.75)) + geom_point(size = 3,position = position_dodge(.75)) +
  scale_color_colorblind() + 
  #scale_color_manual(values = c("#3B5C00","#5D959A","#D4605A")) +
  facet_wrap(~variable)+
  theme(axis.text = element_text(color = "#808080"),axis.ticks.x = element_blank() ) +
  labs(x = "Past period (years)", y = "Effect relative to current conditions (%)") +
  guides(shape="none",color = guide_legend(title="",nrow = 1) ) +
  theme(legend.position = "none")
g

# t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
#            axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
#            strip.text.x = element_text(colour="white"),         
#            legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
#            panel.grid.minor = element_blank(), title = element_text(colour="white")
# )
# 
# ggsave(filename="../../Conferences/110416_BES2016/MammalBodmass_trans.png",plot=g+t,scale=1.1,dpi=400,bg = "transparent")
# 


# ------------------- #
#### Figure 6 - Relative effect - Trophic ####
results <- out %>% dplyr::filter(focus == "Trophic") %>% 
  group_by(Model,py,variable) %>% 
  # First average per PermRun the similar ones
  dplyr::summarise(avg.eff = mean(estimate),
                   avg.se = mean(std.error),
                   N = max(N),
                   NS = max(N.sites)) %>% 
  ungroup()

# Z score loop
results2 <- data.frame()
for(i in unique(results$Model)){
  sub <- as.data.frame( subset(results, Model == i) )
  # Transform and calculate standard errors
  ref.eff = sub$avg.eff[which(sub$py == 0)]; ref.se = sub$avg.se[which(sub$py == 0)]
  sub$Z = (ref.eff - sub$avg.eff) / sqrt( (ref.se^2) + sub$avg.se^2  )
  sub$Z.pval = (pnorm(sub$Z))
  results2 <- rbind(results2,sub)
}
results <- results2

o <- makeRatios(results,variable = "variable",exp=T,ci=F,seComp = T)
o$Time <- F
o$Time[which(o$py=="0")] <- T
o$Model <- factor(o$Model,
                  levels = c("All","AllExceptPlantae","Plantae","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))
o$variable <- droplevels(o$variable)
o$variable <- factor(o$variable, levels = c("Autotroph","Herbivore","Omnivore","Carnivore","Detritivore"))
# Make a label
lab = results %>% dplyr::filter(py==5,Model == "All") %>% 
  group_by(variable) %>% summarise(N = min(N),NS = min(NS))
lab$label <- factor(paste0(lab$variable," ","(",lab$N,"|",lab$NS,")"))
o$label <- lab$label[match(o$variable,lab$variable)]
o$label <- factor(o$label,levels = c("Autotroph (42|775)","Herbivore (61|1133)","Omnivore (31|1034)","Carnivore (15|329)","Detritivore (20|357)") )

g2 <- o %>% dplyr::filter(Model == "All") %>% 
  ggplot(.,aes(x = py, y = fix-100, ymin = fix.low-100, ymax = fix.high-100,
               group=label,color=label)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  scale_x_discrete(labels = c(expression(yr[0]),expression(yr[1]),expression(yr[1:2]),expression(yr[1:3]),expression(yr[1:4]),expression(yr[1:5])) ) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_linerange(position = position_dodge(.4),size=1.25) + geom_point(size = 3,position = position_dodge(.4)) +
  scale_shape_circlefill() + 
#  facet_wrap(~Model,scales = "free_y") +
  scale_color_manual(values =c("#7fbf7b","darkgreen","#ffc425","tomato2","tan4")) +
  theme(axis.text = element_text(color = "black"),axis.ticks.x = element_blank() ) +
  labs(x = "Past period (years)", y = "Effect relative to current conditions (%)") +
  guides(shape="none",color = guide_legend(title="",nrow =2) ) +
  theme(legend.position = "none")
g2
ggsave("Manuscript/Figure7.png",plot=g2 )

# # Transparent for poster
# t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "white",colour = NA),
#            axis.text = element_text(colour="black",size = 22),axis.text.x =element_text(colour="black"),axis.text.y =element_text(colour="black"),axis.ticks = element_line(colour = "black"), axis.line = element_line(colour="black"),
#            axis.title = element_text(colour="black",size = 26),
#            strip.text.x = element_text(colour="black"),strip.background = element_rect(fill = "transparent",colour = NA),
#            legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="black",size=20,hjust = 0),
#            panel.grid.minor = element_blank(), title = element_text(colour="black")
# )
# ggplot2::ggsave("../../Conferences/110416_BES2016/TrophicLevels.png",
#        plot=g2+t+labs(y=""),scale=1.1,dpi=400,bg = "transparent"
#        )



#### Table-1 of R square values per time and group ####
# BCVeg
res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds"))
#res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5yINDIVIDUAL",indic,"_",metric,"_",dism,".rds"))
out <- res.pred$Coeff

# Average
out$py[which(out$Time=="Present")] <- 0
out$py <- factor(out$py)
results <- out %>% dplyr::filter(Time != "PresentPast",py!=6) %>% 
  dplyr::filter(Model == "All") %>% 
  group_by(Model,PermRun,py) %>%  # Group by model and time
  dplyr::summarise(avg.eff = mean(eff),
                   avg.se = mean(se),
                   avg.R2m = mean(R2m),
                   avg.R2c = mean(r2c),
                   avg.pval = mean(pval)#length(which(pval < 0.05)),
                   )
# Maximal number of Study and N.sites included
out %>% dplyr::group_by(Model) %>% summarise(N = max(N), N.sites = max(N.sites))

# Because of permutations, some do vary
out <- res.pred$StudyCoeff
out %>% dplyr::filter(py == 1) %>% summarise(N = n_distinct(SS))
out %>% dplyr::filter(py == 5) %>% summarise(N = n_distinct(SS))
# Which studies
SS1 <- out %>% dplyr::filter(py == 1) %>% dplyr::select(SS) %>% mutate(SS = as.character(SS)) %>% distinct()
SS2 <- out %>% dplyr::filter(py == 5) %>% dplyr::select(SS) %>% mutate(SS = as.character(SS)) %>% distinct()
marfunky::setdiff2(SS1$SS,SS2$SS)


# Difference in effect per group
#o <- makeRatios(results,pred = "py",exp=T,ci=F)

# Relative difference in R2m between sampling conditions and past conditions
results$rel.changeR2 <- (results$avg.R2m - results$avg.R2m[which(results$py==0)]) / results$avg.R2m[which(results$py==0)]

# Output Table
# Overall | TimePeriod | Avg.Effect | Avg.SE | Avg.R2m | Avg.R2C | Relative Diff Present-Past  
results <- results %>% 
  mutate(N = 187,N.sites = 3955) %>% # Sampling studies and sites of present as ref
  #dplyr::filter(Time == "Past") %>% 
  # Clean up and reorder
  dplyr::select(Model,py,avg.eff,avg.se,avg.R2m,avg.R2c,rel.changeR2,N,N.sites) %>% 
  rename(TimePeriod = py) %>% ungroup()
# Round numeric columns
results[,c(4:8)] <- apply(results[,c(4:8)],2,function(x) round(x,3))
# Convert relative change to percentage
results$rel.changeR2 <- results$rel.changeR2 * 100

results$Model <- NULL # No need for Model
write.csv(results,"Manuscript/Table1.csv",sep=",",dec=".",row.names=F)

# ---- #

#### NOT-INCLUDED ~ SI Table-2 Effecter per Terrestial Biome | Studywise ####
library(tidyr)
# BCVeg
res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds"))
out <- res.pred$StudyCoeff
# Calculate average Sampling conditions / Past conditions
results <- out %>% dplyr::filter(TGrouping == "All",y==5) %>%  # Overall for 5 year period
  dplyr::group_by(SS,time) %>%
  dplyr::summarise(avg.value = mean(value,na.rm=T),
                   avg.se = mean(value.se,na.rm=T)
                   ) %>% ungroup() %>% 
  mutate(time = str_replace_all(time," ",".")) %>% # Correcting time
  melt(.) %>% spread(time,value) # Melt & Spread

# Calculate difference between past and sampling
results <- results %>% dplyr::group_by(SS) %>% 
  dplyr::summarise(avg.effect = (Past.conditions[which(variable == "avg.value")] / Sampling.conditions[which(variable == "avg.value")]) *100,
                   avg.se = (Past.conditions[which(variable == "avg.se")] /
                               Sampling.conditions[which(variable == "avg.se")])*100
                   #avg.effect.high = ((Past.conditions[which(variable == "avg.value")] + Past.conditions[which(variable == "avg.se")]) /
                   #                  (Sampling.conditions[which(variable == "avg.value")] + Sampling.conditions[which(variable == "avg.se")]))*100,
                   #avg.effect.low = ((Past.conditions[which(variable == "avg.value")] - Past.conditions[which(variable == "avg.se")]) /
                   #                      (Sampling.conditions[which(variable == "avg.value")] - Sampling.conditions[which(variable == "avg.se")]))*100
                   ) 
# Get the majority biome of sites per study
ag <- aggregate(list(Biome = sites$Biome),list(SS = sites$SS),function(x) names(which.max(table(x))) )
results <- left_join(results,ag)

# Table for effects per biome in 5 year period
# Biome | N | Avg.Effect | diff.se + (Significant difference of interaction ?)
o <- results %>% dplyr::group_by(Biome) %>% 
  dplyr::summarize(N = n(),
                   avg.effect = mean(avg.effect,na.rm=T),
                   avg.se = mean(avg.se,na.rm=T)
                   
#                   avg.effect.high = mean(avg.effect.high,na.rm=T),
#                   avg.effect.low = mean(avg.effect.low,na.rm=T)
                   )
o$avg.effect <- o$avg.effect - 100 
o$avg.se <- o$avg.se - 100 
write.csv(o,"Manuscript/Table2.csv",row.names=F)

# ---- #
#### SI FIGURE 4 - Potential Biases in study-wise coefficients ####
# Parameter above at Figure 3
res.pred <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds"))
out <- res.pred$StudyCoeff
out <- out %>% dplyr::group_by(time,SS,y) %>% 
  dplyr::summarise(avg.dbc = mean(value),
            avg.dbc.se = mean(value.se)) %>% ungroup() %>% 
  mutate(time = str_replace_all(time," ",".")) 

out <- reshape2::melt(out,id.vars=c("SS","time","y")) %>% 
        tidyr::spread(time,value) # Split
# Full sites
sites <- readRDS("AllPredictsSites_full.rds")# Get sites
sites <- subset(sites,SS %in% out$SS)

# Join in Grouping
out <- left_join(out,( sites %>% dplyr::select(SS,TGrouping) %>% distinct() ))

## Potential biases ##
# Per startyear
o <- left_join(out,sites %>% group_by(SS) %>% summarise(st = factor(max(startyear))) ) %>% 
  dplyr::filter(variable == "avg.dbc",y== 5) %>%  mutate(Effect = ((Past.conditions / Sampling.conditions) *100)-100)
# Create labs
la <- as.data.frame(table(o$st)) %>% rename(st = Var1)
g1 <- qplot(st,Effect,geom="boxplot",data=o) + theme_few() +
  scale_y_continuous(breaks = pretty_breaks(7),limits = c(-50,50)) +
  labs(x = "Year of sampling (per study)", y = "Overall effect relative\n to current conditions\n (%)") +
  geom_label(data = la,aes(x =st, y = 49,label=Freq),colour = "black", fontface = "bold",inherit.aes = F)
g1

# Max Lin Extent
o <- left_join(out,sites %>% group_by(SS) %>% summarise(ml = as.numeric(mean(Max_linear_extent))) ) %>%
  dplyr::filter(variable == "avg.dbc",y== 5) %>% 
  mutate(Effect = ((Past.conditions / Sampling.conditions) *100)-100)
g2 <- qplot(ml,Effect,data=o) + theme_few() + 
  scale_y_continuous(breaks = pretty_breaks(7),limits = c(-50,50)) +
  scale_x_continuous(breaks= c(0.1,1,10,100,1000),labels = c(0.1,1,10,100,1000),trans = "log10") +
  labs(x = "Average maximum linear extent per study\n(in m, log10-transf.)", y = "Overall effect relative\n to current conditions\n (%)")
g2

# Per study length
o <- left_join(out,sites %>% group_by(SS) %>% summarise(sl = as.numeric(mean(StudyLength))) ) %>% 
  dplyr::filter(variable == "avg.dbc",y== 5) %>% 
  mutate(Effect = ((Past.conditions / Sampling.conditions) *100)-100)
g3 <- qplot(sl,Effect,data=o) + theme_few() + 
  scale_y_continuous(breaks = pretty_breaks(7),limits = c(-50,50)) +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  labs(x = "Average study sampling duration (in days)", y = "Overall effect relative\n to current conditions\n (%)")
g3

# Per Latitude center
o <- left_join(out,sites %>% group_by(SS) %>% summarise(lat = as.numeric(mean(Latitude))) )%>%
  dplyr::filter(variable == "avg.dbc",y== 5) %>% 
  mutate(Effect = ((Past.conditions / Sampling.conditions) *100)-100)
g4 <- qplot(lat,Effect,data=o) + theme_few(base_family = "Arial") + 
  geom_vline(xintercept = 0,linetype="dotted") +
  # And tropics
  geom_rect(aes(xmin =-23.5,xmax=23.5,ymin=-Inf,ymax=Inf),fill="grey90",alpha=.1)+
  geom_vline(xintercept = 0,linetype="dotted",lwd=1) +
  geom_point() +
  scale_y_continuous(breaks = pretty_breaks(7),limits = c(-50,50)) +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  labs(x = "Average latitude of sites within study (°)", y = "Overall effect relative\n to current conditions\n (%)")
g4

# Now combine them all
library(cowplot)
pg <- plot_grid(g1,g2 + labs(y=""),
                g3,g4 + labs(y=""),
                nrow = 2,labels = "auto")
pg
cowplot::ggsave("Manuscript/SIFigure4.png",plot=pg)




#### SI FIGURE SAR correction - Are observed patterns spatiall autocorrelated?  ####
library(lme4)
library(nlme)
library(arm)
library(pbkrtest)
library(MuMIn)
library(gamm4)
### Construct model and visualize results ###
source("000_HelperFunction.R")
indic = "EVI2"#"EVI2"
dism = "BCVeg" # "SorVeg" #Dissimilarity metric
metric =  "BC" #"_meanSPD"# "correctedBC"
transf = FALSE # Should response be transformed?
ylist <- c(1,2,3,4,5) # Sequence of years

# Get list of 20 permutations (one fifth of the dataset)
f.permut <- sort(list.files(PermPath,pattern = paste0("s_metric",dism,"_permute"),full.names = T))
res.cor <- list() #Output

# Make sure that res.dist exists
stopifnot(exists("res.dist"))

# Sequence of historic years
for(y in ylist){
  myLog("Processing past period - ",y)
  # Load matrices
  res_yearsampling <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"samplingperiod_0y-matrix.rds")) # The year of before start
  res_yearbefore <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"pastperiod_",y,"y-matrix.rds")) # The years before
  for(p in seq(length(f.permut))){
    myLog("Running permutation ",p)
    # Load permutation
    s.metric <- readRDS(f.permut[p])
    
    # Available past data
    focal.period  = ymd("2000.02.18") + (365*(y+1)) # Available past in total
    # Add in startdate
    s.metric$startdate_x <- sites$Sample_start_earliest[match(s.metric$SSBS_x,sites$SSBS)]
    s.metric$startdate_y <- sites$Sample_start_earliest[match(s.metric$SSBS_y,sites$SSBS)]
    # Calculate distance in days and remove sites with greater dist of 3 months (~90days)
    s.metric <- mutate(s.metric,startd = as.numeric( abs(startdate_x - startdate_y))) %>% 
      dplyr::filter(startd <= 90) %>% # 3 months (~ 90 days)
      dplyr::filter(startdate_x > focal.period & startdate_y > focal.period) %>%  # Must have sufficient historical coverage
      # Filter to only those studies that are within the maximum of ylist, thus equalizing nr of studies for all years
      dplyr::filter(startdate_x > (ymd("2000.02.18") + (365*(max(ylist)+1))) & startdate_y > (ymd("2000.02.18") + (365*(max(ylist)+1)))) %>% 
      dplyr::select(-startdate_x,-startdate_y,-startd) # remove again
    
    # Subset and assemble the dataset of the permutations and full matrix
    per.df <- data.frame()
    for(study in unique( s.metric$SS) ){
      sub <- subset(s.metric,SS == study)
      if(!(study %in% names(res_yearbefore))) next()
      # Melt past and sampling cond
      
      m <- reshape2::melt( res_yearsampling[[as.character(study)]] ) %>% rename(SSBS_x = Var1,SSBS_y = Var2,EVIpres = value)
      m2 <- reshape2::melt( res_yearbefore[[as.character(study)]] ) %>% rename(SSBS_x = Var1,SSBS_y = Var2,EVIpast = value)
      mm <- merge.data.frame(m2,m,by =c("SSBS_x","SSBS_y"),all.x =T) # Merge present with past
      rm(m,m2)
      # Merge with permutation to get only the permuteted pairs
      per.df <- rbind(per.df,
                      merge.data.frame(sub,mm,by =c("SSBS_x","SSBS_y"),all.x = T)
      )
    }
    stopifnot(nrow(s.metric)==nrow(per.df)) # Check. Should be equal to the permutation dataset
    
    per.df$py <- y # Add year
    # -------------------------------------------- #
    myLog("Running permutation ",p," --- ","Built model sets")
    # Formatting
    per.df$TGrouping <- sites$TGrouping[match(per.df$SS,sites$SS)] # Get taxonomic grouping
    # Exclude Fungi
    per.df <- subset(per.df,TGrouping != "Fungi")
    # Get complete subset
    per.df <- subset(per.df,complete.cases(per.df))
    if(exists("res.dist")){
      per.df <- merge.data.frame(per.df,res.dist,c("SS","SSBS_x","SSBS_y"),all.x=T) %>% mutate(log10Distance = log10(Distance))
    }
    
    myLog("Fitting models")
    options(na.action = na.omit)
    mlist <- list()
    ## Need to correct to dissimilarity if SorVeg
    if(dism == "SorVeg") per.df$value <- abs(1 - per.df$value)
    # Indiv coefficients
    print("All")
    if(transf) per.df$value <- asin(sqrt(per.df$value))
    
    # Normal models
    mod_all_pr <- glmer(value ~ EVIpres + (1+EVIpres|SS), data=per.df,family=gaussian)
    mod_all_pa <- glmer(value ~ EVIpast + (1+EVIpast|SS), data=per.df,family=gaussian)
    
    # Save residuals
    per.df$pres.resid <- resid(mod_all_pr)
    per.df$past.resid <- resid(mod_all_pa)
    
    res.sub.cor <- data.frame()
    # Now loop through
    for(ss in unique(per.df$SS)){
      # Create sub.data
      sub.data <- subset(per.df,SS == as.character(ss))
      
      # Save output
      cp = try(cor.test(sub.data$log10Distance,sub.data$pres.resid),silent=T)
      if(class(cp)!="try-error"){      
        # Get tidy out put of pearson product moment - PRESENT
        x = tidy(cp) %>% 
          mutate(p = p,py=y,SS=as.character(ss),Time="Present",
                 meanDist = mean(sub.data$Distance,na.rm=T),
                 meanLDist = mean(sub.data$log10Distance,na.rm=T),
                 meandBC = mean(sub.data$EVIpres,na.rm=T))
        # If degrees of freedom equal to one = Insert NA in conf
        if(!("conf.low" %in% names(x))) {
          x$conf.low <- NA
          x$conf.high <- NA
        }
        res.sub.cor <- rbind(res.sub.cor,x)
        
        }
      cp = try(cor.test(sub.data$log10Distance,sub.data$past.resid),silent=T)
      if(class(cp)!="try-error"){      
        # Get tidy out put of pearson product moment - PASR
        x = tidy(cp) %>% 
          mutate(p = p,py=y,SS=as.character(ss),Time="Past",
                 meanDist = mean(sub.data$Distance,na.rm=T),
                 meanLDist = mean(sub.data$log10Distance,na.rm=T),
                 meandBC = mean(sub.data$EVIpast,na.rm=T))
        # If degrees of freedom equal to one = Insert NA in conf
        if(!("conf.low" %in% names(x))) {
          x$conf.low <- NA
          x$conf.high <- NA
        }
        res.sub.cor <- rbind(res.sub.cor,x)
      }
    }
    # Save
    res.cor[["Tests"]] <- rbind(res.cor[["Tests"]], res.sub.cor )
    
    # How many are significant?
    #rm.pres <- res.sub.cor$SS[which(res.sub.cor$p.value < 0.05 & res.sub.cor$Time == "Present")]
    rm.pres <- "DL1_2012__Hernandez 2"
    myLog("Found ",length(rm.pres), " studies with significant autocorrelation. PRESENT")
    # Refit by removing the said studies
    mod_all_pr2 <- glmer(value ~ EVIpres + (1+EVIpres|SS), data=per.df,subset = !(SS%in%rm.pres) ,family=gaussian)
    rm(rm.pres)
    
    #rm.past <- res.sub.cor$SS[which(res.sub.cor$p.value < 0.05 & res.sub.cor$Time == "Past")]
    rm.past <- "DL1_2012__Hernandez 2"
    myLog("Found ",length(rm.past), " studies with significant autocorrelation. PAST")
    # Refit by removing the said studies
    mod_all_pa2 <- glmer(value ~ EVIpast + (1+EVIpast|SS), data=per.df,subset = !(SS%in%rm.past) ,family=gaussian)
    
    
    out <- rbind( 
      tidy(mod_all_pr)[2,] %>% mutate(p=p,py=y,Time="Present",Type="Normal"),
      tidy(mod_all_pr2)[2,] %>% mutate(p=p,py=y,Time="Present",Type="Corrected"),
      tidy(mod_all_pa)[2,] %>% mutate(p=p,py=y,Time="Past",Type="Normal"),
      tidy(mod_all_pa2)[2,] %>% mutate(p=p,py=y,Time="Past",Type="Corrected")
      )
    
    # Save the different effects as well
    res.cor[["Coeff"]] <- rbind(res.cor[["Coeff"]], out )
    rm(out)
    
  }
}
saveRDS(res.cor,paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,"_StudyCorrelations.rds")) # dBC both axis

res.cor <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,"_StudyCorrelations.rds"))$Tests

# Average accross the permutations
out <- res.cor %>% 
  group_by(py,SS,Time) %>% 
  summarise(avg.est = mean(estimate),
            avg.pval = mean(p.value),
            avg.dBC = mean(meandBC))
res.cor %>% group_by(py,Time) %>% summarise(N = n_distinct(SS))
# For p-values
# In order to compensate for this, you must divide your threshold α by the number of tests n,
# so a result is significant when p < α/n 

## Make a cor plot assessment ##
# How many are significant and which are those #
length(which(out$avg.pval<0.05))
out$SS[which(out$avg.pval<0.05)]
# Hernandez seems problematic - "DL1_2012__Hernandez 2"

# Overview figure showing effect with studies having sig. autocorrelation removed
out <- readRDS(paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,"_StudyCorrelations.rds"))$Coeff

out$py[which(out$Time=="Present")] <- 0
out$py <- factor(out$py)
results <- out %>% 
  group_by(py,Time,Type) %>% 
  summarise(avg.eff = mean(estimate),
            avg.se = mean(std.error)) %>% mutate(Model="All",N=NA)

o <- makeRatios(results,variable="Type",exp=T,ci=F)
o$Time <- F
o$Time[which(o$py=="0")] <- T

g1 <- o  %>% 
  ggplot(.,aes(x = py, y = fix-100, ymin = fix.low-100, ymax = fix.high-100,shape=Time,group=variable,color=variable)) +
  theme_few(base_size = 18) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_linerange(position = position_dodge(.75)) + geom_point(size = 3,position = position_dodge(.75)) +
  scale_shape_circlefill() +
  scale_fill_aaas() +
  theme(axis.text = element_text(color = "#808080"),axis.ticks.x = element_blank() ) +
  labs(x = "Past period (years)", y = "Effect relative to current conditions\n (%)") +
  guides(shape="none", color = guide_legend(title = "Spatial\nautocorrelation"))
g1
ggsave("Manuscript/SI-FigureXX_SAR2.png",plot=g1)


# ------------------------ #
#### NOT INCLUDED - SI Figure - All time series overlayed ####
library(zoo)
# Idea: Overlay all time series and calculate average per taxonomic group
# Prop: Need to get rid of time
refy <- readRDS("Center_ReferenceYearMonthlyFilled.rds")
specPast <- readRDS("Center_FullTimeSeriesMonthlyFilled.rds")

res <- data.frame()
# Looop
for(site in names(refy)){
  myLog(site)
  x <- specPast[[site]]$EVI2
  y <- refy[[site]]$EVI2
  z <- rbind(x,y)
  
  # Create ideal time series
  min.start <- (unique( year(y) ) - 5 ) # subset to values 5 years before
  alldates <- seq.Date(as.Date(yearmon(min.start)),as.Date(as.yearmon(paste("Dec",unique(year(y))) )), by="months")
  ideal <- zoo(-9999,order.by = as.yearmon(alldates))
  z <- merge(ideal,z,all=T)
  z <- z[which(z==-9999)] # Subset to ideal time series
  
  res <- rbind(res,
               data.frame(SSBS = site,months = month(z[,2]),values = coredata(z[,2]),nr = seq(length(z[,2])),
                          time = c(rep("past",length(z[,2])-length(y)),rep("present",length(y)))
                          )
  )
}
  

# Now merge back with sites
sites <- readRDS("AllPredictsSites_full.rds") %>% 
  dplyr::select(SS,SSBS,midyear,TGrouping)
res$TGrouping <- sites$TGrouping[match(res$SSBS,sites$SSBS)]
res <- res[which(res$SSBS %in% sites$SSBS[which(sites$midyear>2005)]),]


# Construct average per taxonomic group and time
res.avg <- res %>% group_by(time,TGrouping,nr) %>% 
  dplyr::summarise(val.avg = mean(values,na.rm=T),
                   val.median = median(values,na.rm=T),
                   val.sd = sd(values,na.rm=T),
                   val.var = var(values,na.rm = T),
                   val.cv = marfunky::co.var(values,na.rm = T)
                   )
res.avg$TGrouping <- factor(res.avg$TGrouping,
                           levels = c("All","Plantae","Fungi","Invertebrates","Amphibia","Reptilia","Aves","Mammalia"))

# Plot
# Grey alpha lines per group(site)
# solid line and sd with average per taxonomic group
g <- ggplot(res,aes(x = nr, y = values, group=SSBS)) + theme_bw()
g <- g + geom_vline(xintercept = c(12,24,36,48),linetype="dashed",alpha=.6,size=.95) + # Past 
  geom_vline(xintercept = 60, linetype="solid",alpha=.6,size=.95) # present
g <- g + geom_line(color="lightgrey",alpha=.25,size = .25)
# Add average per taxonomic group
g <- g + geom_line(data=res.avg,aes(x = nr,y = val.median,group=TGrouping,color=TGrouping),size=1.5) + 
  scale_color_brewer(palette="Spectral")
# Add some labels on top
g <- g + annotate("text",x = seq(6,72,by = 12),y=1,size=7,label = c("y5","y4","yr3","y2","y1","yrs"))
g <- g + scale_y_continuous(breaks = pretty_breaks(7)) + scale_x_continuous(expand=c(0,1),labels = NULL, breaks = NULL)
g <- g + theme(legend.position = "bottom") + guides(color = guide_legend("Taxonomic group"))
g <- g + labs(x = "", y = "Enhanced vegetation index (EVI2)")
ggsave(paste0( "PubFigures/SI_AllTSOverlayed_EVI.png"),plot = g,width = 12,height = 10)


#### SI Figure - Temporal spread ####
# GGplot version
f.permut <- readRDS( list.files(".",pattern = "s_metricBC_permute")[1] ) %>% 
  dplyr::select(SS) %>% distinct()
s <- left_join(f.permut,sites,by="SS")

g <- ggplot(s, aes(x=midyear)) + theme_bw() +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=1,
                 colour="white", fill="black") +
#  geom_density(alpha=.2, fill="#FF6666") + # Overlay with transparent density plot
  scale_x_continuous(breaks = pretty_breaks(n=20),limits=c(2000,2014)) +
  scale_y_continuous(expand=c(0,0),breaks= pretty_breaks(10),limits=c(0,.15)) +
  labs(x = "Year of sampling",y = "Density of studies")
g

g2 <- ggplot(s,aes(x=Max_linear_extent)) + theme_bw() +
  geom_histogram(aes(y=..density..),  # Histogram with density instead of count on y-axis
                                                              binwidth=1,
                                                              colour="black", fill="black") +
  scale_x_continuous(limits=c(0,4000)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Maximum linear extent", y = "")
g2

cp = cowplot::plot_grid(g,g2,nrow=1,rel_widths = c(1.5,1),labels = "AUTO")
cowplot::save_plot(paste0(OutputPath,"SI-Figure3_TemporalSpatialDensityOfStudies.png"),base_width = 12,plot=cp)


#----------------------------------------------#
#### SI Figure - Effectv of distance on residuals ####
# Load before 
# res.dist, per.df
# model and residuals
g <- qplot(resid(mod_all_pa),x=log10(per.df$Distance)) + theme_few() +
  scale_x_continuous(breaks=pretty_breaks(10)) + scale_y_continuous(breaks=pretty_breaks(10)) +
  labs(x = substitute(paste("Pairwise ",log10[Distance])), y = substitute(paste("Model residuals  ",dBC[Current])))
g
ggsave(filename = "PubFigures/SI-Figure4-ResidualsPairwiseDistance.png",plot = g)

#### Supplementary Figure 3 ####
# Recreate figure 3 for appendix.
library(ggplot2);library(scales);library(ggthemes)
bc <- function(x,y){  ( sum( abs(x - y)) / (sum(x)+sum(y)) ) }
#
# 46 observations per year
# 6 years

res <- data.frame() # Resulting container
for( k in  seq(46,46*6,1) ){
    x1 = runif(k)
    x2 = runif(k)
    
    o <- bc(x1,x2)
    res <- rbind(res, data.frame(k = k, bc = o) )
}

yrl <- c( "yr[1]","yr[1-2]","yr[1-3]","yr[1-4]","yr[1-5]" )
g <- ggplot(res,aes(x = k, y = bc) ) + 
  theme_few() +
  geom_point() +
  geom_vline(xintercept = seq(46+46,by = 46,length.out = 4) ,linetype="dotted") +
  # Text above each vline
  annotate(geom="text",x=seq(46+23,by = 46,length.out = 5),y=max(res$bc),label=yrl,parse=T,size=6) +
#  geom_text(aes(x = seq(46+23,by = 46,length.out = 5),y = max(res$bc)),label = yrl,parse=T,inherit.aes = F) +
  # Linear regression fit
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = pretty_breaks(10),limits = c(44,280),expand = c(0,0)) + 
  scale_y_continuous(labels = NULL) + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank() ) +
  labs(x = "Number of included MODIS measurements", y = "", title = "Bray-Curtis index calculated on simulated time series with increasing length") +
  theme(plot.title = element_text(hjust = .5))
g
ggsave("PubFigures/SIFigure2.png",plot=g,dpi = 400)

broom::glance( summary( lm(bc~k,data=res) ) )
coef(lm(bc~k,data=res))
