# ------------------------------------------------ #
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
library(lme4)
library(nlme)
library(arm)
library(pbkrtest)
library(MuMIn)
library(broom) # New. Force using it
par.ori <- par(no.readonly = T)
OutputPath = "PubFigures/"
PermPath ="permutations/"
sites <- readRDS("AllPredictsSites_full.rds") # All
sites <- subset(sites, startyear > 2000)

### Construct model and visualize results ###
source("000_HelperFunction.R")
indic = "EVI2"
dism = "BCVeg" # #"SorVeg" #Dissimilarity metric
metric = "BC"# "correctedBC"
transf =FALSE # Should response be transformed?
ylist <- c(1,2,3,4,5) # Sequence of years

tg_BM <- readRDS("resSaves/MassBinStudy.rds")
# Get list of permutations
f.permut <- sort(list.files(PermPath,pattern = paste0("s_metric",dism,"_permute"),full.names = T))
res.pred <- list() # Output for all predictions
# Sequence of historic years
for(y in ylist){
  myLog("Processing past period - ",y)
  # Load matrices
  res_yearsampling <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"samplingperiod_0y-matrix.rds")) # The year of before start
  res_yearbefore <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"pastperiod_",y,"y-matrix.rds")) # The years before
  # Limit to studies with bodymass group
  res_yearsampling <- res_yearsampling[ which(names(res_yearsampling) %in% tg_BM$SS) ]
  res_yearbefore <- res_yearbefore[ which(names(res_yearbefore ) %in% tg_BM$SS) ]
  
  for(p in seq(length(f.permut))){
    myLog("Running permutation ","y",y," - ",p)
    # Load permutation
    s.metric <- readRDS(f.permut[p])
    
    # Limit to studies with bodymass group
    s.metric <- subset(s.metric, SS %in% unique(tg_BM$SS))
    
    # Available past data
    focal.period  = ymd("2000.02.18") + (365*(y+1)) # Available past in total
    # Add in startdate
    s.metric$startdate_x <- sites$Sample_start_earliest[match(s.metric$SSBS_x,sites$SSBS)]
    s.metric$startdate_y <- sites$Sample_start_earliest[match(s.metric$SSBS_y,sites$SSBS)]
    # Calculate distance in days and remove sites with greater dist of 3 months (~90days)
    s.metric <- mutate(s.metric,startd = as.numeric( abs(startdate_x - startdate_y))) %>% 
      dplyr::filter(startd <= 90) %>% # 3 months (~ 90 days)
      dplyr::filter(startdate_x > focal.period & startdate_y > focal.period) %>%  # Must have sufficient historical coverage
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
    # Get only complete studies to avoid seperate Present - Past models
    per.df <- subset(per.df,complete.cases(per.df))
    
    # Merge in bodymass information
    per.df <- merge.data.frame(per.df,tg_BM,by = "SS",all.x = T)
    per.df$Trophic_SS.abd <- factor(per.df$Trophic_SS.abd,levels = c("Autotroph","Herbivore","Carnivore","Detritivore","Omnivore"))
    per.df$MassBin.Tim.abd <- factor(per.df$MassBin.Tim.abd,levels = c(1,2,3,4))
    per.df$MassBin.Own.abd <- factor(per.df$MassBin.Own.abd,levels = c(1,2,3))
    
    # --- #
    myLog("Fitting models")
    options(na.action = na.omit)
    mlist <- list()
    ## Need to correct to dissimilarity if SorVeg
    if(dism == "SorVeg") per.df$value <- abs(1 - per.df$value)
    # Indiv coefficients
    print("All")
    if(transf){ per.df$value <- asin(sqrt(per.df$value)) } # Transform response to approach normality
    # Trophic overall models
    for(trophic in unique(per.df$Trophic_SS.abd)){
      if(is.na(trophic)) next() # Skip for non-classifiyed levels
      mod_all_pr <- glmer(value ~ EVIpres + (1+EVIpres|SS), data=per.df, subset = Trophic_SS.abd == trophic, family=gaussian)
      mod_all_pa <- glmer(value ~ EVIpast + (1+EVIpast|SS), data=per.df, subset = Trophic_SS.abd == trophic,family=gaussian)
      if(class(mod_all_pa)=="try-error") next() # Skip on failure
      # COnvergence
      if(!is.null( mod_all_pr@optinfo$conv$lme4$code)){ myLog(y," - ",p," ",trophic," - ","NOT CONVERGED");next()}
      if(!is.null( mod_all_pa@optinfo$conv$lme4$code)){ myLog(y," - ",p," ",trophic," - ","NOT CONVERGED");next()}
      print(trophic)
      mlist[["All"]][["Trophic"]][[trophic]] <- list(pr=mod_all_pr,pa=mod_all_pa)
    }
    # Bodymass groups overall models
    for(bm in unique(per.df$MassBin.Own.abd)){
      if(is.na(bm)) next() # Skip for non-classifiyed levels
      mod_all_pr <- glmer(value ~ EVIpres + (1+EVIpres|SS), subset = MassBin.Own.abd == bm, data=per.df,family=gaussian)
      mod_all_pa <- glmer(value ~ EVIpast + (1+EVIpast|SS), subset = MassBin.Own.abd == bm, data=per.df,family=gaussian)
      if(class(mod_all_pa)=="try-error") next() # Skip on failure
      # COnvergence
      if(!is.null( mod_all_pr@optinfo$conv$lme4$code)){ myLog(y," - ",p," ",bm," - ","NOT CONVERGED");next()}
      if(!is.null( mod_all_pa@optinfo$conv$lme4$code)){ myLog(y," - ",p," ",bm," - ","NOT CONVERGED");next()}
      print(bm)
      mlist[["All"]][["Bodymass"]][[bm]] <- list(pr=mod_all_pr,pa=mod_all_pa)
    }
    # Overall without plants
    for(bm in unique(per.df$MassBin.Own.abd)){
      if(is.na(bm)) next() # Skip for non-classifiyed levels
      sub <- subset(per.df,TGrouping != "Plantae" &  MassBin.Own.abd == bm)
      mod_all_pr <- glmer(value ~ EVIpres + (1+EVIpres|SS), data=sub,family=gaussian)
      mod_all_pa <- glmer(value ~ EVIpast + (1+EVIpast|SS), data=sub,family=gaussian)
      if(class(mod_all_pa)=="try-error") next() # Skip on failure
      # COnvergence
      if(!is.null( mod_all_pr@optinfo$conv$lme4$code)){ myLog(y," - ",p," ",bm," - ","NOT CONVERGED");next()}
      if(!is.null( mod_all_pa@optinfo$conv$lme4$code)){ myLog(y," - ",p," ",bm," - ","NOT CONVERGED");next()}
      print(bm)
      mlist[["AllExceptPlantae"]][["Bodymass"]][[bm]] <- list(pr=mod_all_pr,pa=mod_all_pa)
    }
    
    # For some individual groups as well
    for(taxa in unique(per.df$TGrouping)){
      sub <- subset(per.df,TGrouping == taxa)
      # Trophic
      for(trophic in unique(sub$Trophic_SS.abd)){
        if(is.na(trophic)) next() # Skip for non-classifiyed levels
        if(length(unique(sub$SS[which(sub$Trophic_SS.abd == trophic)]))<=1){
          mod_all_pr <- glm(value ~ EVIpres , data=sub, subset = Trophic_SS.abd == trophic, family=gaussian)
          mod_all_pa <- glm(value ~ EVIpast , data=sub, subset = Trophic_SS.abd == trophic, family=gaussian)
          # Skip if df of residuals is 0
          if( (df.residual(mod_all_pr) == 0) | (df.residual(mod_all_pa) == 0) ) next()
        } else {
          mod_all_pr <- glmer(value ~ EVIpres + (1+EVIpres|SS), data=sub, subset = Trophic_SS.abd == trophic, family=gaussian)
          mod_all_pa <- glmer(value ~ EVIpast + (1+EVIpast|SS), data=sub, subset = Trophic_SS.abd == trophic, family=gaussian)
          if(class(mod_all_pa)=="try-error") next() # Skip on failure
          # COnvergence
          if(!is.null( mod_all_pr@optinfo$conv$lme4$code)){ myLog(y," - ",p,trophic," - ","NOT CONVERGED");next()}
          if(!is.null( mod_all_pa@optinfo$conv$lme4$code)){ myLog(y," - ",p,trophic," - ","NOT CONVERGED");next()}
        }
        print(trophic)
        mlist[[taxa]][["Trophic"]][[trophic]] <- list(pr=mod_all_pr,pa=mod_all_pa)
      }
      # Bodymass
      for(bm in unique(sub$MassBin.Own.abd)){
        if(is.na(bm)) next() # Skip for non-classifiyed levels
        if(length(unique(sub$SS[which(sub$MassBin.Own.abd==bm)])) <= 1){ # If there is only one study
          mod_all_pr <- glm(value ~ EVIpres, subset = MassBin.Own.abd == bm, data=sub,family=gaussian)
          mod_all_pa <- glm(value ~ EVIpast, subset = MassBin.Own.abd == bm, data=sub,family=gaussian)
          # Skip if df of residuals is 0
          if( (df.residual(mod_all_pr) == 0) | (df.residual(mod_all_pa) == 0) ) next()
        } else {
          mod_all_pr <- glmer(value ~ EVIpres + (1+EVIpres|SS), subset = MassBin.Own.abd == bm, data=sub,family=gaussian)
          mod_all_pa <- glmer(value ~ EVIpast + (1+EVIpast|SS), subset = MassBin.Own.abd == bm, data=sub,family=gaussian)
          if(class(mod_all_pa)=="try-error") next() # Skip on failure
          # COnvergence
          if(!is.null( mod_all_pr@optinfo$conv$lme4$code)){ myLog(y," - ",p,bm," - ","NOT CONVERGED");next()}
          if(!is.null( mod_all_pa@optinfo$conv$lme4$code)){ myLog(y," - ",p,bm," - ","NOT CONVERGED");next()}
        }
        print(bm)
        mlist[[taxa]][["Bodymass"]][[bm]] <- list(pr=mod_all_pr,pa=mod_all_pa)
      }
      myLog(taxa," models done!")
    }
    # ------------------- #
    # # Predict for all
    # results <- data.frame()
    # for(model in names(mlist)){
    #   
    #   # For trophic level
    #   for(size in names(mlist[[model]]$Trophic) ){
    #     ## Present
    #     # Need to access something different if model is glm (studies == 1)
    #     if(any(class(mlist[[model]]$Trophic[[size]]$pr) == "glm")){
    #       pd <- data.frame(EVIpres = seq( min(mlist[[model]]$Trophic[[size]]$pr$data$EVIpres,na.rm = T),max(mlist[[model]]$Trophic[[size]]$pr$data$EVIpres,na.rm = T),length.out = 100 ), SS=NA,
    #                        Trophic_SS.abd = size)
    #       fit = predict( mlist[[model]]$Trophic[[size]]$pr, pd,type= "response")
    #     } else {
    #       pd <- data.frame(EVIpres = seq( min(mlist[[model]]$Trophic[[size]]$pr@frame$EVIpres,na.rm = T),max(mlist[[model]]$Trophic[[size]]$pr@frame$EVIpres,na.rm = T),length.out = 100 ), SS=NA,
    #                        Trophic_SS.abd = size)
    #       fit = predict( mlist[[model]]$Trophic[[size]]$pr, pd,type= "response", re.form =NA)
    #     }
    #     if(transf) fit <- sin(fit)^2 # Backtransform
    #     # Number of studies
    #     N = ifelse(any(class(mlist[[model]]$Trophic[[size]]$pr) == "glm"),1,
    #                length(unique(mlist[[model]]$Trophic[[size]]$pr@frame$SS))
    #                )
    #     N.sites = ifelse(any(class(mlist[[model]]$Trophic[[size]]$pr) == "glm"),
    #                      nrow(mlist[[model]]$Trophic[[size]]$pr$data),
    #                      nrow( mlist[[model]]$Trophic[[size]]$pr@frame)
    #                      )
    #     
    #     results <- rbind(results, data.frame(Model = model, Type = indic, Time = "Sampling conditions",focus = "Trophic",variable = size, py=y, x = pd[,1],xx = seq(1,nrow(pd)),
    #                                          N.studies = N,N.sites = N.sites,
    #                                          fit = as.numeric(fit))  )
    #     ## Past
    #     # Need to access something different if model is glm (studies == 1)
    #     if(any(class(mlist[[model]]$Trophic[[size]]$pa) == "glm")){
    #       pd <- data.frame(EVIpast = seq( min(mlist[[model]]$Trophic[[size]]$pa$data$EVIpast,na.rm = T),max(mlist[[model]]$Trophic[[size]]$pa$data$EVIpast,na.rm = T),length.out = 100 ), SS=NA,
    #                        Trophic_SS.abd = size)
    #       fit = predict( mlist[[model]]$Trophic[[size]]$pa, pd,type= "response")
    #     } else {
    #       pd <- data.frame(EVIpast = seq( min(mlist[[model]]$Trophic[[size]]$pa@frame$EVIpast,na.rm = T),max(mlist[[model]]$Trophic[[size]]$pa@frame$EVIpast,na.rm = T),length.out = 100 ), SS=NA,
    #                        Trophic_SS.abd = size)
    #       fit = predict( mlist[[model]]$Trophic[[size]]$pa, pd,type= "response", re.form =NA)
    #     }
    #     if(transf) fit <- sin(fit)^2 # Backtransform
    #     # Number of studies
    #     N = ifelse(any(class(mlist[[model]]$Trophic[[size]]$pa) == "glm"),1,
    #                length(unique(mlist[[model]]$Trophic[[size]]$pa@frame$SS))
    #     )
    #     N.sites = ifelse(any(class(mlist[[model]]$Trophic[[size]]$pa) == "glm"),
    #                      nrow(mlist[[model]]$Trophic[[size]]$pa$data),
    #                      nrow( mlist[[model]]$Trophic[[size]]$pa@frame)
    #     )
    #     
    #     results <- rbind(results, data.frame(Model = model, Type = indic, Time = "Past conditions",focus = "Trophic",variable = size, py=y, x = pd[,1],xx = seq(1,nrow(pd)),
    #                                          N.studies = N,N.sites = N.sites,
    #                                          fit = as.numeric(fit))  )
    #   }
    #   # For each BM group
    #   for(size in names(mlist[[model]]$Bodymass)){
    #     ## Present
    #     # Need to access something different if model is glm (studies == 1)
    #     if(any(class(mlist[[model]]$Bodymass[[size]]$pr) == "glm")){
    #       pd <- data.frame(EVIpres = seq( min(mlist[[model]]$Bodymass[[size]]$pr$data$EVIpres,na.rm = T),max(mlist[[model]]$Bodymass[[size]]$pr$data$EVIpres,na.rm = T),length.out = 100 ), SS=NA,
    #                        MassBin.Tim.abd = size)
    #       fit = predict( mlist[[model]]$Bodymass[[size]]$pr, pd,type= "response")
    #     } else {
    #       pd <- data.frame(EVIpres = seq( min(mlist[[model]]$Bodymass[[size]]$pr@frame$EVIpres,na.rm = T),max(mlist[[model]]$Bodymass[[size]]$pr@frame$EVIpres,na.rm = T),length.out = 100 ), SS=NA,
    #                        MassBin.Tim.abd = size)
    #       fit = predict( mlist[[model]]$Bodymass[[size]]$pr, pd,type= "response", re.form =NA)
    #     }
    #     if(transf) fit <- sin(fit)^2 # Backtransform
    #     # Number of studies
    #     N = ifelse(any(class(mlist[[model]]$Bodymass[[size]]$pr) == "glm"),1,
    #                length(unique(mlist[[model]]$Bodymass[[size]]$pr@frame$SS))
    #     )
    #     N.sites = ifelse(any(class(mlist[[model]]$Bodymass[[size]]$pr) == "glm"),
    #                      nrow(mlist[[model]]$Bodymass[[size]]$pr$data),
    #                      nrow( mlist[[model]]$Bodymass[[size]]$pr@frame)
    #     )
    #     
    #     results <- rbind(results, data.frame(Model = model, Type = indic, Time = "Sampling conditions",focus = "BM",variable = size, py=y, x = pd[,1],xx = seq(1,nrow(pd)),
    #                                          N.studies = N,N.sites = N.sites,
    #                                          fit = as.numeric(fit))  )
    #     ## Past
    #     # Need to access something different if model is glm (studies == 1)
    #     if(any(class(mlist[[model]]$Bodymass[[size]]$pa) == "glm")){
    #       pd <- data.frame(EVIpast = seq( min(mlist[[model]]$Bodymass[[size]]$pa$data$EVIpast,na.rm = T),max(mlist[[model]]$Bodymass[[size]]$pa$data$EVIpast,na.rm = T),length.out = 100 ), SS=NA,
    #                        MassBin.Tim.abd = size)
    #       fit = predict( mlist[[model]]$Bodymass[[size]]$pa, pd,type= "response")
    #     } else {
    #       pd <- data.frame(EVIpast = seq( min(mlist[[model]]$Bodymass[[size]]$pa@frame$EVIpast,na.rm = T),max(mlist[[model]]$Bodymass[[size]]$pa@frame$EVIpast,na.rm = T),length.out = 100 ), SS=NA,
    #                        MassBin.Tim.abd = size)
    #       fit = predict( mlist[[model]]$Bodymass[[size]]$pa, pd,type= "response", re.form =NA)
    #     }
    #     if(transf) fit <- sin(fit)^2 # Backtransform
    #     # Number of studies
    #     N = ifelse(any(class(mlist[[model]]$Bodymass[[size]]$pa) == "glm"),1,
    #                length(unique(mlist[[model]]$Bodymass[[size]]$pa@frame$SS))
    #     )
    #     N.sites = ifelse(any(class(mlist[[model]]$Bodymass[[size]]$pa) == "glm"),
    #                      nrow(mlist[[model]]$Bodymass[[size]]$pa$data),
    #                      nrow( mlist[[model]]$Bodymass[[size]]$pa@frame)
    #     )
    #     
    #     results <- rbind(results, data.frame(Model = model, Type = indic, Time = "Past conditions",focus = "BM",variable = size, py=y, x = pd[,1],xx = seq(1,nrow(pd)),
    #                                          N.studies = N,N.sites = N.sites,
    #                                          fit = as.numeric(fit))  ) 
    #   }
    #   
    #   }
    # 
    # # Save in list the prediction
    # res.pred[["PermRuns"]] <- rbind(res.pred[["PermRuns"]],
    #                                 results %>% mutate(PermRun = p) )
    
    # # Save coefficients
    co <- data.frame()
    for(model in names(mlist)){
      ## Trophic
      for( size in names(mlist[[model]]$Trophic)){
        mod <- mlist[[model]]$Trophic[[size]]$pr
        if(any(class(mod) == "glm")){
          # Special if a glm was fitted (single study)
          N = 1
          N.sites = nrow(mod$data)
          # Present
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(coef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-p.value),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Sampling conditions",focus = "Trophic",variable = size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites)
            )
          }
          # Past
          mod <- mlist[[model]]$Trophic[[size]]$pa
          N = 1
          N.sites = nrow(mod$data)
          
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(coef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-p.value),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Past conditions",focus = "Trophic",variable=size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites )
            )
          }
        } else {
          N = length(unique(mod@frame$SS))
          N.sites = nrow(mod@frame)
          # Present
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(fixef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-group),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Sampling conditions",focus = "Trophic",variable = size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites)
            )
          }
          # Past
          mod <- mlist[[model]]$Trophic[[size]]$pa
          N = length(unique(mod@frame$SS))
          N.sites = nrow(mod@frame)
          
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(fixef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-group),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Past conditions",focus = "Trophic",variable=size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites )
            )
          }
        }
        
        
      }
      ## Bodymass
      for( size in names(mlist[[model]]$Bodymass)){
        mod <- mlist[[model]]$Bodymass[[size]]$pr
        if(any(class(mod) == "glm")){
          # Special if a glm was fitted (single study)
          N = 1
          N.sites = nrow(mod$data)
          # Present
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(coef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-p.value),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Sampling conditions",focus = "BM",variable = size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites)
            )
          }
          # Past
          mod <- mlist[[model]]$Bodymass[[size]]$pa
          N = 1
          N.sites = nrow(mod$data)
          
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(coef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-p.value),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Past conditions",focus = "BM",variable=size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites )
            )
          }
        } else {
          N = length(unique(mod@frame$SS))
          N.sites = nrow(mod@frame)
          # Present
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(fixef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-group),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Sampling conditions",focus = "BM",variable = size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites)
            )
          }
          # Past
          mod <- mlist[[model]]$Bodymass[[size]]$pa
          N = length(unique(mod@frame$SS))
          N.sites = nrow(mod@frame)
          
          if(!is.null(mod)){
            ci <- broom::confint_tidy(mod,method = "Wald") %>% subset(.,complete.cases(.)) %>%  mutate(term = names(fixef(mod)))
            eff = left_join(ci, broom::tidy(mod) %>% dplyr::select(-group),by="term" ) %>% 
              dplyr::filter(term != "(Intercept)")# Remove intercept
            
            co <- rbind(co,
                        data.frame(Model = model,PermRun = p, py = y,Time="Past conditions",focus = "BM",variable=size,
                                   R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),
                                   eff,
                                   N = N,N.sites = N.sites )
            )
          }
        }
        
      }
      
    }
    
    
    
    # ---------------------------------- #
    
    
    res.pred[["Coeff"]] <- rbind( res.pred[["Coeff"]], co )
    rm(results,mlist,pd,fit,mod_all_pr,mod_all_pa,mod)
  }
}
# --- #
saveRDS(res.pred,"resSaves/outBCpermute_minStartEVI2_BC_BodyMassINT_OWN.rds")
# DONE
# --- #
