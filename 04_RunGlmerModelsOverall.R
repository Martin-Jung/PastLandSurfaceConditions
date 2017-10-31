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
library(gamm4)
library(broom) # New. Force using it
par.ori <- par(no.readonly = T)
OutputPath = "PubFigures/"
PermPath ="permutations/"
sites <- readRDS("AllPredictsSites_full.rds") # All
sites <- subset(sites, startyear > 2000)

### Construct model and visualize results ###
source("000_HelperFunction.R")
source("0001_GLMERFunction.R")
indic = "EVI2"#
dism = "BCVeg"# #Dissimilarity metric
metric =  "BC" #"_meanSPD"# "correctedBC"
transf = F # Should response be transformed?
ylist <- c(1,2,3,4,5) # Sequence of years
# Get list of permutations
f.permut <- sort(list.files(PermPath,pattern = paste0("s_metric",dism,"_permute"),full.names = T))
res.pred <- list() # Output for all predictions

stopifnot(exists("extractCoef")) # load the function first

# Sequence of historic years
for(y in ylist){
  myLog("Processing past period - ",y)
  # Load matrices
  #res_yearsampling <- readRDS(paste0("res_yearbefored",metric,"_minStart_",indic,"samplingperiod_INDIVIDUAL_0y-matrix.rds")) # The year of before start
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
      
      m <- reshape2::melt( res_yearsampling[[as.character(study)]] ) %>% dplyr::rename(SSBS_x = Var1,SSBS_y = Var2,EVIpres = value)
      m2 <- reshape2::melt( res_yearbefore[[as.character(study)]] ) %>% dplyr::rename(SSBS_x = Var1,SSBS_y = Var2,EVIpast = value)
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
    if(transf) { per.df$value <- asin(sqrt(per.df$value)) }
    
    # Normal models
    mod_all_pr <- glmer(value ~ EVIpres + (1+EVIpres|SS), data=per.df,family=gaussian)
    mod_all_pa <- glmer(value ~ EVIpast + (1+EVIpast|SS), data=per.df,family=gaussian)
    
    if(exists("res.dist")){
      # Added covariate with logDistance # Interaction was insignificant
      mod_all_pr.ld <- glmer(value ~ EVIpres+log10Distance + (1+EVIpres|SS), data=per.df,family=gaussian)
      mod_all_pa.ld <- glmer(value ~ EVIpast+log10Distance + (1+EVIpast|SS), data=per.df,family=gaussian)
      
      mlist[["All"]] <- list(evi=list(pr = mod_all_pr, pa=mod_all_pa) ,
                             ld = list(pr = mod_all_pr.ld, pa=mod_all_pa.ld))
    } else {
      mlist[["All"]] <- list(evi=list(pr = mod_all_pr, pa=mod_all_pa))
      
    }
    
    for(taxa in unique(per.df$TGrouping)){
      test = subset(per.df,TGrouping==taxa)
      mod_pr <- try( glmer(value ~ EVIpres + (1+EVIpres|SS), data=na.omit(test), family=gaussian),silent=T )
      if(class(mod_pr)=="try-error") next() # Skip on failure
      if(!is.null( mod_pr@optinfo$conv$lme4$code)){ myLog(y," - ",p,"-",taxa," - ","NOT CONVERGED");next()}
      test = subset(per.df,TGrouping==taxa) %>% dplyr::select(-EVIpres)
      mod_pa <- try( glmer(value ~ EVIpast + (1+EVIpast|SS), data=na.omit(test), family=gaussian),silent=T )
      if(!is.null( mod_pa@optinfo$conv$lme4$code)){ myLog(y," - ",p,"-",taxa," - ","NOT CONVERGED");next()}
      if(exists("res.dist")){
        # Interaction with logDistance
        test = subset(per.df,TGrouping==taxa) %>% dplyr::select(-EVIpast)
        mod_pr.ld <- glmer(value ~ EVIpres+log10Distance + (1+EVIpres|SS), data=na.omit(test),family=gaussian)
        test = subset(per.df,TGrouping==taxa) %>% dplyr::select(-EVIpres)
        mod_pa.ld <- glmer(value ~ EVIpast+log10Distance + (1+EVIpast|SS), data=na.omit(test),family=gaussian)
      }      
      print(taxa)
      if(exists("res.dist")){
        mlist[[taxa]] <- list(evi=list(pr = mod_pr, pa = mod_pa),ld = list(pr = mod_pr.ld, pa=mod_pa.ld))
      } else {
        mlist[[taxa]] <- list(evi=list(pr = mod_pr, pa = mod_pa))
      }
    }
    
    # Predict for all
    results <- data.frame()
    for(model in names(mlist)){
      
      pd <- data.frame(EVIpres = seq( min(mlist[[model]]$evi$pr@frame$EVIpres),max(mlist[[model]]$evi$pr@frame$EVIpres),length.out = 1000 ), SS=NA)
      fit = predict( mlist[[model]]$evi$pr, pd,se.fit =T, re.form =NA)
      if(transf) { fit$fit <- sin(fit$fit)^2  }# Backtransform
      
      results <- rbind(results, data.frame(Model = model, Type = indic, Time = "Present",py=y, x = pd[,1],xx = seq(1,nrow(pd)),
                                           N.studies = length(unique(mlist[[model]]$evi$pr@frame$SS)),
                                           fit = fit$fit, se = fit$se.fit)  )
      
      pd <- data.frame(EVIpast = seq( min(mlist[[model]]$evi$pa@frame$EVIpast),max(mlist[[model]]$evi$pa@frame$EVIpast),length.out = 1000 ), SS=NA)
      fit = predict( mlist[[model]]$evi$pa, pd,se.fit =T, re.form =NA)
      if(transf) { fit$fit <- sin(fit$fit)^2 } # Backtransform
      
      results <- rbind(results, data.frame(Model = model, Type = indic, Time = "Past",py=y, x = pd[,1],xx = seq(1,nrow(pd)),
                                           N.studies = length(unique(mlist[[model]]$evi$pa@frame$SS)),
                                           fit = fit$fit, se = fit$se.fit)  )
      
      # pd <- data.frame(EVIpres = seq( min(mlist[[model]]$evi$prpa@frame$EVIpres),max(mlist[[model]]$evi$prpa@frame$EVIpres),length.out = 1000 ),
      #                  EVIpast = seq( min(mlist[[model]]$evi$prpa@frame$EVIpast),max(mlist[[model]]$evi$prpa@frame$EVIpast),length.out = 1000 ), SS=NA)
      # fit = predict( mlist[[model]]$evi$prpa, pd,se.fit =T, re.form =NA)
      # if(transf) fit$fit <- sin(fit$fit)^2 # Backtransform
      # 
      # results <- rbind(results, data.frame(Model = model, Type = indic, Time = "PresentPast",py=y, x = pd[,1],xx = seq(1,nrow(pd)),
      #                                      N.studies = length(unique(mlist[[model]]$evi$prpa@frame$SS)),
      #                                      fit = fit$fit, se = fit$se.fit)  )
    }
    
    # Save in list the prediction
    res.pred[["PermRuns"]] <- rbind(res.pred[["PermRuns"]],
                                    results %>% mutate(PermRun = p) )
    
    # Save coefficients
    co <- data.frame()
    for(model in names(mlist)){
      # Present
      mod <- mlist[[model]]$evi$pr
      ci <- confint.merMod(mod,method = "Wald")
      # Use a kenwards roger approximation for testing the p-value
      kr <- KRmodcomp.mer(mod,glmer(value~1 + (1|SS),data=mod@frame))
      
      co <- rbind(co,
                  data.frame(Model = model,PermRun = p, py = y,Time="Present",eff = as.numeric(fixef(mod)[2]),se = as.numeric(se.fixef(mod))[2],
                             R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),pval = getKR(kr,"p.value"),
                             ci.low = ci[4,1],ci.high = ci[4,2],N = length(unique(mod@frame$SS)),N.sites = nrow(mod@frame) )
      )
      # Past
      mod <- mlist[[model]]$evi$pa
      ci <- confint.merMod(mod,method = "Wald")
      # Use a kenwards roger approximation for testing the p-value
      kr <- KRmodcomp.mer(mod,glmer(value~1 + (1|SS),data=mod@frame))
      co <- rbind(co,
                  data.frame(Model = model,PermRun = p, py = y,Time="Past",eff = as.numeric(fixef(mod)[2]),se = as.numeric(se.fixef(mod))[2],
                             R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),pval = getKR(kr,"p.value"),
                             ci.low = ci[4,1],ci.high = ci[4,2],N = length(unique(mod@frame$SS)),N.sites = nrow(mod@frame) )
      )
      
      # PresentPast
      # mod <- mlist[[model]]$evi$prpa
      # ci <- confint.merMod(mod,method = "Wald")
      # # Use a kenwards roger approximation for testing the p-value
      # kr <- KRmodcomp.mer(mod,glmer(value~1 + (1|SS),data=mod@frame))
      # co <- rbind(co,
      #             data.frame(Model = model,PermRun = p, py = y,Time="PresentPast",eff = as.numeric(fixef(mod)[2]),se = as.numeric(se.fixef(mod))[2],
      #                        R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),pval = getKR(kr,"p.value"),
      #                        ci.low = ci[4,1],ci.high = ci[4,2],N = length(unique(mod@frame$SS)),N.sites = nrow(mod@frame) )
      # )
    }
    res.pred[["Coeff"]] <- rbind( res.pred[["Coeff"]], co )
    
    # Save coefficients for interactive model
    if(exists("res.dist")){
      co <- data.frame()
      for(model in names(mlist)){
        # Present
        mod <- mlist[[model]]$ld$pr
        ci <- confint.merMod(mod,method = "Wald")
        # Use a kenwards roger approximation for testing the p-value
        kr <- KRmodcomp.mer(mod,glmer(value~EVIpres + (1+EVIpres|SS),data=mod@frame))
        
        co <- rbind(co,
                    data.frame(Model = model,PermRun = p, py = y,Time="Present",eff = as.numeric(fixef(mod)[2]),se = as.numeric(se.fixef(mod)[2]),
                               R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),pval = getKR(kr,"p.value"),
                               ci.low = ci[6,1],ci.high = ci[6,2],N = length(unique(mod@frame$SS)),N.sites = nrow(mod@frame) )
        )
        # Past
        mod <- mlist[[model]]$ld$pa
        ci <- confint.merMod(mod,method = "Wald")
        # Use a kenwards roger approximation for testing the p-value
        kr <- KRmodcomp.mer(mod,glmer(value~EVIpast + (1+EVIpast|SS),data=mod@frame))
        co <- rbind(co,
                    data.frame(Model = model,PermRun = p, py = y,Time="Past",eff = as.numeric(fixef(mod)[2]),se = as.numeric(se.fixef(mod)[2]),
                               R2m = as.numeric(r.squaredGLMM(mod)[1]), r2c = as.numeric(r.squaredGLMM(mod)[2]),pval = getKR(kr,"p.value"),
                               ci.low = ci[6,1],ci.high = ci[6,2],N = length(unique(mod@frame$SS)),N.sites = nrow(mod@frame) )
        )
        
      }
      res.pred[["LDCoeff"]] <- rbind( res.pred[["LDCoeff"]], co )
    }
    
    # Save coefficients per studies for each model
    res_co <- data.frame()
    for(model in names(mlist)){
      mod_pr <- mlist[[model]]$evi$pr
      mod_pa <- mlist[[model]]$evi$pa
      prc = extractCoef(mod_pr,group = "SS")
      pac = extractCoef(mod_pa,group = "SS")
      pr = data.frame(y=y,py=p,TGrouping = model,time = "Sampling conditions",SS = row.names(prc), value = prc$slope, value.se = prc$slope_se)
      pa = data.frame(y=y,py=p,TGrouping = model,time = "Past conditions",SS = row.names(pac), value = pac$slope, value.se = pac$slope_se)
      res_co <- rbind(res_co, rbind(pr,pa) )
      rm(mod_pr,mod_pa,prc,pac,pr,pa)
    }
    res.pred[["StudyCoeff"]] <- rbind(res.pred[["StudyCoeff"]], res_co )
    
    rm(results,mlist,pd,fit,mod_all_pr,mod_all_pa,mod_all_prpa)
  }
}
# Individual with new one
saveRDS(res.pred,paste0("resSaves/outBCpermute_minStart5y",indic,"_",metric,"_",dism,".rds")) # dBC both axis
#saveRDS(res.pred,paste0("resSaves/outBCpermute_minStart5yINDIVIDUAL",indic,"_",metric,"_",dism,".rds")) # dBC both axis
#saveRDS(res.pred,paste0("resSaves/outBCpermute_",indic,"_",dism,".rds")) # dBC both axis
myLog("Done!!!")
rm(res.pred)
