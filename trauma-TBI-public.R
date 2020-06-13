############################## Description ##############################

## This is a commented file containing code to perform secondary analyses of PAMPer trial data and TBI biomarker data (Sperry et al., 2018, NEJM)
## Project Support: NIH T32 in Trauma and Sepsis, University of Pittsburgh Department of Surgery
## PI: Dr. Jason Sperry
## Written by Danielle Gruen, October 2019
## Last updated: June 2020

############################## Workspace Setup ##############################
## Set Up Workspace
    rm(list=ls()) # Clear all variables; delete all the objects in the workspace and start with a clean slate
    setwd("/Users/dgruen/Documents/Academia/Pitt_T32/projects/pamper/data/") # set the working directory
    set.seed(386) # make code reproducible

## Dependencies (these libraries must be loaded for code and analyses)
    # library(gmodels) # cross tabs for chi square
    library(yardstick) # ROC curves & AUC
    library(haven)  # reads in SPSS files
    library(gdata)  # reads excel files and more
    library(RColorBrewer) # colors for graphs
    library(scales) # pretty plots
    library(reshape2) # used for plotting with ggplot2
    library(survminer) # creates survival curves
    library(survival) # required for survminer & cox hazard
    library(caTools) # classification & splitting data for validation etc. 
    library(broom) # way of dealing with model output results in a tidy way, need for GLM
    library(dendextend) # dendrogram building and visualization
    library(RCurl) # cross validation
    library(prettyR) # cross validation
    library(lsmeans) # Least-squares means
    library(multcompView) # Summarize multiple paired comparisons
    library(rcompanion) # Variety of Functions to Support Extension Education Program Evaluation
    library(nlme) # Anova
    library(devtools) # developer tools
    library(tableone) # to make a summmary table
    library(np) # nonparametric linear regression
    library(xtable) # make things go to latex e.g. tables
    library(jtools) # good for stats summaries
    library(ggstance) # required for jtools / summ
    library(sandwich) # robust standard errors for linear regression
    library(lmtest) # for lin reg
    library(multiwayvcov)  # for lin reg
    library(wesanderson) # color palette for ggplot
    library(tidyverse)  # for nice visuals (also includes things like haven and ggplot2)
    library(gmodels)
    library(boot) # for bootstrapping
    library(ranger)
    library(ggfortify)
    library(rms)
    library(pec)
    library(corrplot)
    library(mlbench)
    library(gee)
    library(geepack)
    library(cmprsk)
    library(frailtySurv)
    library(fragilityindex)

############################## File Setup ##############################
    
############################## Clean Data ##############################
    
############################## Data Inspection ##############################
    
############################## Tables and Figures ##############################   
    
## CONSORT diagram
    df <- df_master_filtered_wide
    tally_df <- df %>%
      group_by(FFP, TBI) %>%
      tally()
    tally_df
    
## Table 1: Characteristics of TBI vs. No TBI
    df <- df_master_filtered_wide
    df <- filter(df, !is.na(FFP))
    df <- filter(df, !is.na(TBI))
    
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    df$gcs_cat <- NA
    df$gcs_cat[which(df$total_gcs_score < 8)] = 0
    df$gcs_cat[which(df$total_gcs_score >= 8 & df$total_gcs_score<13)] = 1
    df$gcs_cat[which(df$total_gcs_score >= 13)] = 2
    df$gcs_cat <- as.factor(df$gcs_cat)
    
    df$ED_gcs_cat <- NA
    df$ED_gcs_cat[which(df$ed_initial_gcs < 8)] = 0
    df$ED_gcs_cat[which(df$ed_initial_gcs >= 8 & df$ed_initial_gcs<13)] = 1
    df$ED_gcs_cat[which(df$ed_initial_gcs >= 13)] = 2
    df$ED_gcs_cat <- as.factor(df$ED_gcs_cat)
    
    
    myVars <- c("age", "gender",  "race", "vitk", "antiplt",
                
                "iss","ais_head", "ais_chest", "ais_abdomen", "ais_extremity",
                "gcs_cat", "ED_gcs_cat", "PH_sbp_70", "TBI", "spinal_cord_injury", "any_blunt",
                "blunt_mechanism",
                
                "PH_intubation", "PH_CPR", "plasma_on_helicopter", "PH_crystalloid",
                "PH_prbc", "PH_blood", "PH_time", "transfer",
                
                "INR", "ACT", "alpha", "K", "MA", "LY30", 
                "transfusion_24h", "prbc_24h", "plasma_24h",
                "platelets_24h", "crystalloid_24h", "vaso_24h", "prbc_10_24h",
                "crani_24",
                
                "t_30d_censor_h", "mortality_24h", "MOF", "icu_los", "hospital_los", "mech_vent_days")
    
    catVars <- c("gender","race", "vitk", "antiplt", "t_30d_censor_h", "TBI", "severe_head", 
                 "spinal_cord_injury", "MOF", "PH_intubation", "PH_sbp_70", "any_blunt", "blunt_mechanism",
                 "PH_CPR", "plasma_on_helicopter", "transfer", "PH_blood", "mortality_24h", "vaso_24h", "prbc_10_24h",
                 "crani_24", "PH_time_high", "gcs_cat", "ED_gcs_cat")
    
    # with both p values and standardized mean differences
    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "TBI", factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    View(tab1)
    
## Table 2: Within TBI, Characteristics of prehospital plasma vs. standard care
    df <- df_master_filtered_wide
    df <- filter(df, !is.na(FFP))
    df <- filter(df, !is.na(TBI))
    df <- filter(df, TBI==1)
    
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    df$gcs_cat <- NA
    df$gcs_cat[which(df$total_gcs_score < 8)] = 0
    df$gcs_cat[which(df$total_gcs_score >= 8 & df$total_gcs_score<13)] = 1
    df$gcs_cat[which(df$total_gcs_score >= 13)] = 2
    df$gcs_cat <- as.factor(df$gcs_cat)
    
    df$ED_gcs_cat <- NA
    df$ED_gcs_cat[which(df$ed_initial_gcs < 8)] = 0
    df$ED_gcs_cat[which(df$ed_initial_gcs >= 8 & df$ed_initial_gcs<13)] = 1
    df$ED_gcs_cat[which(df$ed_initial_gcs >= 13)] = 2
    df$ED_gcs_cat <- as.factor(df$ED_gcs_cat)
    
    myVars <- c("age", "gender",  "race", "vitk", "antiplt",
                
                "iss","ais_head", "ais_chest", "ais_abdomen", "ais_extremity",
                "gcs_cat", "ED_gcs_cat", "PH_sbp_70", "TBI", "spinal_cord_injury", "any_blunt",
                "blunt_mechanism",
                
                "PH_intubation", "PH_CPR", "plasma_on_helicopter", "PH_crystalloid",
                "PH_prbc", "PH_blood", "PH_time", "transfer",
                
                "INR", "ACT", "alpha", "K", "MA", "LY30", 
                "transfusion_24h", "prbc_24h", "plasma_24h",
                "platelets_24h", "crystalloid_24h", "vaso_24h", "prbc_10_24h",
                "crani_24",
                
                "t_30d_censor_h", "mortality_24h", "MOF", "icu_los", "hospital_los", "mech_vent_days")
    
    catVars <- c("gender","race", "vitk", "antiplt", "t_30d_censor_h", "TBI", "severe_head", 
                 "spinal_cord_injury", "MOF", "PH_intubation", "PH_sbp_70", "any_blunt", "blunt_mechanism",
                 "PH_CPR", "plasma_on_helicopter", "transfer", "PH_blood", "mortality_24h", "vaso_24h", "prbc_10_24h",
                 "crani_24", "PH_time_high", "gcs_cat", "ED_gcs_cat")
    
    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "FFP", factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    View(tab1)
    
    
#### TBI Subtype #### 
    df <- Subj_Initial_Injury 
    df <- select(df, "stnum", "traumatic_brain_injury_types")
    df <- filter(df, !is.na(traumatic_brain_injury_types))
    
    df_tally <- df %>%
      group_by(traumatic_brain_injury_types) %>%
      tally()
    df_tally    
    
    # TBI Subtypes (SUBJ_Initial_Injury	traumatic_brain_injury_types) (1|2|3|4	Subdural hematoma/hemorrhage|Epidural hematoma/hemorrhage|Contusion|Other)
    
    df_sub <- df %>% separate(traumatic_brain_injury_types, sep='\\|', into = c('a','b','c','d'), remove=T) #%>%
    #filter(a != '') #%>% spread(a)
    for (i in 2:5) df_sub[,i] <- as.integer(df_sub[,i])
    df_sub$injury1 <- 0; df_sub$injury2 <- 0; df_sub$injury3 <- 0; df_sub$injury4 <- 0
    for (i in 1:nrow(df_sub)){
      if (length(which(df_sub[i,c(2:5)] == 1)) > 0) df_sub$injury1[i] <- 1
      if (length(which(df_sub[i,c(2:5)] == 2)) > 0) df_sub$injury2[i] <- 1
      if (length(which(df_sub[i,c(2:5)] == 3)) > 0) df_sub$injury3[i] <- 1
      if (length(which(df_sub[i,c(2:5)] == 4)) > 0) df_sub$injury4[i] <- 1
      
    }
    
    df_new <- df_sub %>% select(-c('a','b','c','d'))
    
    
    ## Tally each type
    df <- df_master_filtered_wide
    df <- select(df, "stnum", "TBI", "FFP", "alive_at_30")
    
    df_merged <- merge(df, df_new, by="stnum")
    
    df_tally <- df_merged %>%
      group_by(FFP) %>%
      tally()
    df_tally    
    
    
    df <- df_merged
    df <-filter(df, TBI==1)
    
    myVars <- c("injury1", "injury2", "injury3", "injury4")
    catVars <- c("injury1", "injury2", "injury3", "injury4")
    
    tab1 <- CreateTableOne(vars = myVars, data = df, factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    View(tab1)
    
    
#### Cause of Death #### 
    df <- df_master_filtered_wide
    df <- filter(df, TBI==1)
    df <- filter(df, alive_at_30==2)
    
    myVars <- c("cause_of_death")
    catVars <- c("cause_of_death")
    
    tab1 <- CreateTableOne(vars = myVars, data = df, factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    View(tab1)
    
## Analysis: Generalized Estimating Equations (GEE) (account for site clusters)
    ## Data frame
    df <- df_master_filtered_wide
    df$ais_head <- as.numeric(as.character(df$ais_head))
    
    df <- select(df, stnum, t_30d_censor_h, FFP, TBI, alive_at_30,
                 iss, MOI, SiteID, PH_time_high,
                 PH_prbc, PH_crystalloid, PH_blood,
                 PH_sbp_70, age, gender, initial_GCS_8, ais_head)
    
    df <- na.omit(df)

    fit.gee <- geeglm(t_30d_censor_h ~ TBI * FFP + iss + ais_head +
                        PH_prbc + PH_crystalloid + PH_blood + 
                        MOI + PH_sbp_70  + PH_time_high +
                        age + gender + initial_GCS_8
                        ,
                      family = binomial,
                      data = df, 
                      id = SiteID,
                      corstr = "exchangeable")
    summary(fit.gee)
    
    # Note also test iss*FFP interaction term
    

## Figure: Kaplan Meier Survival
    ## No TBI
    df <- df_master_filtered_wide
    df <- filter(df, TBI==2)
    
    fit <- survfit(Surv(t_30d_mort_h, t_30d_censor_h) ~ plasma_on_helicopter, data = df) # Fit survival curves 
    
    ggsurvplot(fit,
               legend = "bottom", 
               legend.title = "Prehospital Plasma",
               break.time.by = 250, # break time axis by 168 hours (make it weeks)
               legend.labs = c("Plasma", "No Plasma"),
               # conf.int = TRUE, # Add confidence interval
               pval = TRUE,
               risk.table = TRUE,
               cumcensor = TRUE) # Add p-value
    
    ## TBI
    df <- df_master_filtered_wide
    df <- filter(df, TBI==1)
    
    fit <- survfit(Surv(t_30d_mort_h, t_30d_censor_h) ~ plasma_on_helicopter, data = df) # Fit survival curves 
    
    ggsurvplot(fit,
               legend = "bottom", 
               legend.title = "Prehospital Plasma",
               break.time.by = 250, # break time axis by 168 hours (make it weeks)
               legend.labs = c("Plasma", "No Plasma"),
               # conf.int = TRUE, # Add confidence interval
               pval = TRUE,
               risk.table = TRUE,
               cumcensor = TRUE) # Add p-value
    

## Analysis: Cox Hazard Model TBI and No TBI

    df <- df_master_filtered_wide
    df <- filter(df, TBI==1)

    df <- select(df, "stnum", "TBI", "t_30d_mort_h", "t_30d_censor_h",
                 "age", "total_gcs_score", "PH_crystalloid", "FFP", "PH_prbc", "PH_intubation", "PH_time_high",
                 "gender", "transfer", "SiteID", "MOI", "iss", "PH_sbp_70", "initial_GCS_8",
                 "ais_head", "ais_face", "ais_chest", "ais_abdomen", "ais_extremity", "ais_external")
    

    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    df$SiteID <- as.factor(df$SiteID)
    
    # category for other body regions >3
    df$AIS_cat <- NA
    df$AIS_cat[which(df$ais_chest<3 & df$ais_abdomen<3 & df$ais_extremity<3)] = 0 
    df$AIS_cat[which(df$ais_chest>=3 | df$ais_abdomen>=3 | df$ais_extremity>=3)] = 1 
    df$AIS_cat <- as.factor(df$AIS_cat)
    
    tally_df <- df %>%
      group_by(AIS_cat) %>%
      tally()
    tally_df
    
    res.frailty <- coxph(Surv(t_30d_mort_h, t_30d_censor_h) ~ 
                       age + 
                       total_gcs_score +
                       AIS_cat +
                       ais_head +
                       PH_crystalloid + 
                       FFP + 
                       PH_prbc +
                       PH_intubation +
                       PH_time_high +
                       gender +
                       PH_sbp_70 +
                       transfer +
                       frailty(SiteID, distribution="gamma"),
                     data =  df)
    
    summary(res.frailty)
    
    vif(res.frailty)
    
    test.ph <- cox.zph(res.frailty)
    test.ph
   
    
## Figure: box plot SC vs FFP
    df <- df_master_filtered_long
    df <- filter(df, !is.na(hour)) # get rid of NA
    df <- filter(df, Biomarker=="GFAP" | Biomarker=="UCH_L1")
    ggplot(df %>% 
             unite(twoby, hour, FFP, remove = F), 
           aes(x = twoby,
               y = Concentration,
               fill=as.factor(FFP))) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width=0.3, alpha=0.2) +
      facet_wrap(~Biomarker, scales = "free") +
      theme_classic() +
      theme_classic() +
      scale_fill_manual(values = c("#ffffff", "#969696")) +
      theme(aspect.ratio=1)
    
## Figure: box plot TBI vs. No TBI
    df <- df_master_filtered_long
    df <- filter(df, !is.na(hour)) # get rid of NA
    df <- filter(df, Biomarker=="GFAP" | Biomarker=="UCH_L1")
    ggplot(df %>% 
             unite(twoby, hour, TBI, remove = F), 
           aes(x = twoby,
               y = Concentration,
               fill=as.factor(TBI))) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width=0.3, alpha=0.2) +
      facet_wrap(~Biomarker, scales = "free") +
      theme_classic() +
      theme_classic() +
      scale_fill_manual(values = c("#ffffff", "#969696")) +
      theme(aspect.ratio=1)

  
## Figure: Kaplan Meier early vs late plasma; Scene PP vs SC vs Transfer PP vs SC
    ## Transfer
    df <- df_master_filtered_wide
    df <- filter(df, transfer==1)
    df <- filter(df, TBI==1)
    
    fit <- survfit(Surv(t_30d_mort_h, t_30d_censor_h) ~ plasma_on_helicopter, data = df) # Fit survival curves 
    
    ggsurvplot(fit,
               legend = "bottom", 
               legend.title = "Prehospital Plasma",
               break.time.by = 250, # break time axis by 168 hours (make it weeks)
               legend.labs = c("Plasma", "No Plasma"),
               # conf.int = TRUE, # Add confidence interval
               pval = TRUE,
               risk.table = TRUE,
               cumcensor = TRUE) # Add p-value
    
    ## Scene
    df <- df_master_filtered_wide
    df <- filter(df, transfer==0)
    df <- filter(df, TBI==1)
    
    fit <- survfit(Surv(t_30d_mort_h, t_30d_censor_h) ~ plasma_on_helicopter, data = df) # Fit survival curves 
    
    ggsurvplot(fit,
               legend = "bottom", 
               legend.title = "Prehospital Plasma",
               break.time.by = 250, # break time axis by 168 hours (make it weeks)
               legend.labs = c("Plasma", "No Plasma"),
               # conf.int = TRUE, # Add confidence interval
               pval = TRUE,
               risk.table = TRUE,
               cumcensor = TRUE) # Add p-value
  
## Forest plot
    ## GCS Groups
    df <- df_master_filtered_wide
    df$head_groups_dsg <- NA
    df$head_groups_dsg[which(df$TBI == 2)] = 0
    df$head_groups_dsg[which(df$TBI == 1 & df$total_gcs_score > 8)] = 1 
    df$head_groups_dsg[which(df$TBI == 1 & df$total_gcs_score <= 8)] = 2 
    df$head_groups_dsg <- as.factor(df$head_groups_dsg) # make it a factor
    df_master_filtered_wide_head_groups <- df
    
    ## Polytrauma Groups
    df <- df_master_filtered_wide
    df$trauma_groups_dsg <- NA
    df$trauma_groups_dsg[which(df$TBI == 2)] = 0
    df$trauma_groups_dsg[which(df$TBI == 1 & df$ais_abdomen<3 & df$ais_chest<3 & df$ais_extremity<3)] = 1
    df$trauma_groups_dsg[which(df$TBI == 1 & (df$ais_abdomen>=3 | df$ais_chest>=3 | df$ais_extremity>=3))] = 2
    df$trauma_groups_dsg <- as.factor(df$trauma_groups_dsg) # make it a factor
    df_master_filtered_wide_head <- df
    
    ## Generate HR
    df <- df_master_filtered_wide_head
    df <- filter(df, !is.na(trauma_groups_dsg))
    df <- filter(df, !is.na(alive_at_30))
    
    df <- df_master_filtered_wide_head
    df <- filter(df, trauma_groups_dsg==2)
    
    res.cox <- coxph(Surv(t_30d_mort_h, t_30d_censor_h) ~ FFP, data =  df)
    
    test <- summary(res.cox)
    coef(summary(res.cox))
    
    ggforest(res.cox, data = df)
   
    
    ## Forest Plot
    RR_data <- read.xls("cox_output_all_groups.xlsx")
    p = ggplot(data=RR_data,
               aes(x = Group,y = RiskRatio, ymin = LowerLimit, ymax = UpperLimit ))+
      geom_pointrange(aes(col=Group))+
      geom_hline(aes(fill=Group),yintercept =1, linetype=1)+
      geom_hline(aes(fill=Group),yintercept =0.64, linetype=2)+
      xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
      geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.5,cex=1)+ 
      facet_wrap(~Condition,strip.position="left",nrow=9,scales = "free_y") +
      theme(plot.title=element_text(size=16,face="bold"),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_text(face="bold"),
            axis.title=element_text(size=12,face="bold"),
            strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
      coord_flip() +
      scale_y_log10(limits=c(0.01,12))
    p

