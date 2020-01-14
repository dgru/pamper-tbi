############################## Description ##############################

## This is a commented file containing code to perform secondary analyses of PAMPer trial data and TBI biomarker data (Sperry et al., 2018, NEJM)
## Project Support: NIH T32 in Trauma and Sepsis, University of Pittsburgh Department of Surgery
## PI: Dr. Jason Sperry
## Written by Danielle Gruen, October 2019
## Last updated: January 2020

############################## Workspace Setup ##############################
## Set Up Workspace
    rm(list=ls()) # Clear all variables; delete all the objects in the workspace and start with a clean slate
    setwd("/Users/dgruen/Documents/Academia/Pitt_T32/projects/pamper/data/") # set the working directory
    set.seed(386) # make code reproducible

## Dependencies
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
    library(Amelia)
    library(mlbench)
    library(gee)
    library(geepack)

############################## File Setup ##############################
############################## Clean Data ##############################
############################## Data Inspection ##############################
############################## Main Tables and Figures ##############################   
## Fig: CONSORT diagram
    df <- df_master_filtered_wide
    
    tally_df <- df %>%
      group_by(FFP, TBI) %>%   # group by Biomarker
      tally()
    tally_df
    
## Table: Characteristics of TBI vs. No TBI
    # note: this was cross/spot checked with SPSS 
    
    df <- df_master_filtered_wide
    # dput(df)
    
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    myVars <- c("age", "gender",  "race",
                
                "initial_GCS", "initial_GCS_8", "total_gcs_score", "iss",
                "PH_sbp_70", "TBI", "chest_injury", "abdominal_injury", 
                "extremity_injury", "spinal_cord_injury", "any_blunt", "blunt_mechanism",
                "any_penetrating", "penetrating_mechanism",
                "severe_head", "ais_head", "ais_chest", "ais_abdomen", "ais_extremity", "vitals_hr", "vitals_sbp",
                
                "PH_intubation", "PH_CPR", "plasma_on_helicopter", "PH_crystalloid",
                "PH_prbc", "PH_blood", "PH_time", "PH_time_high", "transfer",
                
                "INR", "alpha", "K", "MA", "LY30", "hyperfibrinolysis",
                "transfusion_24h", "prbc_24h", "plasma_24h", "lactate_val",
                "platelets_24h", "crystalloid_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h",
                "ex_lap_24", "crani_24", "IR_V_E_24", "ortho_24", "proc_other_24",
                
                "alive_at_30", "mortality_24h", "ed_coagulopathy", 
                "MOF", "ALI", "NI", "icu_los", "hospital_los", "mech_vent_days",
                
                "GFAP_0", "GFAP_24", "GFAP_72", "UCH_L1_0", "UCH_L1_24", "UCH_L1_72")
    
    catVars <- c("gender","race", "alive_at_30", "TBI", "severe_head", "ed_coagulopathy", "chest_injury",
                 "abdominal_injury", "extremity_injury", "spinal_cord_injury",
                 "ALI", "NI", "MOF", "hyperfibrinolysis", "PH_intubation", "PH_sbp_70",
                 "any_blunt", "blunt_mechanism", "any_penetrating", "penetrating_mechanism",
                 "PH_CPR", "initial_GCS_8", "plasma_on_helicopter", "transfer",
                 "PH_blood", "ed_trali", "mortality_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h",
                 "ex_lap_24", "crani_24", "IR_V_E_24", "ortho_24", "proc_other_24", "PH_time_high"
                 )
    
    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "TBI", factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    # xtable(tab1)
    View(tab1)
    
## Table: Within TBI, Characteristics of prehospital plasma vs. standard care
    df <- df_master_filtered_wide
    df <- filter(df, TBI==1)

    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    myVars <- c("age", "gender", "race",
                
                "initial_GCS", "initial_GCS_8", "total_gcs_score", "iss",
                "PH_sbp_70", "TBI", "chest_injury", "abdominal_injury", 
                "extremity_injury", "spinal_cord_injury", "any_blunt", "blunt_mechanism", 
                "any_penetrating", "penetrating_mechanism",
                "severe_head",  "ais_head", "ais_chest", "ais_abdomen", "ais_extremity", "vitals_hr", "vitals_sbp",
                
                "PH_intubation", "PH_CPR", "plasma_on_helicopter", "PH_crystalloid",
                "PH_prbc", "PH_blood", "PH_time", "PH_time_high", "transfer",
                
                "INR", "alpha", "K", "MA", "LY30", "hyperfibrinolysis",
                "transfusion_24h", "prbc_24h", "plasma_24h", "lactate_val",
                "platelets_24h", "crystalloid_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h",
                "ex_lap_24", "crani_24", "IR_V_E_24", "ortho_24", "proc_other_24",
                
                "alive_at_30", "mortality_24h", "ed_coagulopathy", 
                "MOF", "ALI", "NI", "icu_los", "hospital_los", "mech_vent_days",
                
                "GFAP_0", "GFAP_24", "GFAP_72", "UCH_L1_0", "UCH_L1_24", "UCH_L1_72")
    
    catVars <- c("gender", "race", "alive_at_30", "TBI", "severe_head", "ed_coagulopathy", "chest_injury",
                 "abdominal_injury", "extremity_injury", "spinal_cord_injury",
                 "ALI", "NI", "MOF", "hyperfibrinolysis", "PH_intubation", "PH_sbp_70",
                 "any_blunt", "blunt_mechanism",  "any_penetrating", "penetrating_mechanism",
                 "PH_CPR", "initial_GCS_8", "plasma_on_helicopter", "transfer",
                 "PH_blood", "ed_trali", "mortality_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h",
                 "ex_lap_24", "crani_24", "IR_V_E_24", "ortho_24", "proc_other_24", "PH_time_high"
                 )
    
    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "FFP", factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, exact = "crani_24", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    # xtable(tab1)
    View(tab1)

## Analysis: Generalized Estimating Equations (GEE) (account for site clusters)
    ## Using a generalized estimating equations model to account for trial cluster effects and multiple confounders, 
    ## and testing an interaction between TBI and randomization group (plasma vs control):
  
    ## Data frame
    df <- df_master_filtered_wide
    df <- filter(df, alive_at_30==1 | alive_at_30==2) # filter to only keep patients with outcomes
    df$ais_head <- as.numeric(as.character(df$ais_head))
    
    ## get rid of any missing data in the fields used to run code
    df <- select(df, stnum, t_30d_censor_h, FFP, TBI, 
                 iss, MOI, SiteID, PH_time_high,
                 PH_prbc, PH_crystalloid, PH_blood,
                 PH_sbp_70, ais_head)
    
    df <- na.omit(df)

    fit.gee <- geeglm(t_30d_censor_h ~ TBI * FFP + iss + 
                        PH_prbc + PH_crystalloid + PH_blood + 
                        MOI + PH_sbp_70 + ais_head + PH_time_high,
                      family = binomial,
                      data = df, 
                      id = SiteID,
                      corstr = "exchangeable")
    summary(fit.gee)
    
    ## Conclusion: after adjustment for multiple confounders and site clustering, 
    ## TBI significantly modified the effect of plasma vs control on death
    ## with the benefit effect of plasma seen only in TBI trauma patients. 	
    ## Can split for K-M survival analysis

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
    ## COX within TBI
    df <- df_master_filtered_wide
    # dput(df)
    df <- filter(df, TBI==1)
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    res.cox <- coxph(Surv(t_30d_mort_h, t_30d_censor_h) ~ FFP + age +
                       transfer + ais_head + PH_crystalloid + PH_blood + PH_prbc +
                       iss + PH_intubation,
                     data =  df)
    
    summary(res.cox)
    
    vif(res.cox) # to test multicolinearity
    
    test.ph <- cox.zph(res.cox)
    test.ph
    
    # From the output above, the test is not statistically significant for each of the covariates,
    # and the global test is also not statistically significant.
    # Therefore, we can assume the proportional hazards.
    
    ## COX within NO TBI
    df <- df_master_filtered_wide
    # dput(df)
    df <- filter(df, TBI==2)
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    res.cox <- coxph(Surv(t_30d_mort_h, t_30d_censor_h) ~ FFP + age +
                       transfer + PH_crystalloid + PH_blood + PH_prbc +
                       iss + PH_intubation,
                     data =  df)
    
    summary(res.cox)
    
    test.ph <- cox.zph(res.cox)
    test.ph
    

## TBI biomarkers / SC vs PP
    df <- df_master_filtered_wide
    df <- filter(df, !is.na(TBI))
    df <- filter(df, !is.na(FFP))
    
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    myVars <- c("GFAP_0", "GFAP_24", "GFAP_72", "UCH_L1_0", "UCH_L1_24", "UCH_L1_72")
    
    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "TBI")
    tab1 <- print(tab1, nonnormal = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    # xtable(tab1)
    View(tab1)
    
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

    
## Figure: Forest Plot
    ## Generate Hazard Ratios: head injury severity
    df <- df_master_filtered_wide
    df$head_groups_dsg <- NA
    
    ## Assign head injury severity groups
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    df$head_groups_dsg[which(df$TBI == 2)] = 0
    df$head_groups_dsg[which(df$TBI == 1 & df$ais_head<=3)] = 1
    df$head_groups_dsg[which(df$TBI == 1 & (df$ais_head>3))] = 2
    df$head_groups_dsg <- as.factor(df$head_groups_dsg) # make it a factor
    df_master_filtered_wide_head_groups <- df
    
    df <- df_master_filtered_wide_head_groups
    df <- filter(df, !is.na(head_groups_dsg))
    df <- filter(df, !is.na(alive_at_30))

    df <- df_master_filtered_wide_head_groups
    df <- filter(df, head_groups_dsg==0)
    
    res.cox <- coxph(Surv(t_30d_mort_h, t_30d_censor_h) ~ FFP, data =  df)
    
    test <- summary(res.cox)
    coef(summary(res.cox))
    
    ggforest(res.cox, data = df)
    
    # Put this output into df
    
    
    ## Generate Hazard Ratios: polytrauma
    df <- df_master_filtered_wide
    df$trauma_groups_dsg <- NA
    
    ## Make polytrauma groups
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    df$trauma_groups_dsg[which(df$TBI == 2)] = 0
    df$trauma_groups_dsg[which(df$TBI == 1 & df$ais_abdomen<3 & df$ais_chest<3 & df$ais_extremity<3)] = 1
    df$trauma_groups_dsg[which(df$TBI == 1 & (df$ais_abdomen>=3 | df$ais_chest>=3 | df$ais_extremity>=3))] = 2
    df$trauma_groups_dsg <- as.factor(df$trauma_groups_dsg) # make it a factor
    
    df_master_filtered_wide_head <- df
    df <- df_master_filtered_wide_head
    df <- filter(df, !is.na(trauma_groups_dsg))
    df <- filter(df, !is.na(alive_at_30))
    
    df <- df_master_filtered_wide_head
    df <- filter(df, trauma_groups_dsg==0)
    
    res.cox <- coxph(Surv(t_30d_mort_h, t_30d_censor_h) ~ FFP, data =  df)
    
    test <- summary(res.cox)
    coef(summary(res.cox))
    
    ggforest(res.cox, data = df)
    
    # Put this output into df
    
    
    ## Final Plot
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
  
    
## Figure: Kaplan Meier early vs late plasma; Scene PP vs SC vs Transfer PP vs SC
    ## Transfer
    df <- df_master_filtered_wide
    df <- filter(df, transfer==1)
    
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
    
############################## Save ##############################   
## Optional: Save a table to .csv
    # write.table(biomarkers_long,
    #             file = "PAMPer_biomarkers.csv",
    #             append = FALSE, quote = TRUE, sep = ",",
    #             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
    #             col.names = TRUE, qmethod = c("escape", "double"),
    #             fileEncoding = "")
 