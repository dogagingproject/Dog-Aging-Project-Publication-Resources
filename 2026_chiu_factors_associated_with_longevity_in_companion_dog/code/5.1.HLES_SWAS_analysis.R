library(haven)
library(lubridate)
library(survival)
library(flexsurv)
library(survminer)
library(dplyr)
library(viridisLite)
library(stringr)
source("ggforest2.R")
datfolder <- "../data"
resultsfolder <- "results"
figfolder <- "figures"
load(file.path(datfolder,"SurvivalData.RData"))

###### Screening of variables
codebook <- read.csv(file.path(datfolder,"DAP_2024_CODEBOOK_v1.0.csv"))
codebook$Variable[codebook$Variable=="ss_household_dog_count"] <- 
  "od_household_dog_count" # Make part of owner contact
codebook$Variable[codebook$Variable=="ss_vet_frequency"] <- 
  "hs_vet_frequency" # Make part of health status
vargroups <- c(dd="Dog Demographics",
               oc="Owner Contact",
               pa="Physical Activity",
               de="Dog Environment",
               db="Dog Behavior",
               df="Diet Survey",
               dt="Comprehensive Diet Survey",
               mp="Medicines and Preventatives",
               hs="Health Status",
               od="Owner Demographics",
               fs="Additional Studies"
)
vars.to.keep <- codebook %>% 
  filter((DataFile == "HLES_dog_owner") &
           (Values != "YYYY-MM-DD") &
           (Values != "Text")) %>%
  filter(str_detect(Variable,
                    "dd_insurance|dd_activities|dd_alternate|pa_|de_|db_|df_|dt_|mp_|ss_|hs_|fs_|od_|oc_household")
  )
vars.to.keep$variable_group <- str_split(vars.to.keep$Variable,"_",simplify = TRUE)[,1]
vars.to.keep$variable_group_name <- vargroups[vars.to.keep$variable_group]
vars.to.keep$variable_group <- factor(vars.to.keep$variable_group,
                                      levels=unique(vars.to.keep$variable_group))
vars.to.keep$variable_group_name <- factor(vars.to.keep$variable_group_name,
                                           levels=unique(vars.to.keep$variable_group_name))

vars.to.keep$num.values <- 1+str_count(vars.to.keep$Values,"\\|")
vars.to.keep <- subset(vars.to.keep,num.values <= 51) # remove variables with many factors
vars.to.keep$make.factor <- TRUE
vars.to.keep$make.factor[vars.to.keep$Values=="Numeric"] <- FALSE
vars.to.keep$make.factor[str_detect(vars.to.keep$SurveyText,"Number|number|Duration|Fraction|Whole")] <- FALSE
vars.to.keep$make.factor[str_detect(vars.to.keep$SurveyText,"How many|how many") &
                           vars.to.keep$num.values > 3] <- FALSE
write.csv(vars.to.keep,file.path(resultsfolder,"Supp.HLES_Cox-variables.csv"),row.names = FALSE)

###### Extract variables from HLES ### 
load(file.path(datfolder,"DAP_2024_HLES_dog_owner_v1.0.RData"))
names(HLES_dog_owner)[names(HLES_dog_owner)=="ss_household_dog_count"] <- 
  "od_household_dog_count" # Make part of owner demographics
names(HLES_dog_owner)[names(HLES_dog_owner)=="ss_vet_frequency"] <- 
  "hs_vet_frequency" # Make part of health status
d <- HLES_dog_owner[,vars.to.keep$Variable]
vars.to.keep$num.na <- 0
for (j in 1:ncol(d)) {
  if (vars.to.keep$make.factor[j]) {
    labels.now <- trimws(strsplit(vars.to.keep$ValueLabels[j],"\\|")[[1]])
    names(labels.now) <- trimws(strsplit(vars.to.keep$Values[j],"\\|")[[1]])
    if (sum(duplicated(labels.now))>0) {
      labels.now <- paste(labels.now,names(labels.now),sep="-")
    }
    d[[j]] <- factor(labels.now[trimws(d[[j]])],
                     levels=as.character(labels.now))
    vars.to.keep$num.na[j]<-sum(is.na(d[[j]]))
  } else {
    d[[j]] <- as.numeric(d[[j]])
    # convert 97 to 0.5 for these variables (text is "<1")
    if(str_detect(vars.to.keep$Variable[j],"How many|how many|Number|number")) {
      d[[j]][d[[j]]==97] <- 0.5
    } 
    # convert 99 to NA for these variables (text is "Don't know")
    if(vars.to.keep$Values[j]!="Numeric") {
      if(str_detect(vars.to.keep$ValueLabels[j],"Don't know")) {
        d[[j]][d[[j]]==99] <- NA
      }
    }
    vars.to.keep$num.na[j]<-sum(is.na(d[[j]]))
  }
}
d <- cbind(data.frame(dog_id=HLES_dog_owner$dog_id),d)
d <- subset(d, dog_id %in% survivalData$dog_id)
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')
remove(HLES_dog_owner)

###### Cox modeling for each variable
pvals <- cbind(vars.to.keep,
               data.frame(Type=NA,Ref=NA,score.pval=NA,zph.pval=NA,score.pval.xrand=NA))
for (j in 1:nrow(pvals)) {
  var.now <- as.character(vars.to.keep$Variable[j])
  print(paste(j,var.now,"----------"))
  # Join with survival data
  data.tmp <- left_join(survivalData,d[,c("dog_id",var.now)])
  names(data.tmp)[ncol(data.tmp)]<-"x"
  set.seed(3.14159)
  data.tmp$xrand <- sample(data.tmp$x)
  pvals$Type[j] <- class(data.tmp$x)
  if (class(data.tmp$x) == "factor") { # re-order factors to reference is largest category
    orig_levels <- levels(data.tmp$x)
    maxlevel <- names(which.max(table(data.tmp$x)))
    pvals$Ref[j] <- maxlevel
    data.tmp$x <- factor(as.character(data.tmp$x),
                         levels=c(orig_levels[orig_levels==maxlevel],orig_levels[orig_levels!=maxlevel])
    )
  }
  try( {
    remove(cox.tmp); # clean up
    cox.tmp <- coxph(surv ~ x + strata(Size_Class_at_HLES,Breed_Class,Sex),
                     data=data.tmp);
    pvals$score.pval[j] <- summary(cox.tmp)$sctest["pvalue"];
    if (pvals$score.pval[j] < 1) {
      pvals$zph.pval[j] <- cox.zph(cox.tmp)$table["x","p"]
    }
    remove(cox.tmp.xrand); # clean up
    cox.tmp.xrand <- coxph(surv ~ xrand + strata(Size_Class_at_HLES,Breed_Class,Sex),
                           data=data.tmp);
    pvals$score.pval.xrand[j] <- summary(cox.tmp.xrand)$sctest["pvalue"];
  }
  )
}

# Use p.adjust to do FDR using BY approach due to potential dependency

pvals$Variable <- factor(pvals$Variable,levels=pvals$Variable)

pvals$score.pval.adjust <- p.adjust(pvals$score.pval,method="BY")

# Randomized adjusted p-values are all = 1
pvals$score.pval.xrand.adjust <- p.adjust(pvals$score.pval.xrand,method="BY")
pvals$zph.pval.adjust <- p.adjust(pvals$zph.pval,method="BY")

write.csv(pvals,file.path(resultsfolder,"Supp.HLES_Cox_pvals.csv"))

variable_group_tab <- table(pvals$variable_group_name)
variable_group_tab.df <- data.frame(xmin=1+c(0,as.numeric(cumsum(variable_group_tab)[-length(variable_group_tab)])),
                                    xmax=as.numeric(cumsum(variable_group_tab)),
                                    xmid=as.numeric(cumsum(variable_group_tab)-variable_group_tab/2),
                                    variable_group_name=unique(pvals$variable_group_name))

vars.to.plot <- subset(pvals, score.pval.adjust < 0.05)
vars.to.plot <- vars.to.plot[order(vars.to.plot$score.pval.adjust),]
cox.coef <- data.frame()
for (j in 1:nrow(vars.to.plot)) {
  var.now <- as.character(vars.to.plot$Variable[j])
  type.now <- vars.to.plot$Type[j]
  print(paste(j,type.now,var.now,"----------"))
  data.tmp <- left_join(survivalData,d[,c("dog_id",var.now)])
  names(data.tmp)[ncol(data.tmp)]<-"x"
  if (class(data.tmp$x) == "factor") { # re-order factors to reference is largest category
    orig_levels <- levels(data.tmp$x)
    maxlevel <- names(which.max(table(data.tmp$x)))
    data.tmp$x <- factor(as.character(data.tmp$x),
                         levels=c(orig_levels[orig_levels==maxlevel],orig_levels[orig_levels!=maxlevel])
    )
  }
  cox.tmp <- coxph(surv ~ x + strata(Size_Class_at_HLES,Breed_Class,Sex),
                   data=data.tmp);
  cox.tmp.sum <- summary(cox.tmp)
  cox.coef.tmp <- cbind(data.frame(Variable=var.now,Type=type.now,y=rownames(cox.tmp.sum$coefficients)),
                        cbind(as.data.frame(cox.tmp.sum$coefficients),as.data.frame(cox.tmp.sum$conf.int)))
  cox.coef <- rbind(cox.coef,cox.coef.tmp[,!duplicated(names(cox.coef.tmp))])
}


cox.coef$Variable.y <- paste(cox.coef$Variable,cox.coef$y,sep="|")
cox.coef$Variable.y <- gsub("\\|x","\\|",cox.coef$Variable.y)
cox.coef$Variable.y[!(cox.coef$Type %in% c("factor","character"))] <- 
  gsub("\\|","",cox.coef$Variable.y[!(cox.coef$Type %in% c("factor","character"))])
cox.coef <- cox.coef[order(cox.coef$Variable.y),]
cox.coef <- left_join(cox.coef,pvals)
write.csv(cox.coef,file.path(resultsfolder,"Supp.HLES_Cox_signif_results.csv"))

save(cox.coef,vars.to.keep,vars.to.plot,d,pvals,vargroups,variable_group_tab,
     variable_group_tab.df,
     file=file.path(resultsfolder,"Supp.HLES_Cox.Rdata"))

###### Sensitivity analysis - Mature Adult LifeStage only

survivalData.MAdult <- subset(survivalData,LifeStage_Class_at_HLES == "Mature Adult")

## Survival for Mature Adult LifeStage only
surv.MAdult <- Surv(time=survivalData.MAdult$first.age,
                    time2=survivalData.MAdult$last.age,
                    event=survivalData.MAdult$event, 
                    type="counting")

## Cox modeling for each variable
pvals.MAdult <- cbind(vars.to.keep,
                      data.frame(Type=NA,Ref=NA,score.pval=NA,zph.pval=NA,score.pval.xrand=NA))
for (j in 1:nrow(pvals.MAdult)) {
  var.now <- as.character(vars.to.keep$Variable[j])
  print(paste(j,var.now,"----------"))
  # Join with survival data
  data.tmp <- left_join(survivalData.MAdult,d[,c("dog_id",var.now)])
  names(data.tmp)[ncol(data.tmp)]<-"x"
  if (class(data.tmp$x) == "factor") { # re-order factors to reference is largest category
    orig_levels <- levels(data.tmp$x)
    maxlevel <- names(which.max(table(data.tmp$x)))
    pvals.MAdult$Ref[j] <- maxlevel
    data.tmp$x <- factor(as.character(data.tmp$x),
                         levels=c(orig_levels[orig_levels==maxlevel],orig_levels[orig_levels!=maxlevel])
    )
  }
  set.seed(3.14159)
  data.tmp$xrand <- sample(data.tmp$x)
  pvals.MAdult$Type[j] <- class(data.tmp$x)
  try( {
    remove(cox.tmp); # clean up
    cox.tmp <- coxph(surv.MAdult ~ x + strata(Size_Class_at_HLES,Breed_Class,Sex),
                     data=data.tmp);
    pvals.MAdult$score.pval[j] <- summary(cox.tmp)$sctest["pvalue"];
    if (pvals.MAdult$score.pval[j] < 1) {
      pvals.MAdult$zph.pval[j] <- cox.zph(cox.tmp)$table["x","p"]
    }
    remove(cox.tmp.xrand); # clean up
    cox.tmp.xrand <- coxph(surv.MAdult ~ xrand + strata(Size_Class_at_HLES,Breed_Class,Sex),
                           data=data.tmp);
    pvals.MAdult$score.pval.xrand[j] <- summary(cox.tmp.xrand)$sctest["pvalue"];
  }
  )
}

# Use p.adjust to do FDR using BY approach due to potential dependency

pvals.MAdult$Variable <- factor(pvals.MAdult$Variable,levels=pvals.MAdult$Variable)

pvals.MAdult$score.pval.adjust <- p.adjust(pvals.MAdult$score.pval,method="BY")

# No false discoveries
pvals.MAdult$score.pval.xrand.adjust <- p.adjust(pvals.MAdult$score.pval.xrand,method="BY")
pvals.MAdult$zph.pval.adjust <- p.adjust(pvals.MAdult$zph.pval,method="BY")

write.csv(pvals.MAdult,file.path(resultsfolder,"Supp.HLES_Cox_pvals.MAdult.csv"))

variable_group_tab <- table(pvals.MAdult$variable_group_name)
variable_group_tab.df <- data.frame(xmin=1+c(0,as.numeric(cumsum(variable_group_tab)[-length(variable_group_tab)])),
                                    xmax=as.numeric(cumsum(variable_group_tab)),
                                    xmid=as.numeric(cumsum(variable_group_tab)-variable_group_tab/2),
                                    variable_group_name=unique(pvals.MAdult$variable_group_name))

vars.to.plot.MAdult <- subset(pvals.MAdult, score.pval.adjust < 0.05)
vars.to.plot.MAdult <- vars.to.plot.MAdult[order(vars.to.plot.MAdult$score.pval.adjust),]
cox.coef.MAdult <- data.frame()
for (j in 1:nrow(vars.to.plot.MAdult)) {
  var.now <- as.character(vars.to.plot.MAdult$Variable[j])
  type.now <- vars.to.plot.MAdult$Type[j]
  print(paste(j,type.now,var.now,"----------"))
  data.tmp <- left_join(survivalData.MAdult,d[,c("dog_id",var.now)])
  names(data.tmp)[ncol(data.tmp)]<-"x"
  if (class(data.tmp$x) == "factor") { # re-order factors to reference is largest category
    orig_levels <- levels(data.tmp$x)
    maxlevel <- names(which.max(table(data.tmp$x)))
    data.tmp$x <- factor(as.character(data.tmp$x),
                         levels=c(orig_levels[orig_levels==maxlevel],orig_levels[orig_levels!=maxlevel])
    )
  }
  cox.tmp <- coxph(surv.MAdult ~ x + strata(Size_Class_at_HLES,Breed_Class,Sex),
                   data=data.tmp);
  cox.tmp.sum <- summary(cox.tmp)
  cox.coef.tmp <- cbind(data.frame(Variable=var.now,Type=type.now,y=rownames(cox.tmp.sum$coefficients)),
                        cbind(as.data.frame(cox.tmp.sum$coefficients),as.data.frame(cox.tmp.sum$conf.int)))
  cox.coef.MAdult <- rbind(cox.coef.MAdult,cox.coef.tmp[,!duplicated(names(cox.coef.tmp))])
}

cox.coef.MAdult$Variable.y <- paste(cox.coef.MAdult$Variable,cox.coef.MAdult$y,sep="|")
cox.coef.MAdult$Variable.y <- gsub("\\|x","\\|",cox.coef.MAdult$Variable.y)
cox.coef.MAdult$Variable.y[!(cox.coef.MAdult$Type %in% c("factor","character"))] <- 
  gsub("\\|","",cox.coef.MAdult$Variable.y[!(cox.coef.MAdult$Type %in% c("factor","character"))])
cox.coef.MAdult <- cox.coef.MAdult[order(cox.coef.MAdult$Variable.y),]
cox.coef.MAdult <- left_join(cox.coef.MAdult,pvals.MAdult)
write.csv(cox.coef.MAdult,file.path(resultsfolder,"Supp.HLES_Cox_signif_results.MAdult.csv"))

save(cox.coef.MAdult,vars.to.keep,vars.to.plot.MAdult,d,pvals.MAdult,vargroups,variable_group_tab,
     variable_group_tab.df,
     file=file.path(resultsfolder,"Supp.HLES_Cox.MAdult.Rdata"))

###### Sensitivity analysis - top16breeds

top16breeds <-head(sort(table(subset(survivalData,
                                     Breed_Class=="AKC-Recognized Breed")$Breed),
                        decreasing = TRUE),16)
survivalData.top16 <- subset(survivalData,Breed_Class=="AKC-Recognized Breed" & 
                               Breed %in% names(top16breeds))
survivalData.top16$Breed <- factor(survivalData.top16$Breed,
                                   levels=names(top16breeds))
print(top16breeds)

## Survival for 16 Individual Breeds
surv.top16 <- Surv(time=survivalData.top16$first.age,
                   time2=survivalData.top16$last.age,
                   event=survivalData.top16$event, 
                   type="counting")

## Cox modeling for each variable
pvals.top16 <- cbind(vars.to.keep,
               data.frame(Type=NA,Ref=NA,score.pval=NA,zph.pval=NA,score.pval.xrand=NA))
for (j in 1:nrow(pvals.top16)) {
  var.now <- as.character(vars.to.keep$Variable[j])
  print(paste(j,var.now,"----------"))
  # Join with survival data
  data.tmp <- left_join(survivalData.top16,d[,c("dog_id",var.now)])
  names(data.tmp)[ncol(data.tmp)]<-"x"
  if (class(data.tmp$x) == "factor") { # re-order factors to reference is largest category
    orig_levels <- levels(data.tmp$x)
    maxlevel <- names(which.max(table(data.tmp$x)))
    pvals.top16$Ref[j] <- maxlevel
    data.tmp$x <- factor(as.character(data.tmp$x),
                         levels=c(orig_levels[orig_levels==maxlevel],orig_levels[orig_levels!=maxlevel])
    )
  }
  set.seed(3.14159)
  data.tmp$xrand <- sample(data.tmp$x)
  pvals.top16$Type[j] <- class(data.tmp$x)
  try( {
    remove(cox.tmp); # clean up
    cox.tmp <- coxph(surv.top16 ~ x + strata(Breed,Sex),
                     data=data.tmp);
    pvals.top16$score.pval[j] <- summary(cox.tmp)$sctest["pvalue"];
    if (pvals.top16$score.pval[j] < 1) {
      pvals.top16$zph.pval[j] <- cox.zph(cox.tmp)$table["x","p"]
    }
    remove(cox.tmp.xrand); # clean up
    cox.tmp.xrand <- coxph(surv.top16 ~ xrand + strata(Breed,Sex),
                           data=data.tmp);
    pvals.top16$score.pval.xrand[j] <- summary(cox.tmp.xrand)$sctest["pvalue"];
  }
  )
}

# Use p.adjust to do FDR using BY approach due to potential dependency

pvals.top16$Variable <- factor(pvals.top16$Variable,levels=pvals.top16$Variable)

pvals.top16$score.pval.adjust <- p.adjust(pvals.top16$score.pval,method="BY")

# Only 1 randomized p-value < 0.05, and only 3 < 1, so very small false discovery rate
pvals.top16$score.pval.xrand.adjust <- p.adjust(pvals.top16$score.pval.xrand,method="BY")
pvals.top16$zph.pval.adjust <- p.adjust(pvals.top16$zph.pval,method="BY")

write.csv(pvals.top16,file.path(resultsfolder,"Supp.HLES_Cox_pvals.top16.csv"))

variable_group_tab <- table(pvals.top16$variable_group_name)
variable_group_tab.df <- data.frame(xmin=1+c(0,as.numeric(cumsum(variable_group_tab)[-length(variable_group_tab)])),
                                    xmax=as.numeric(cumsum(variable_group_tab)),
                                    xmid=as.numeric(cumsum(variable_group_tab)-variable_group_tab/2),
                                    variable_group_name=unique(pvals.top16$variable_group_name))

vars.to.plot.top16 <- subset(pvals.top16, score.pval.adjust < 0.05)
vars.to.plot.top16 <- vars.to.plot.top16[order(vars.to.plot.top16$score.pval.adjust),]
cox.coef.top16 <- data.frame()
for (j in 1:nrow(vars.to.plot.top16)) {
  var.now <- as.character(vars.to.plot.top16$Variable[j])
  type.now <- vars.to.plot.top16$Type[j]
  print(paste(j,type.now,var.now,"----------"))
  data.tmp <- left_join(survivalData.top16,d[,c("dog_id",var.now)])
  names(data.tmp)[ncol(data.tmp)]<-"x"
  if (class(data.tmp$x) == "factor") { # re-order factors to reference is largest category
    orig_levels <- levels(data.tmp$x)
    maxlevel <- names(which.max(table(data.tmp$x)))
    pvals.top16$Ref[j] <- maxlevel
    data.tmp$x <- factor(as.character(data.tmp$x),
                         levels=c(orig_levels[orig_levels==maxlevel],orig_levels[orig_levels!=maxlevel])
    )
  }
  cox.tmp <- coxph(surv.top16 ~ x + strata(Breed,Sex),
                   data=data.tmp);
  cox.tmp.sum <- summary(cox.tmp)
  cox.coef.tmp <- cbind(data.frame(Variable=var.now,Type=type.now,y=rownames(cox.tmp.sum$coefficients)),
                        cbind(as.data.frame(cox.tmp.sum$coefficients),as.data.frame(cox.tmp.sum$conf.int)))
  cox.coef.top16 <- rbind(cox.coef.top16,cox.coef.tmp[,!duplicated(names(cox.coef.tmp))])
}

cox.coef.top16$Variable.y <- paste(cox.coef.top16$Variable,cox.coef.top16$y,sep="|")
cox.coef.top16$Variable.y <- gsub("\\|x","\\|",cox.coef.top16$Variable.y)
cox.coef.top16$Variable.y[!(cox.coef.top16$Type %in% c("factor","character"))] <- 
  gsub("\\|","",cox.coef.top16$Variable.y[!(cox.coef.top16$Type %in% c("factor","character"))])
cox.coef.top16 <- cox.coef.top16[order(cox.coef.top16$Variable.y),]
cox.coef.top16 <- left_join(cox.coef.top16,pvals.top16)
write.csv(cox.coef.top16,file.path(resultsfolder,"Supp.HLES_Cox_signif_results.top16.csv"))

save(cox.coef.top16,vars.to.keep,vars.to.plot.top16,d,pvals.top16,vargroups,variable_group_tab,
     variable_group_tab.df,
     file=file.path(resultsfolder,"Supp.HLES_Cox.top16.Rdata"))

