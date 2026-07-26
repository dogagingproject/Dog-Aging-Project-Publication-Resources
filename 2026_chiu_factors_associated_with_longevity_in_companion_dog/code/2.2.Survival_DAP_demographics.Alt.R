library(lubridate)
library(survival)
library(flexsurv)
library(survminer)
library(dplyr)
library(viridisLite)
datfolder <- "../data"
resultsfolder <- "results"
figfolder <- "figures"
load(file.path(datfolder,"SurvivalData.RData"))
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')

# Alternative analysis using Weight_Class_10KGBin_at_HLES

fit.a.s <- survfit(surv ~ Weight_Class_10KGBin_at_HLES + Breed_Class + Sex,data=survivalData)

fit.a.s.sum <- summary(fit.a.s)
fit.a.s.sum.df <- as.data.frame(fit.a.s.sum$table)
fit.a.s.quant <- quantile(fit.a.s)
fit.a.s.sum.df$q25 <- fit.a.s.quant$quantile[,1]
fit.a.s.sum.df$q25.lcl <- fit.a.s.quant$lower[,1]
fit.a.s.sum.df$q25.ucl <- fit.a.s.quant$upper[,1]
fit.a.s.sum.df$q75 <- fit.a.s.quant$quantile[,3]
fit.a.s.sum.df$q75.lcl <- fit.a.s.quant$lower[,3]
fit.a.s.sum.df$q75.ucl <- fit.a.s.quant$upper[,3]
fit.a.s.sum.df$stratum <- rownames(fit.a.s.sum$table)
strata<-trimws(unlist(strsplit(unlist(
  strsplit(rownames(fit.a.s.sum$table),",")),"=")))
strata <- strata[seq(2,length(strata),2)]
strata <- matrix(strata,ncol=3,byrow=TRUE)
fit.a.s.sum.df$Weight_Class_10KGBin_at_HLES <-
  factor(strata[,1],levels=levels(survivalData$Weight_Class_10KGBin_at_HLES))
fit.a.s.sum.df$Breed_Class <- 
  factor(strata[,2],levels=levels(survivalData$Breed_Class))
fit.a.s.sum.df$Sex <- 
  factor(strata[,3],levels=levels(survivalData$Sex))

write.csv(fit.a.s.sum.df,file.path(resultsfolder,"AltSupp.Survival_DAP_demographics.SumStats.csv"))

## Weight Class effect by breed*sex strata

cox.a.s <- coxph(surv ~ Weight_Class_10KGBin_at_HLES + Breed_Class + Sex,data=survivalData)
cox.a.s.pvalues <- data.frame()
for (Breed_Class_now in levels(survivalData$Breed_Class)) {
  for (Sex_now in levels(survivalData$Sex)) {
    cox.a.tmp <- coxph(surv ~ Weight_Class_10KGBin_at_HLES,data=survivalData,
                       subset=Breed_Class == Breed_Class_now &
                         Sex == Sex_now)
    cox.a.s.pvalues <- rbind(cox.a.s.pvalues,
                             data.frame(Breed_Class=Breed_Class_now,
                                        Sex=Sex_now,
                                        log.rank.pvalue=summary(cox.a.tmp)$sctest["pvalue"]))
  }
}
# Make factors
cox.a.s.pvalues$Breed_Class <- factor(cox.a.s.pvalues$Breed_Class,
                                      levels=levels(survivalData$Breed_Class))
cox.a.s.pvalues$Sex <- factor(cox.a.s.pvalues$Sex,
                              levels=levels(survivalData$Sex))
print(cox.a.s.pvalues)
write.csv(cox.a.s.pvalues,file.path(resultsfolder,"AltSupp.Survival_DAP_demographics_weighteffect-pvals.csv"),
          row.names = FALSE)
# Reverse levels for figure
levels(survivalData$Weight_Class_10KGBin_at_HLES) <- rev(levels(survivalData$Weight_Class_10KGBin_at_HLES))
plt.a.km.s<-ggsurvplot_facet(fit.a.s,data=survivalData,
                             facet.by=c("Sex","Breed_Class"),
                             censor=FALSE,xlab="Age (yr)",
                             surv.median.line = "hv",short.panel.labs=TRUE
)+scale_color_viridis_d(option="turbo",end=0.85)+
  geom_text(aes(x=20,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.a.s.pvalues,hjust=0)

print(plt.a.km.s)
ggsave(file.path(figfolder,"AltSupp.Fig.Survival_DAP_demographics.KM.pdf"),
       plt.a.km.s,
       height=4.5,width=6,scale=1.2)
# Un-Reverse levels 
levels(survivalData$Weight_Class_10KGBin_at_HLES) <- rev(levels(survivalData$Weight_Class_10KGBin_at_HLES))

######

cox.a.breed.pvalues <- data.frame()
for (Weight_Class_10KGBin_at_HLES_now in levels(survivalData$Weight_Class_10KGBin_at_HLES)) {
  for (Sex_now in levels(survivalData$Sex)) {
    cox.a.tmp <- coxph(surv ~ Breed_Class,data=survivalData,
                       subset=Weight_Class_10KGBin_at_HLES == Weight_Class_10KGBin_at_HLES_now &
                         Sex == Sex_now)
    cox.a.breed.pvalues <- rbind(cox.a.breed.pvalues,
                                 data.frame(Weight_Class_10KGBin_at_HLES=Weight_Class_10KGBin_at_HLES_now,
                                            Sex=Sex_now,
                                            log.rank.pvalue=summary(cox.a.tmp)$sctest["pvalue"]))
  }
}
# Make factors
cox.a.breed.pvalues$Weight_Class_10KGBin_at_HLES <- factor(cox.a.breed.pvalues$Weight_Class_10KGBin_at_HLES,
                                                 levels=levels(survivalData$Weight_Class_10KGBin_at_HLES))
cox.a.breed.pvalues$Sex <- factor(cox.a.breed.pvalues$Sex,
                                  levels=levels(survivalData$Sex))
print(cox.a.breed.pvalues)
write.csv(cox.a.breed.pvalues,
          file.path(resultsfolder,"AltSupp.Survival_DAP_demographics_breedeffect-pvals.csv"),
          row.names = FALSE)
plt.a.km.breed<-ggsurvplot_facet(fit.a.s,data=survivalData,
                                 facet.by=c("Sex","Weight_Class_10KGBin_at_HLES"),
                                 censor=FALSE,xlab="Age (yr)",
                                 surv.median.line = "hv",short.panel.labs=TRUE
)+scale_color_viridis_d(option="magma",end = 0.8)+
  geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.a.breed.pvalues,hjust=0)
print(plt.a.km.breed)
ggsave(file.path(figfolder,"AltSupp.Fig.Survival_DAP_demographics.breed_effect.pdf"),
       plt.a.km.breed,
       height=4,width=7,scale=1.6)

######

cox.a.sex.pvalues <- data.frame()
for (Weight_Class_10KGBin_at_HLES_now in levels(survivalData$Weight_Class_10KGBin_at_HLES)) {
  for (Breed_Class_now in levels(survivalData$Breed_Class)) {
    cox.a.tmp <- coxph(surv ~ Sex,data=survivalData,
                       subset=Weight_Class_10KGBin_at_HLES == Weight_Class_10KGBin_at_HLES_now &
                         Breed_Class == Breed_Class_now)
    cox.a.sex.pvalues <- rbind(cox.a.sex.pvalues,
                               data.frame(Weight_Class_10KGBin_at_HLES=Weight_Class_10KGBin_at_HLES_now,
                                          Breed_Class=Breed_Class_now,
                                          log.rank.pvalue=summary(cox.a.tmp)$sctest["pvalue"]))
  }
}
cox.a.sex.pvalues$Weight_Class_10KGBin_at_HLES <- factor(cox.a.sex.pvalues$Weight_Class_10KGBin_at_HLES,
                                               levels=levels(survivalData$Weight_Class_10KGBin_at_HLES))
cox.a.sex.pvalues$Breed_Class <- factor(cox.a.sex.pvalues$Breed_Class,
                                        levels=levels(survivalData$Breed_Class))

print(cox.a.sex.pvalues)
write.csv(cox.a.sex.pvalues,
          file.path(resultsfolder,"AltSupp.Survival_DAP_demographics_sexeffect-pvals.csv"),
          row.names = FALSE)
plt.a.km.sex<-ggsurvplot_facet(fit.a.s,data=survivalData,
                               facet.by=c("Breed_Class","Weight_Class_10KGBin_at_HLES"),
                               censor=FALSE,xlab="Age (yr)",
                               surv.median.line = "hv",short.panel.labs=TRUE
)+scale_color_viridis_d(option="magma",end = 0.8)+
  geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.a.sex.pvalues,hjust=0)

print(plt.a.km.sex)
ggsave(file.path(figfolder,"AltSupp.Fig.Survival_DAP_demographics.sex_effect.pdf"),
       plt.a.km.sex,
       height=4,width=7,scale=1.6)


#### 

plt.forest.median.a.s<-ggplot(fit.a.s.sum.df,aes(y=interaction(Sex,Breed_Class)))+
  geom_boxplot(aes(xmin=q25,xlower=q25,xmiddle=median,xupper=q75,xmax=q75,
                   color=Weight_Class_10KGBin_at_HLES),fill=NA,width=0.5,stat="identity")+
  xlab("Median & IQR lifespan")+
  scale_color_viridis_d(option="turbo",end=0.85)+
  guides(color="none")+theme_bw()+theme(legend.position="bottom")+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  scale_x_continuous(minor_breaks=seq(5,20))+
  facet_wrap(~Weight_Class_10KGBin_at_HLES,ncol=1)
print(plt.forest.median.a.s)


plt.forest.mean.a.s<-ggplot(fit.a.s.sum.df)+
  geom_errorbarh(aes(xmin=rmean-1.96*`se(rmean)`,
                     xmax=rmean+1.96*`se(rmean)`,
                     y=interaction(Sex,Breed_Class),color=Weight_Class_10KGBin_at_HLES),height=0)+
  geom_point(aes(x=rmean,y=interaction(Sex,Breed_Class),color=Weight_Class_10KGBin_at_HLES))+
  xlab("Mean lifespan [CI]")+scale_color_viridis_d(option="turbo",end=0.85)+
  guides(color="none")+theme_bw()+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  scale_x_continuous(minor_breaks=seq(5,20))+
  facet_wrap(~Weight_Class_10KGBin_at_HLES ,ncol=1)
print(plt.forest.mean.a.s)

ggsave(file.path(figfolder,"AltSupp.Fig.Survival_DAP_demographics.forest.median.IQR.pdf"),
       plt.forest.median.a.s)
ggsave(file.path(figfolder,"AltSupp.Fig.Survival_DAP_demographics.forest.mean.pdf"),
       plt.forest.mean.a.s)
