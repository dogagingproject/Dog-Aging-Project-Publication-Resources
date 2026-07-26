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

# Second Alternative analysis removing non-AKC single breed dogs

survivalData <- subset(survivalData,!(Breed_Status=="Purebred" & Breed_Class=="Non-AKC-Recognized or Mixed Breed"))
survivalData$Breed_Class <- gsub("Non-AKC-Recognized or ","",survivalData$Breed_Class)
survivalData$Breed_Class <- factor(survivalData$Breed_Class,levels=
                                     c("Mixed Breed","AKC-Recognized Breed"))
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')

fit.c.s <- survfit(surv ~ Size_Class_at_HLES + Breed_Class + Sex,data=survivalData)

fit.c.s.sum <- summary(fit.c.s)
fit.c.s.sum.df <- as.data.frame(fit.c.s.sum$table)
fit.c.s.quant <- quantile(fit.c.s)
fit.c.s.sum.df$q25 <- fit.c.s.quant$quantile[,1]
fit.c.s.sum.df$q25.lcl <- fit.c.s.quant$lower[,1]
fit.c.s.sum.df$q25.ucl <- fit.c.s.quant$upper[,1]
fit.c.s.sum.df$q75 <- fit.c.s.quant$quantile[,3]
fit.c.s.sum.df$q75.lcl <- fit.c.s.quant$lower[,3]
fit.c.s.sum.df$q75.ucl <- fit.c.s.quant$upper[,3]
fit.c.s.sum.df$stratum <- rownames(fit.c.s.sum$table)
strata<-trimws(unlist(strsplit(unlist(
  strsplit(rownames(fit.c.s.sum$table),",")),"=")))
strata <- strata[seq(2,length(strata),2)]
strata <- matrix(strata,ncol=3,byrow=TRUE)
fit.c.s.sum.df$Size_Class_at_HLES <-
  factor(strata[,1],levels=levels(survivalData$Size_Class_at_HLES))
fit.c.s.sum.df$Breed_Class <- 
  factor(strata[,2],levels=levels(survivalData$Breed_Class))
fit.c.s.sum.df$Sex <- 
  factor(strata[,3],levels=levels(survivalData$Sex))

write.csv(fit.c.s.sum.df,file.path(resultsfolder,"Alt2Supp.Survival_DAP_demographics.SumStats.csv"))

fit.b.s.sum.df <- read.csv(file.path(resultsfolder,"Survival_DAP_demographics.SumStats.csv"))

plt.compare.mean <-
ggplot(subset(data.frame(x=fit.b.s.sum.df$rmean,
                  y=fit.c.s.sum.df$rmean,
                  Sex=fit.b.s.sum.df$Sex,
                  Size_Class_at_HLES=fit.b.s.sum.df$Size_Class_at_HLES,
                  Breed_Class=fit.b.s.sum.df$Breed_Class),
                  Breed_Class!="AKC-Recognized Breed"),aes(x,y))+
  geom_point(aes(color=Size_Class_at_HLES,shape=Sex),size=3) + geom_abline(slope=1,intercept=0) + 
  scale_color_viridis_d(option="turbo",end=0.85,direction=-1)+
  xlab("Mean survival yrs for Mixed Breed\n(including non-AKC single breed)") +
  ylab("Mean survival yrs for Mixed Breed\n(excluding non-AKC single breed)") +
  theme_bw()+theme(legend.position = "top",legend.box = "vertical")

## Size effect by breed*sex strata

cox.c.s <- coxph(surv ~ Size_Class_at_HLES + Breed_Class + Sex,data=survivalData)
cox.c.s.pvalues <- data.frame()
for (Breed_Class_now in levels(survivalData$Breed_Class)) {
  for (Sex_now in levels(survivalData$Sex)) {
    cox.c.tmp <- coxph(surv ~ Size_Class_at_HLES,data=survivalData,
                       subset=Breed_Class == Breed_Class_now &
                         Sex == Sex_now)
    cox.c.s.pvalues <- rbind(cox.c.s.pvalues,
                             data.frame(Breed_Class=Breed_Class_now,
                                        Sex=Sex_now,
                                        log.rank.pvalue=summary(cox.c.tmp)$sctest["pvalue"]))
  }
}
# Make factors
cox.c.s.pvalues$Breed_Class <- factor(cox.c.s.pvalues$Breed_Class,
                                      levels=levels(survivalData$Breed_Class))
cox.c.s.pvalues$Sex <- factor(cox.c.s.pvalues$Sex,
                              levels=levels(survivalData$Sex))
print(cox.c.s.pvalues)
write.csv(cox.c.s.pvalues,file.path(resultsfolder,"Alt2Supp.Survival_DAP_demographics_sizeeffect-pvals.csv"),
          row.names = FALSE)
# Reverse levels for figure
levels(survivalData$Size_Class_at_HLES) <- rev(levels(survivalData$Size_Class_at_HLES))
plt.c.km.s<-ggsurvplot_facet(fit.c.s,data=survivalData,
                             facet.by=c("Sex","Breed_Class"),
                             censor=FALSE,xlab="Age (yr)",
                             surv.median.line = "hv",short.panel.labs=TRUE
)+scale_color_viridis_d(option="turbo",end=0.85)+
  geom_text(aes(x=20,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.c.s.pvalues,hjust=0)+theme(legend.position = "bottom")

print(plt.c.km.s)
ggsave(file.path(figfolder,"Alt2Supp.Fig.Survival_DAP_demographics.KM.pdf"),
       plt.c.km.s,
       height=4.5,width=6,scale=1.2)
# Un-Reverse levels 
levels(survivalData$Size_Class_at_HLES) <- rev(levels(survivalData$Size_Class_at_HLES))

######

cox.c.breed.pvalues <- data.frame()
for (Size_Class_at_HLES_now in levels(survivalData$Size_Class_at_HLES)) {
  for (Sex_now in levels(survivalData$Sex)) {
    cox.c.tmp <- coxph(surv ~ Breed_Class,data=survivalData,
                       subset=Size_Class_at_HLES == Size_Class_at_HLES_now &
                         Sex == Sex_now)
    cox.c.breed.pvalues <- rbind(cox.c.breed.pvalues,
                                 data.frame(Size_Class_at_HLES=Size_Class_at_HLES_now,
                                            Sex=Sex_now,
                                            log.rank.pvalue=summary(cox.c.tmp)$sctest["pvalue"]))
  }
}
# Make factors
cox.c.breed.pvalues$Size_Class_at_HLES <- factor(cox.c.breed.pvalues$Size_Class_at_HLES,
                                                 levels=levels(survivalData$Size_Class_at_HLES))
cox.c.breed.pvalues$Sex <- factor(cox.c.breed.pvalues$Sex,
                                  levels=levels(survivalData$Sex))
print(cox.c.breed.pvalues)
write.csv(cox.c.breed.pvalues,
          file.path(resultsfolder,"Alt2Supp.Survival_DAP_demographics_breedeffect-pvals.csv"),
          row.names = FALSE)
plt.c.km.breed<-ggsurvplot_facet(fit.c.s,data=survivalData,
                                 facet.by=c("Sex","Size_Class_at_HLES"),
                                 censor=FALSE,xlab="Age (yr)",
                                 surv.median.line = "hv",short.panel.labs=TRUE
)+scale_color_viridis_d(option="magma",end = 0.8)+
  geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.c.breed.pvalues,hjust=0)
print(plt.c.km.breed)
ggsave(file.path(figfolder,"Alt2Supp.Fig.Survival_DAP_demographics.breed_effect.pdf"),
       plt.c.km.breed,
       height=4,width=7,scale=1.6)

cox.b.breed.pvalues <- read.csv(file.path(resultsfolder,"Survival_DAP_demographics_breedeffect-pvals.csv"))

plt.compare.breed.pvalues <-
  ggplot(data.frame(x=cox.b.breed.pvalues$log.rank.pvalue,
                           y=cox.c.breed.pvalues$log.rank.pvalue,
                           Sex=cox.b.breed.pvalues$Sex,
                           Size_Class_at_HLES=cox.b.breed.pvalues$Size_Class_at_HLES),
         aes(x,y))+
  geom_point(aes(color=Size_Class_at_HLES,shape=Sex),size=3) + 
  geom_abline(slope=1,intercept=0) + 
  geom_hline(yintercept=0.05,linetype="dotted")+
  geom_vline(xintercept=0.05,linetype="dotted")+
  scale_x_log10(limits=c(NA,1),breaks=10^seq(-12,0,2)) + 
  scale_y_log10(limits=c(NA,1),breaks=10^seq(-12,0,2))+
  scale_color_viridis_d(option="turbo",end=0.85,direction=-1)+
  xlab("Breed effect log rank p-value\n(including non-AKC single breed)") +
  ylab("Breed effect log rank p-value\n(excluding non-AKC single breed)") +
  theme_bw()+theme(legend.position = "top",legend.box = "vertical")

######

cox.c.sex.pvalues <- data.frame()
for (Size_Class_at_HLES_now in levels(survivalData$Size_Class_at_HLES)) {
  for (Breed_Class_now in levels(survivalData$Breed_Class)) {
    cox.c.tmp <- coxph(surv ~ Sex,data=survivalData,
                       subset=Size_Class_at_HLES == Size_Class_at_HLES_now &
                         Breed_Class == Breed_Class_now)
    cox.c.sex.pvalues <- rbind(cox.c.sex.pvalues,
                               data.frame(Size_Class_at_HLES=Size_Class_at_HLES_now,
                                          Breed_Class=Breed_Class_now,
                                          log.rank.pvalue=summary(cox.c.tmp)$sctest["pvalue"]))
  }
}
cox.c.sex.pvalues$Size_Class_at_HLES <- factor(cox.c.sex.pvalues$Size_Class_at_HLES,
                                               levels=levels(survivalData$Size_Class_at_HLES))
cox.c.sex.pvalues$Breed_Class <- factor(cox.c.sex.pvalues$Breed_Class,
                                        levels=levels(survivalData$Breed_Class))

print(cox.c.sex.pvalues)
write.csv(cox.c.sex.pvalues,
          file.path(resultsfolder,"Alt2Supp.Survival_DAP_demographics_sexeffect-pvals.csv"),
          row.names = FALSE)
plt.c.km.sex<-ggsurvplot_facet(fit.c.s,data=survivalData,
                               facet.by=c("Breed_Class","Size_Class_at_HLES"),
                               censor=FALSE,xlab="Age (yr)",
                               surv.median.line = "hv",short.panel.labs=TRUE
)+scale_color_viridis_d(option="magma",end = 0.8)+
  geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.c.sex.pvalues,hjust=0)

print(plt.c.km.sex)
ggsave(file.path(figfolder,"Alt2Supp.Fig.Survival_DAP_demographics.sex_effect.pdf"),
       plt.c.km.sex,
       height=4,width=7,scale=1.6)


#### 

plt.forest.median.c.s<-ggplot(fit.c.s.sum.df,aes(y=interaction(Sex,Breed_Class)))+
  geom_boxplot(aes(xmin=q25,xlower=q25,xmiddle=median,xupper=q75,xmax=q75,
                   color=Size_Class_at_HLES),fill=NA,width=0.5,stat="identity")+
  xlab("Median & IQR lifespan")+
  scale_color_viridis_d(option="turbo",end=0.85)+
  guides(color="none")+theme_bw()+theme(legend.position="bottom")+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  scale_x_continuous(minor_breaks=seq(5,20))+
  facet_wrap(~Size_Class_at_HLES ,ncol=1)
print(plt.forest.median.c.s)

plt.forest.mean.c.s<-ggplot(fit.c.s.sum.df)+
  geom_errorbarh(aes(xmin=rmean-1.96*`se(rmean)`,
                     xmax=rmean+1.96*`se(rmean)`,
                     y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES),height=0)+
  geom_point(aes(x=rmean,y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES))+
  xlab("Mean lifespan [CI]")+scale_color_viridis_d(option="turbo",end=0.85)+
  guides(color="none")+theme_bw()+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  scale_x_continuous(minor_breaks=seq(5,20))+
  facet_wrap(~Size_Class_at_HLES ,ncol=1)
print(plt.forest.mean.c.s)

ggsave(file.path(figfolder,"Alt2Supp.Fig.Survival_DAP_demographics.forest.median.IQR.pdf"),
       plt.forest.median.c.s)
ggsave(file.path(figfolder,"Alt2Supp.Fig.Survival_DAP_demographics.forest.mean.pdf"),
       plt.forest.mean.c.s)

ggsave(file.path(figfolder,"Alt2Supp.Fig.Compare_with_main_analysis.pdf"),
       ggarrange(plt.compare.mean,plt.compare.breed.pvalues,
                 common.legend = TRUE),height=4,width=7,scale=1.2)
