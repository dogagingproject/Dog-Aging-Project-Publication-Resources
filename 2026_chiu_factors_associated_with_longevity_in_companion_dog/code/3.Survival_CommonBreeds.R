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
top16breeds <-head(sort(table(subset(survivalData,
                                     Breed_Class=="AKC-Recognized Breed")$Breed),
                        decreasing = TRUE),16)
survivalData.top16 <- subset(survivalData,Breed_Class=="AKC-Recognized Breed" & 
                               Breed %in% names(top16breeds))
survivalData.top16$Breed <- factor(survivalData.top16$Breed,
                                   levels=names(top16breeds))
print(top16breeds)

## Summary Statistics for 16 Individual Breeds
surv.top16 <- Surv(time=survivalData.top16$first.age,
                   time2=survivalData.top16$last.age,
                   event=survivalData.top16$event, 
                   type="counting")
fit.top16.s <- survfit(surv.top16 ~ Breed+Sex,data=survivalData.top16)
fit.top16.s.sum <- summary(fit.top16.s)
fit.top16.s.sum.df <- as.data.frame(fit.top16.s.sum$table)
fit.top16.s.quant <- quantile(fit.top16.s)
fit.top16.s.sum.df$q25 <- fit.top16.s.quant$quantile[,1]
fit.top16.s.sum.df$q25.lcl <- fit.top16.s.quant$lower[,1]
fit.top16.s.sum.df$q25.ucl <- fit.top16.s.quant$upper[,1]
fit.top16.s.sum.df$q75 <- fit.top16.s.quant$quantile[,3]
fit.top16.s.sum.df$q75.lcl <- fit.top16.s.quant$lower[,3]
fit.top16.s.sum.df$q75.ucl <- fit.top16.s.quant$upper[,3]
fit.top16.s.sum.df$stratum <- rownames(fit.top16.s.sum$table)
strata<-trimws(unlist(strsplit(unlist(
  strsplit(rownames(fit.top16.s.sum$table),",")),"=")))
strata <- strata[seq(2,length(strata),2)]
strata <- matrix(strata,ncol=2,byrow=TRUE)
fit.top16.s.sum.df$Breed <-
  factor(strata[,1],levels=levels(survivalData.top16$Breed))
fit.top16.s.sum.df$Sex <- 
  factor(strata[,2],levels=levels(survivalData.top16$Sex))
write.csv(fit.top16.s.sum.df,file.path(resultsfolder,
                                       "Survival_Top16_demographics.SumStats.csv"))

plt.forest.median.top16.s<-ggplot(fit.top16.s.sum.df)+
  geom_errorbarh(aes(xmin=`0.95LCL`,xmax=`0.95UCL`,
                     y=Sex),height=0.1)+
  geom_errorbarh(aes(xmin=`q25.lcl`,xmax=`q25.ucl`,
                     y=Sex),height=0,alpha=0.5)+
  geom_errorbarh(aes(xmin=`q75.lcl`,xmax=`q75.ucl`,
                     y=Sex),height=0,alpha=0.5)+
  geom_point(aes(x=median,y=Sex,shape="median"))+
  geom_point(aes(x=q25,y=Sex,shape="IQR"))+
  geom_point(aes(x=q75,y=Sex,shape="IQR"))+
  xlab("Median & IQR lifespan [CI]")+
  scale_shape_discrete("",limits=rev)+
  theme_bw()+theme(legend.position = "bottom",
                   strip.text.x = element_text(margin = margin(1,0,1,0)))+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  facet_wrap(~Breed,ncol=1)
print(plt.forest.median.top16.s)

plt.forest.median.top16.s<-ggplot(fit.top16.s.sum.df,aes(y=Sex))+
  geom_boxplot(aes(xmin=q25,xlower=q25,xmiddle=median,xupper=q75,xmax=q75),
               fill=NA,width=0.5,stat="identity")+
  xlab("Median & IQR lifespan")+
  theme_bw()+theme(legend.position="bottom")+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  theme(legend.position = "bottom",
        strip.text.x = element_text(margin = margin(1,0,1,0)))+
  scale_x_continuous(minor_breaks=seq(5,20))+
  facet_wrap(~Breed,ncol=1)
print(plt.forest.median.top16.s)

plt.forest.mean.top16.s<-ggplot(fit.top16.s.sum.df)+
  geom_errorbarh(aes(xmin=rmean-1.96*`se(rmean)`,
                     xmax=rmean+1.96*`se(rmean)`,
                     y=Sex),height=0)+
  geom_point(aes(x=rmean,y=Sex,shape="mean"))+
  xlab("Mean lifespan [CI]")+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  theme_bw()+theme(legend.position = "none",
                   strip.text.x = element_text(margin = margin(1,0,1,0)))+
  facet_wrap(~Breed,ncol=1)
print(plt.forest.mean.top16.s)

ggsave(file.path(figfolder,"Supp.Fig.Survival_Top16_demographics.forest.pdf"),
       ggarrange(plt.forest.median.top16.s,plt.forest.mean.top16.s,nrow=1),
       height=7,width=7,scale=1.2)

## Test Sex Effect in Each 16 Individual Breeds
cox.top16.s <- coxph(surv.top16 ~ Breed+Sex,data=survivalData.top16)
cox.top16.s.pvalues <- data.frame()
for (Breed_now in levels(survivalData.top16$Breed)) {
  cox.top16.tmp <- coxph(surv.top16 ~ Sex,data=survivalData.top16,
                         subset=Breed == Breed_now)
  cox.top16.s.pvalues <- rbind(cox.top16.s.pvalues,
                               data.frame(Breed=Breed_now,
                                          log.rank.pvalue=summary(cox.top16.tmp)$sctest["pvalue"]))
}
cox.top16.s.pvalues$Breed <- factor(cox.top16.s.pvalues$Breed,
                                        levels=levels(survivalData.top16$Breed))
write.csv(cox.top16.s.pvalues,
          file.path(resultsfolder,"Survival_Top16_demographics_sexeffect-pvals.csv"))
plt.km.top16.s<-ggsurvplot_facet(fit.top16.s,data=survivalData.top16,
                                 facet.by=c("Breed"),
                                 censor=FALSE,xlab="Age (yr)",
                                 surv.median.line = "hv",short.panel.labs=TRUE)+
  scale_color_viridis_d(option="magma",end=0.8)+
  geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.top16.s.pvalues,hjust=0)
print(plt.km.top16.s)
ggsave(file.path(figfolder,"Supp.Fig.Survival_Top16_demographics.sex_effect.KM.pdf"),
       plt.km.top16.s,
       height=4,width=7,scale=1.6)

## Summary Statistics for 16 Individual Breeds sexes combined
surv.top16 <- Surv(time=survivalData.top16$first.age,
                   time2=survivalData.top16$last.age,
                   event=survivalData.top16$event, 
                   type="counting")
fit.top16 <- survfit(surv.top16 ~ Breed,data=survivalData.top16)
fit.top16.sum <- summary(fit.top16)
fit.top16.sum.df <- as.data.frame(fit.top16.sum$table)
fit.top16.quant <- quantile(fit.top16)
fit.top16.sum.df$q25 <- fit.top16.quant$quantile[,1]
fit.top16.sum.df$q25.lcl <- fit.top16.quant$lower[,1]
fit.top16.sum.df$q25.ucl <- fit.top16.quant$upper[,1]
fit.top16.sum.df$q75 <- fit.top16.quant$quantile[,3]
fit.top16.sum.df$q75.lcl <- fit.top16.quant$lower[,3]
fit.top16.sum.df$q75.ucl <- fit.top16.quant$upper[,3]
fit.top16.sum.df$stratum <- rownames(fit.top16.sum$table)
strata<-trimws(unlist(strsplit(unlist(
  strsplit(rownames(fit.top16.sum$table),",")),"=")))
strata <- strata[seq(2,length(strata),2)]
strata <- matrix(strata,ncol=1,byrow=TRUE)
fit.top16.sum.df$Breed <-
  factor(strata[,1],levels=levels(survivalData.top16$Breed))
write.csv(fit.top16.sum.df,file.path(resultsfolder,
                                       "Survival_Top16_demographics_bothsexes.SumStats.csv"))


## Test if Survival of Each 16 Individual Breeds is Consistent with its Assigned Size Class
top16breeds.labels <- paste0(names(top16breeds)," (n=",top16breeds,")")
names(top16breeds.labels) <- names(top16breeds)
cox.top16.strata.pvalues <- data.frame()
top16.strata.plt <- list()
for (j in 1:length(top16breeds)) {
  Breed.now <- names(top16breeds)[j]
  Size_Class.now <- first(subset(survivalData.top16,Breed==Breed.now)$Size_Class_at_HLES)
  survivalData.tmp <- subset(survivalData,Size_Class_at_HLES==Size_Class.now)
  survivalData.tmp$CompBreed <- Breed.now
  survivalData.tmp$isBreed <- ifelse(survivalData.tmp$Breed == Breed.now,
                                     Breed.now,Size_Class.now)
  surv.tmp <- Surv(time=survivalData.tmp$first.age,
                   time2=survivalData.tmp$last.age,
                   event=survivalData.tmp$event,
                   type='counting')
  fit.tmp <- survfit(surv.tmp ~ isBreed+Sex,data=survivalData.tmp)
  for (Sex.now in levels(survivalData.tmp$Sex)) {
    survivalData.tmp.sex <- subset(survivalData.tmp,Sex==Sex.now)
    surv.tmp.sex <- Surv(time=survivalData.tmp.sex$first.age,
                     time2=survivalData.tmp.sex$last.age,
                     event=survivalData.tmp.sex$event,
                     type='counting')
    cox.tmp <- coxph(surv.tmp.sex ~ isBreed,data=survivalData.tmp.sex)
    cox.top16.strata.pvalues <- rbind(cox.top16.strata.pvalues,
                                      data.frame(Breed=Breed.now,
                                                 Size_Class_at_HLES=Size_Class.now,
                                                 Sex=Sex.now,
                                                 log.rank.pvalue=summary(cox.tmp)$sctest["pvalue"]
                                                 ))
  }
  ggplt<-ggsurvplot_facet(fit.tmp,data=survivalData.tmp,
                          facet.by=c("Sex","CompBreed"),
                    censor=FALSE,xlab="Age (yr)",
                    short.panel.labs=TRUE,surv.median.line = "hv")+
                    scale_color_viridis_d(option="magma",end=0.8,
                        labels=c(paste0("Other ",Size_Class.now," (n=",
                                        nrow(survivalData.tmp)-
                                          top16breeds[Breed.now],")"),
                                 paste(top16breeds.labels[Breed.now])))+
        geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=rename(subset(cox.top16.strata.pvalues,Breed==Breed.now),
                        CompBreed=Breed),hjust=0)+
    scale_x_continuous(limits=c(0,25),breaks=seq(0,25,5))+
    coord_cartesian(xlim=c(0,25))+
    theme(legend.position = "bottom")+
        guides(color = guide_legend(reverse = TRUE,nrow=2))
  top16.strata.plt[[j]] <- ggplt
}
write.csv(cox.top16.strata.pvalues,
          file.path(resultsfolder,"Survival_Top16_demographics.size_class-pvals.csv"))
ggsave(file.path(figfolder,"Supp.Fig.Survival_Top16_demographics.size_class.KM.pdf"),
       ggarrange(plotlist=top16.strata.plt,ncol=4,nrow=4),
       height=8,width=6,scale=2.75)
