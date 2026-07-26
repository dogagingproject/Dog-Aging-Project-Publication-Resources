library(ggplot2)
library(tidyverse)
library(stringr)
library(cowplot)
figfolder <- "Figures"
resultsfolder <- "Results"
load(file.path(resultsfolder,"Supp.HLES_Cox.Rdata"))


## Table 4
cox.coef$abs.coef <- abs(cox.coef$coef)
cox.coef.maxeffect <- aggregate(abs.coef ~ Variable,data=cox.coef,FUN=max)
cox.coef.maxeffect <- left_join(cox.coef.maxeffect,cox.coef)
cox.coef.maxeffect <- cox.coef.maxeffect[order(
  cox.coef.maxeffect$score.pval.adjust),]


write.csv(cox.coef.maxeffect,
          file.path(resultsfolder,"Supp.HLES_Cox_maxeffects.csv"))

## Manhattan plot

ymax <- max(-log10(pvals$score.pval.adjust),na.rm=TRUE)
pvals.tmp<-subset(pvals,!is.na(score.pval.adjust))
indx.min <- which.min(pvals.tmp$score.pval.adjust)
pvals.tmp$score.pval.adjust[indx.min] <- 0
plt.man.hles<-
  ggplot(pvals.tmp)+
  geom_point(aes(y=log10(1/score.pval.adjust),x=as.numeric(Variable),
                 color=variable_group_name,size=score.pval.adjust<0.05))+
  geom_segment(aes(y=log10(1/score.pval.adjust),x=as.numeric(Variable),
                   color=variable_group_name,yend=0,xend=as.numeric(Variable)))+
  geom_label(aes(y=57,x=as.numeric(Variable),color=variable_group_name,label="//"),
             label.size=0,label.padding = unit(0, "lines"),angle=90,
             data=subset(pvals.tmp,score.pval.adjust==0))+
  geom_rect(aes(xmin=xmin,xmax=xmax,ymin=-3,ymax=-1,fill=variable_group_name),
            data=variable_group_tab.df)+
  scale_fill_viridis_d("",option="turbo",labels=gsub(" ","\n",
                                     as.character(variable_group_tab.df$variable_group_name)))+
  scale_size_manual("p.adjust<0.05",values=c(0.01,1))+
  ylab("-Log10(FDR-adjusted p value)")+
  theme_bw()+
  xlab("HLES Variable")+theme(axis.text.x = element_blank(),
                              axis.ticks = element_blank(),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              legend.position = "inside",
                              legend.position.inside = c(0.5,0.01),
                              legend.justification = c("center","bottom"),
                              legend.text = element_text(size = 9),
                              panel.border=element_blank())+
  scale_y_continuous(breaks=seq(0,50,10),limits=c(-75,ymax))+
  coord_cartesian(ylim=c(-25,59),xlim=c(-1,1)+range(as.numeric(pvals.tmp$Variable)),
                  clip="off",expand=FALSE)+
  scale_color_viridis_d(option="turbo")+
  geom_hline(yintercept=-log10(c(0.05,1e-6)),linetype=c("dotted","dashed"),color="grey50")+
  guides(color="none",size="none",fill=guide_legend(nrow=1))
print(plt.man.hles)

ggsave(file.path(figfolder,"Fig.3.Manh.pdf"),plt.man.hles,height=2,width=5.75,scale=2)

## Plots of categorical variables
plt.cox.categorical.list <- list()
pdf(file.path(figfolder,"HLES_Cox_Details","HLES_HR.cox.categorical.pdf"),height=4.5,width=6)
for (j in 1:nrow(vars.to.plot)) {
  var.now <- as.character(vars.to.plot$Variable[j])
  pval.now <- vars.to.plot$score.pval.adjust[j]
  cox.coef.tmp <- subset(cox.coef,Variable == var.now)
  if (cox.coef.tmp$Type[1] == "factor") {
    print(paste(j,var.now,"----------"))
    labels.now <- trimws(strsplit(vars.to.plot$ValueLabels[j],"\\|")[[1]])
    names(labels.now) <- trimws(strsplit(vars.to.plot$Values[j],"\\|")[[1]])
    if (sum(duplicated(labels.now))>0) {
      labels.now <- paste(labels.now,names(labels.now),sep="-")
    }
    cox.coef.tmp$y <- factor(sub(".","",cox.coef.tmp$y),
                             levels=labels.now)
    cox.coef.tmp <- cox.coef.tmp[order(cox.coef.tmp$y,decreasing = TRUE),]
    cox.coef.tmp[nrow(cox.coef.tmp)+1,] <- cox.coef.tmp[1,]
    n <- nrow(cox.coef.tmp)
    rownames(cox.coef.tmp)[n] <- paste0("x",cox.coef.tmp$Ref[1]) #labels.now[1])
    cox.coef.tmp$y[n] <-  cox.coef.tmp$Ref[1] # labels.now[1]
    cox.coef.tmp$Variable.y[n] <- paste(cox.coef.tmp$Variable[n],cox.coef.tmp$y[n],sep="|")
    cox.coef.tmp[n,c("coef","exp(coef)","se(coef)","z","Pr(>|z|)",
                     "exp(-coef)","lower .95","upper .95")]<-
      c(0,1,0,0,1,1,1,1)
    var.lab <- paste0(cox.coef.tmp$Variable,"\n",str_wrap(cox.coef.tmp$SurveyText,width=70))
    names(var.lab) <- cox.coef.tmp$Variable
    plt.cox.categorical.tmp <-
      ggplot(cox.coef.tmp)+
      geom_vline(xintercept=1,linetype="dotted")+
      geom_point(aes(x=`exp(coef)`,y=y,shape=(y==Ref)),size=3,color="grey50")+
      geom_errorbarh(aes(xmin=`lower .95`,xmax=`upper .95`,y=y),height=0)+
      geom_text(aes(x=1/10,y=y,
                    label=paste0(round(`exp(coef)`,2),
                                 " [",round(`lower .95`,2),"-",
                                 round(`upper .95`,2),"]")),
                hjust=0,vjust=0.5,size=3,
                data=cox.coef.tmp[-nrow(cox.coef.tmp),])+
      geom_text(aes(x=1/10,y=y,
                    label="Reference"),
                hjust=0,vjust=0.5,size=3,
                data=last(cox.coef.tmp))+
      geom_text(aes(x=10,y=y),
                label=" p(FDR-adj) =",
                hjust=1,vjust=-1,size=3,
                data=last(cox.coef.tmp))+
      geom_text(aes(x=10,y=y),
                label=paste(signif(pval.now,3)),
                hjust=1,vjust=0.5,size=3,
                data=last(cox.coef.tmp))+
      xlab("Hazard Ratio")+ylab("")+
      theme_bw()+annotation_logticks(sides="b")+
      scale_shape_manual(values=c(15,1))+
      facet_wrap(~Variable,labeller = labeller(Variable = var.lab))+guides(shape="none")+
      scale_y_discrete(labels = function(x) str_wrap(x, width = 25),limits=rev)+
      scale_x_log10()+coord_cartesian(xlim=c(0.1,10))
    print(plt.cox.categorical.tmp)
    plt.cox.categorical.list[[j]] <- plt.cox.categorical.tmp
    ggsave(file.path(figfolder,"HLES_Cox_Details",paste0("HLES_HR_cox.",var.now,".pdf")),
           plt.cox.categorical.tmp,height=4.5,width=6)
  }
}
dev.off()

## Plots of numeric variables
cox.coef.num <- subset(cox.coef,(Type == "numeric"))
cox.coef.num <- cox.coef.num[order(cox.coef.num$coef),]
cox.coef.num <- left_join(cox.coef.num,vars.to.plot)
cox.coef.num$Variable <- factor(cox.coef.num$Variable,
                                levels=cox.coef.num$Variable)
plt.cox.numeric <-
  ggplot(cox.coef.num)+
  geom_vline(xintercept=1,linetype="dotted")+
  geom_point(aes(x=`exp(coef)`,y=Variable,color=variable_group_name),shape=15,size=3)+
  geom_errorbarh(aes(xmin=`lower .95`,xmax=`upper .95`,y=Variable),height=0)+
  geom_text(aes(x=1/4,y=Variable,
                label=paste0(round(`exp(coef)`,3),
                             " [",round(`lower .95`,3),"-",
                             round(`upper .95`,3),"]")),
            hjust=0,vjust=0.5,size=2.5,data=cox.coef.num)+
  geom_text(aes(x=3,y=Variable,
                label=paste0("p(FDR-adj)=")),
            hjust=0.5,vjust=-1,size=2.5,data=last(cox.coef.num))+
  geom_text(aes(x=3,y=Variable,
                label=paste0(" ",signif(score.pval.adjust,3))),
            hjust=0.5,vjust=0.5,size=2.5,data=cox.coef.num)+
  xlab("Hazard Ratio per unit increase variable")+ylab("")+
  scale_color_discrete("")+
  theme_bw()+annotation_logticks(sides="b")+theme(legend.position = "bottom")+
  scale_x_log10()+coord_cartesian(xlim=c(1/4,4))
print(plt.cox.numeric)
ggsave(file.path(figfolder,"HLES_Cox_Details","HLES_HR.cox.numeric.pdf"),plt.cox.numeric,height=4,width=6)

cox.coef.num.brisk.slow <-cox.coef.num[grepl("brisk",cox.coef.num$Variable) |
                                         grepl("slow",cox.coef.num$Variable),]
cox.coef.num.brisk.slow$SurveyText <- str_wrap(cox.coef.num.brisk.slow$SurveyText,
                                               width=25)
cox.coef.num.brisk.slow <- cox.coef.num.brisk.slow[order(cox.coef.num.brisk.slow$SurveyText),]
plt.cox.numeric.brisk.slow <-
  ggplot(cox.coef.num.brisk.slow)+
  geom_vline(xintercept=1,linetype="dotted")+
  geom_point(aes(x=`exp(coef)`,y=SurveyText),shape=15,size=3,color="grey50")+
  geom_errorbarh(aes(xmin=`lower .95`,xmax=`upper .95`,y=SurveyText),height=0)+
  geom_text(aes(x=1/4,y=SurveyText,
                label=paste0(round(`exp(coef)`,2),
                             " [",round(`lower .95`,2),"-",
                             round(`upper .95`,2),"]")),
            hjust=0,vjust=0.5,size=3,data=cox.coef.num.brisk.slow)+
  geom_text(aes(x=3,y=SurveyText,
                label=paste0("p(FDR-adj)=")),
            hjust=0.5,vjust=-1,size=3,data=last(cox.coef.num.brisk.slow))+
  geom_text(aes(x=3,y=SurveyText,
                label=paste0(" ",signif(score.pval.adjust,3))),
            hjust=0.5,vjust=0.5,size=3,data=cox.coef.num.brisk.slow)+
  xlab("Hazard Ratio for change from 0%-100%")+ylab("")+
  theme_bw()+annotation_logticks(sides="b")+theme(legend.position = "bottom")+
  scale_x_log10()+coord_cartesian(xlim=c(1/4,4))+
  facet_wrap(~variable_group_name)
print(plt.cox.numeric.brisk.slow)
ggsave(file.path(figfolder,"HLES_Cox_Details","HLES_HR.cox.numeric.brisk.slow.pdf"),
       plt.cox.numeric.brisk.slow,height=4.5,width=6)

## Make Figure 3

ggsave(file.path(figfolder,"Fig.4.HLES_Cox_Exp_Resp.pdf"),
       ggarrange(plt.cox.categorical.list[[1]],
                 plt.cox.categorical.list[[5]],
                 plt.cox.categorical.list[[21]],
                 plt.cox.numeric.brisk.slow,ncol=2,nrow=2,labels="auto"),
       scale=1.5,height=4,width=6)
