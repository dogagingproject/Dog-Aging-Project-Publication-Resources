# SWAS figures — sensitivity: restricted to Mature Adult lifestage.
library(ggplot2)
library(tidyverse)
library(stringr)
library(cowplot)
figfolder <- "Figures"
resultsfolder <- "Results"
load(file.path(resultsfolder,"Supp.HLES_Cox.Rdata"))
load(file.path(resultsfolder,"Supp.HLES_Cox.MAdult.Rdata"))
pvals.comp<-full_join(pvals,pvals.MAdult,by=names(pvals)[c(1:9,11)])
# Spearman cor = 0.68
cor(pvals.comp$score.pval.adjust.x,pvals.comp$score.pval.adjust.y,
    method="spearman",use="complete.obs")
plot(gplots::venn(list(MatureAdultDogs=vars.to.plot.MAdult$Variable,AllDogs=vars.to.plot$Variable)))
mtext("q < 0.05",side=3)

plot(gplots::venn(list(MatureAdultDogs=subset(vars.to.plot.MAdult,score.pval.adjust<=3e-4)$Variable,
                       AllDogs=subset(vars.to.plot,score.pval.adjust<=3e-4)$Variable)))
mtext("q < 0.0003",side=3)

vars.both <-
  vars.to.plot.MAdult$Variable[(vars.to.plot.MAdult$Variable %in% vars.to.plot$Variable)]
print(vars.both)

# [1] hs_general_health                                         
# [2] hs_chronic_condition_present                              
# [3] pa_activity_level                                         
# [4] pa_on_leash_walk_slow_pace_pct                            
# [5] pa_physical_games_frequency                               
# [6] db_energetic_frequency                                    
# [7] hs_health_conditions_endocrine                            
# [8] pa_avg_activity_intensity                                 
# [9] hs_health_conditions_cancer                               
# [10] pa_swim                                                   
# [11] pa_off_leash_walk_slow_pace_pct                           
# [12] hs_new_condition_diagnosed_recently                       
# [13] df_primary_diet_component_change_health_condition_specific
# [14] db_chases_squirrels_frequency                             
# [15] de_recreational_spaces                                    
# [16] hs_health_conditions_cardiac                              
# [17] df_appetite_change_last_year                              
# [18] hs_chronic_condition_recently_changed_or_treated          
# [19] pa_on_leash_walk_brisk_pace_pct                           
# [20] db_training_obeys_sit_command_frequency                   
# [21] db_training_obeys_stay_command_frequency                  
# [22] hs_health_conditions_neurological                         
# [23] df_ever_underweight                                       
# [24] db_playful_frequency                                      
# [25] de_routine_toys                                           
# [26] de_daytime_sleep_avg_hours                                
# [27] db_pulls_leash_frequency                                  
# [28] df_ever_overweight                                        
# [29] db_chases_birds_frequency                                 
# [30] mp_recent_non_prescription_meds_upset_stomach             
# [31] pa_other_aerobic_activity_avg_intensity                   
# [32] de_eats_grass_frequency                                   
# [33] df_weight_change_last_year                                
# [34] de_drinks_outdoor_water                                   
# [35] hs_recent_hospitalization                                 
# [36] de_stairs_in_home                                         
# [37] de_home_years_lived_in                                    
# [38] pa_off_leash_walk_brisk_pace_pct  
## Manhattan plot

ymax <- max(-log10(pvals.MAdult$score.pval.adjust),na.rm=TRUE)
plt.man.hles<-
  ggplot(subset(pvals.MAdult,!is.na(score.pval.adjust)))+
  geom_point(aes(y=-log10(score.pval.adjust),x=as.numeric(Variable),
                 color=variable_group_name,size=score.pval.adjust<0.05))+
  geom_segment(aes(y=-log10(score.pval.adjust),x=as.numeric(Variable),
                   color=variable_group_name,yend=0,xend=as.numeric(Variable)))+
  geom_rect(aes(xmin=xmin,xmax=xmax,ymin=-10,ymax=-1,fill=variable_group_name),
            data=variable_group_tab.df)+
  scale_fill_discrete("",labels=gsub(" ","\n",
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
                              legend.text = element_text(size = 9))+
  scale_y_continuous(breaks=seq(0,200,50),limits=c(-75,ymax))+
  geom_hline(yintercept=-log10(c(0.05,1e-6)),linetype=c("dotted","dashed"),color="grey50")+
  guides(color="none",size="none",fill=guide_legend(nrow=1))+
  ggtitle("Restricted to Dogs in Mature Adult Lifestage at Baseline")
print(plt.man.hles)

plt.man.hles.zoom<-
  ggplot(subset(pvals.MAdult,!is.na(score.pval.adjust)))+
  geom_point(aes(y=-log10(score.pval.adjust),x=as.numeric(Variable),
                 color=variable_group_name,size=score.pval.adjust<0.05))+
  geom_segment(aes(y=-log10(score.pval.adjust),x=as.numeric(Variable),
                   color=variable_group_name,yend=0,xend=as.numeric(Variable)))+
  geom_rect(aes(xmin=xmin,xmax=xmax,ymin=-1,ymax=0,fill=variable_group_name),
            data=variable_group_tab.df)+
  scale_fill_discrete("",labels=gsub(" ","\n",
                                     as.character(variable_group_tab.df$variable_group_name)))+
  scale_size_manual("p.adjust<0.05",values=c(0.01,1))+
  ylab("")+
  theme_bw()+
  xlab("HLES Variable")+theme(axis.text.x = element_blank(),
                              axis.ticks = element_blank(),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              legend.position = "inside",
                              legend.position.inside = c(0.5,0.01),
                              legend.justification = c("center","bottom"),
                              legend.text = element_text(size = 9))+
  coord_cartesian(ylim=c(-1,12))+scale_y_continuous(breaks=seq(0,12,2))+
  geom_hline(yintercept=-log10(c(0.05,1e-6)),linetype=c("dotted","dashed"),color="grey50")+
  guides(color="none",size="none",fill="none")
print(plt.man.hles.zoom)

suppfig.with.inset<-ggdraw() + draw_plot(plt.man.hles) +
  draw_plot(plt.man.hles.zoom,x=0.05,y=0.6,width=0.8,height=0.3)

ggsave(file.path(figfolder,"Supp.Fig.Manh.MAdult.pdf"),suppfig.with.inset,height=2,width=6,scale=2)

## Plots of categorical variables
plt.cox.categorical.list <- list()
pdf(file.path(figfolder,"HLES_Cox_Details_SensAn","HLES_HR.cox.categorical.MAdult.pdf"),height=4.5,width=6)
for (j in 1:nrow(vars.to.plot.MAdult)) {
  var.now <- as.character(vars.to.plot.MAdult$Variable[j])
  pval.now <- vars.to.plot.MAdult$score.pval.adjust[j]
  cox.coef.tmp <- subset(cox.coef.MAdult,Variable == var.now)
  if (cox.coef.tmp$Type[1] == "factor") {
    print(paste(j,var.now,"----------"))
    labels.now <- trimws(strsplit(vars.to.plot.MAdult$ValueLabels[j],"\\|")[[1]])
    names(labels.now) <- trimws(strsplit(vars.to.plot.MAdult$Values[j],"\\|")[[1]])
    if (sum(duplicated(labels.now))>0) {
      labels.now <- paste(labels.now,names(labels.now),sep="-")
    }
    cox.coef.tmp$y <- factor(sub(".","",cox.coef.tmp$y),
                             levels=labels.now)
    cox.coef.tmp <- cox.coef.tmp[order(cox.coef.tmp$y,decreasing = TRUE),]
    cox.coef.tmp[nrow(cox.coef.tmp)+1,] <- cox.coef.tmp[1,]
    n <- nrow(cox.coef.tmp)
    rownames(cox.coef.tmp)[n] <- paste0("x",cox.coef.tmp$Ref[1]) # labels.now[1])
    cox.coef.tmp$y[n] <- cox.coef.tmp$Ref[1] # labels.now[1]
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
    ggsave(file.path(figfolder,"HLES_Cox_Details_SensAn",paste0("HLES_HR_cox.",var.now,".MAdult.pdf")),
           plt.cox.categorical.tmp,height=4.5,width=6)
  }
}
dev.off()

## Plots of numeric variables
cox.coef.MAdult.num <- subset(cox.coef.MAdult,(Type == "numeric"))
cox.coef.MAdult.num <- cox.coef.MAdult.num[order(cox.coef.MAdult.num$coef),]
cox.coef.MAdult.num <- left_join(cox.coef.MAdult.num,vars.to.plot.MAdult)
cox.coef.MAdult.num$Variable <- factor(cox.coef.MAdult.num$Variable,
                                      levels=cox.coef.MAdult.num$Variable)
plt.cox.numeric <-
  ggplot(cox.coef.MAdult.num)+
  geom_vline(xintercept=1,linetype="dotted")+
  geom_point(aes(x=`exp(coef)`,y=Variable,color=variable_group_name),shape=15,size=3)+
  geom_errorbarh(aes(xmin=`lower .95`,xmax=`upper .95`,y=Variable),height=0)+
  geom_text(aes(x=1/4,y=Variable,
                label=paste0(round(`exp(coef)`,3),
                             " [",round(`lower .95`,3),"-",
                             round(`upper .95`,3),"]")),
            hjust=0,vjust=0.5,size=2.5,data=cox.coef.MAdult.num)+
  geom_text(aes(x=3,y=Variable,
                label=paste0("p(FDR-adj)=")),
            hjust=0.5,vjust=-1,size=2.5,data=last(cox.coef.MAdult.num))+
  geom_text(aes(x=3,y=Variable,
                label=paste0(" ",signif(score.pval.adjust,3))),
            hjust=0.5,vjust=0.5,size=2.5,data=cox.coef.MAdult.num)+
  xlab("Hazard Ratio per unit increase variable")+ylab("")+
  scale_color_discrete("")+
  theme_bw()+annotation_logticks(sides="b")+theme(legend.position = "bottom")+
  scale_x_log10()+coord_cartesian(xlim=c(1/4,4))
print(plt.cox.numeric)
ggsave(file.path(figfolder,"HLES_Cox_Details_SensAn","HLES_HR.cox.numeric.MAdult.pdf"),plt.cox.numeric,height=4,width=6)
