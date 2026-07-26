library(haven)
library(lubridate)
library(survival)
library(flexsurv)
library(survminer)
library(dplyr)
library(viridisLite)
library(usdata)
library(wCorr)
library(emmeans)
source("ggforest2.R")
datfolder <- "../data"
resultsfolder <- "results"
figfolder <- "figures"

## Join geographic variables
load(file.path(datfolder,"SurvivalData.RData"))
load(file.path(datfolder,"DAP_2024_HLES_dog_owner_v1.0.RData"))
survivalData <- left_join(survivalData,
                          HLES_dog_owner[,c("dog_id",
                                            "de_home_area_type",
                                            "oc_primary_residence_census_division",
                                            "oc_primary_residence_state")])
remove(HLES_dog_owner)
level.list <- list()
# Largest N is reference for each one
level.list[[1]] <- c("Suburban","Rural","Urban") 
level.list[[2]] <- rev(c("New England","Middle Atlantic",
                         "East North Central","West North Central",
                         "South Atlantic","East South Central","West South Central",
                         "Mountain","Pacific")) 
level.list[[3]] <- sort(unique(survivalData$oc_primary_residence_state))
ref.state <- which(level.list[[3]]=="WA")
level.list[[3]] <- c(level.list[[3]][ref.state],level.list[[3]][-ref.state])

survivalData$de_home_area_type <- factor(
  c("Urban","Suburban","Rural")[survivalData$de_home_area_type],
  levels=level.list[[1]])

survivalData$oc_primary_residence_census_division <- factor( 
  c("New England","Middle Atlantic",
    "East North Central","West North Central",
    "South Atlantic","East South Central","West South Central",
    "Mountain","Pacific")[
      survivalData$oc_primary_residence_census_division],
  levels=level.list[[2]])
survivalData$oc_primary_residence_state <- factor(
  survivalData$oc_primary_residence_state,
  levels=level.list[[3]])

survivalData <- rename(survivalData,
                       home_area = de_home_area_type,
                       census_division = oc_primary_residence_census_division,
                       state = oc_primary_residence_state)
  
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')



###### Unadjusted

fit.geo.unadj.ha <- survfit(surv ~ home_area,data=survivalData)
cox.geo.unadj.ha <- coxph(surv ~ home_area,data=survivalData)
fit.geo.unadj.cd <- survfit(surv ~ census_division,data=survivalData)
cox.geo.unadj.cd <- coxph(surv ~ census_division,data=survivalData)
fit.geo.unadj.st <- survfit(surv ~ state,data=survivalData)
cox.geo.unadj.st <- coxph(surv ~ state,data=survivalData)

cox.geo.unadj.pvalues <- data.frame(Analysis="Unadjusted",
                                    Variable=c("home_area","census_division","state"),
                                    log.rank.pvalue=c(summary(cox.geo.unadj.ha)$sctest["pvalue"],
                                                      summary(cox.geo.unadj.cd)$sctest["pvalue"],
                                                      summary(cox.geo.unadj.st)$sctest["pvalue"]),
                                    cox.zph.pvalue=c(
                                      cox.zph(cox.geo.unadj.ha)$table[1,"p"],
                                      cox.zph(cox.geo.unadj.cd)$table[1,"p"],
                                      cox.zph(cox.geo.unadj.st)$table[1,"p"])
)
print(cox.geo.unadj.pvalues)

###### Distribution of size, breed, sex

## Home area 

plt.size.by.area.ha <- ggplot(survivalData, 
                           aes(x = home_area, fill = Size_Class_at_HLES)) +
  scale_fill_viridis_d()+ggtitle("Demographics by Home Area Type")+
  scale_x_discrete(limits=c("Rural","Suburban","Urban"))+
  geom_bar(position="fill")+ylab("proportion")
plt.breed.by.area.ha <- ggplot(survivalData, 
                            aes(x = home_area, fill = Breed_Class)) +
  scale_fill_viridis_d(end=0.8)+
  scale_x_discrete(limits=c("Rural","Suburban","Urban"))+
  geom_bar(position="fill")+ylab("proportion")
plt.sex.by.area.ha <- ggplot(survivalData, 
                          aes(x = home_area, fill = Sex)) +
  scale_fill_viridis_d(end=0.8)+
  scale_x_discrete(limits=c("Rural","Suburban","Urban"))+
  geom_bar(position="fill")+ylab("proportion")

ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_DemographicDist_homearea.pdf"),
       ggarrange(plt.size.by.area.ha,
                plt.breed.by.area.ha,
                plt.sex.by.area.ha,ncol=1),
       height=7,width=5,scale=1.2)

## Census Division
plt.size.by.area.cd <- ggplot(survivalData, 
                           aes(x = as.factor(as.integer(census_division)), fill = Size_Class_at_HLES)) +
  geom_bar(position="fill")+ggtitle("Demographics by Census Division")+
  scale_fill_viridis_d()+ylab("proportion")+xlab("Census Division #")
plt.breed.by.area.cd <- ggplot(survivalData, 
                            aes(x = as.factor(as.integer(census_division)), fill = Breed_Class)) +
  geom_bar(position="fill")+scale_fill_viridis_d(end=0.8)+ylab("proportion")+xlab("Census Division #")
plt.sex.by.area.cd <- ggplot(survivalData, 
                          aes(x = as.factor(as.integer(census_division)), fill = Sex)) +
  geom_bar(position="fill")+scale_fill_viridis_d(end=0.8)+ylab("proportion")+xlab("Census Division #")

ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_DemographicDist_censusdivision.pdf"),
       ggarrange(plt.size.by.area.cd,
                 plt.breed.by.area.cd,
                 plt.sex.by.area.cd,ncol=1),
       height=7,width=5,scale=1.2)

## State

plt.size.by.area.st <- ggplot(survivalData, 
                           aes(y = state, fill = Size_Class_at_HLES)) +
  geom_bar(position="fill")+ggtitle("Demographics by State")+
  xlab("proportion")+ylab("state")+
  scale_y_discrete(limits=rev(sort(unique(as.character(survivalData$state)))))+
  scale_fill_viridis_d(guide = guide_legend(reverse = TRUE,position="bottom",nrow=2))


plt.breed.by.area.st <- ggplot(survivalData, 
                            aes(y = state, fill = Breed_Class)) +
  geom_bar(position="fill")+xlab("proportion")+ylab("state")+ggtitle("")+
  scale_y_discrete(limits=rev(sort(unique(as.character(survivalData$state)))))+
  scale_fill_viridis_d(end=0.8,guide = guide_legend(reverse = TRUE,position="bottom",nrow=2))

plt.sex.by.area.st <- ggplot(survivalData, 
                          aes(y = state, fill = Sex)) +
  geom_bar(position="fill")+xlab("proportion")+ylab("state")+ggtitle("")+
  scale_y_discrete(limits=rev(sort(unique(as.character(survivalData$state)))))+
  scale_fill_viridis_d(end=0.8,guide = guide_legend(reverse = TRUE,position="bottom",nrow=2))


ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_DemographicDist_state.pdf"),
       ggarrange(plt.size.by.area.st,
                 plt.breed.by.area.st,
                 plt.sex.by.area.st,ncol=3,
                 widths=c(2,1,1)),
       height=7,width=8.5,scale=1.5)

###### Adjusted Analyses

## Common geo effect


fit.geo.adj.ha <- survfit(surv ~ home_area + 
                                 Size_Class_at_HLES + Breed_Class + Sex,
                               data=survivalData)
fit.geo.adj.cd <- survfit(surv ~ census_division + 
                                 Size_Class_at_HLES + Breed_Class + Sex,
                               data=survivalData)
fit.geo.adj.st <- survfit(surv ~ state + 
                                 Size_Class_at_HLES + Breed_Class +Sex,
                               data=survivalData)

cox.geo.adj.ha <- coxph(surv ~ home_area +
                                 strata(Size_Class_at_HLES,Breed_Class,Sex),
                             data=survivalData)
cox.geo.adj.cd <- coxph(surv ~ census_division +
                               strata(Size_Class_at_HLES,Breed_Class,Sex),
                             data=survivalData)
cox.geo.adj.st <- coxph(surv ~ state +
                               strata(Size_Class_at_HLES,Breed_Class,Sex),
                             data=survivalData)

cox.geo.adj.pvalues <- data.frame(Analysis="Adjusted",
                                  Variable=c("home_area","census_division","state"),
                                    log.rank.pvalue=c(summary(cox.geo.adj.ha)$sctest["pvalue"],
                                                      summary(cox.geo.adj.cd)$sctest["pvalue"],
                                                      summary(cox.geo.adj.st)$sctest["pvalue"]),
                                  cox.zph.pvalue=c(
                                    cox.zph(cox.geo.adj.ha)$table[1,"p"],
                                    cox.zph(cox.geo.adj.cd)$table[1,"p"],
                                    cox.zph(cox.geo.adj.st)$table[1,"p"]
                                  )
)
print(cox.geo.adj.pvalues)

hrplt.ha <- ggarrange(
  ggforest2(cox.geo.unadj.ha,data=survivalData,main="Unadjusted Hazard Ratio",xrange=1.5),
  ggforest2(cox.geo.adj.ha,data=survivalData,main="Adjusted Hazard Ratio",xrange=1.5),
  ncol=1)

hrplt.cd <- ggarrange(
  ggforest2(cox.geo.unadj.cd,data=survivalData,main="Unadjusted Hazard Ratio",xrange=1.5,
            noDigits = 3),
  ggforest2(cox.geo.adj.cd,data=survivalData,main="Adjusted Hazard Ratio",xrange=1.5,
            noDigits = 3),
  ncol=1)

hrplt.st <- ggarrange(
  ggforest2(cox.geo.unadj.st,data=survivalData,main="Unadjusted Hazard Ratio\n(relative to mean across states)",xrange=2,
            reference_free = TRUE),
  ggforest2(cox.geo.adj.st,data=survivalData,main="Adjusted Hazard Ratio\n(relative to mean across states)",xrange=2,
            reference_free = TRUE),
  ncol=2)

ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_homearea.HR.pdf"),
       hrplt.ha,
       height=5,width=5,scale=1.25)
ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_censusdistrict.HR.pdf"),
       hrplt.cd,
       height=7,width=5,scale=1.35)
ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_state.HR.pdf"),
       hrplt.st,
       height=7,width=7,scale=2)
write.csv(rbind(cox.geo.unadj.pvalues,
                cox.geo.adj.pvalues),
          file=file.path(resultsfolder,
                         "Geoeffect-pvals.csv"),
          row.names = FALSE)

## Separate geo effect for each stratum

cox.geo.pvalues <- data.frame()
for (Size_Class_now in levels(survivalData$Size_Class_at_HLES)) {
  for (Breed_Class_now in levels(survivalData$Breed_Class)) {
    for (Sex_now in levels(survivalData$Sex)) {
      cox.geo.tmp.ha <- coxph(surv ~ home_area,data=survivalData,
                              subset=Size_Class_at_HLES==Size_Class_now &
                                Breed_Class == Breed_Class_now &
                                Sex == Sex_now)
      cox.geo.tmp.cd <- coxph(surv ~ census_division,data=survivalData,
                              subset=Size_Class_at_HLES==Size_Class_now &
                                Breed_Class == Breed_Class_now &
                                Sex == Sex_now)
      cox.geo.tmp.st <- coxph(surv ~ state,data=survivalData,
                              subset=Size_Class_at_HLES==Size_Class_now &
                                Breed_Class == Breed_Class_now &
                                Sex == Sex_now)
      zp.ha <- NA; zp.ha <- as.numeric(try(cox.zph(cox.geo.tmp.ha)$table[1,"p"]))
      zp.cd <- NA; zp.cd <- as.numeric(try(cox.zph(cox.geo.tmp.cd)$table[1,"p"]))
      zp.st <- NA; zp.st <- as.numeric(try(cox.zph(cox.geo.tmp.st)$table[1,"p"]))
      
      cox.geo.pvalues <- rbind(cox.geo.pvalues,
                               data.frame(Size_Class_at_HLES=Size_Class_now,
                                          Breed_Class=Breed_Class_now,
                                          Sex=Sex_now,
                                          log.rank.pvalue.ha=
                                            summary(cox.geo.tmp.ha)$sctest["pvalue"],
                                          log.rank.pvalue.cd=
                                            summary(cox.geo.tmp.cd)$sctest["pvalue"],
                                          log.rank.pvalue.st=
                                            summary(cox.geo.tmp.st)$sctest["pvalue"],
                                          cox.zph.pvalue.ha=zp.ha,
                                          cox.zph.pvalue.cd=zp.cd,
                                          cox.zph.pvalue.st=zp.st)
                               )
    }
  }
}
print(cox.geo.pvalues)

write.csv(cox.geo.pvalues,
          file=file.path(resultsfolder,
                         "Geoeffect-bystratum-pvals.csv"),
          row.names = FALSE)

save(survivalData,
     cox.geo.adj.ha,
     cox.geo.adj.cd,
     cox.geo.adj.st,
     cox.geo.unadj.pvalues,
     cox.geo.adj.pvalues,
     cox.geo.pvalues,
     file=file.path(resultsfolder,"Geoeffect-Cox-results.Rdata"))
