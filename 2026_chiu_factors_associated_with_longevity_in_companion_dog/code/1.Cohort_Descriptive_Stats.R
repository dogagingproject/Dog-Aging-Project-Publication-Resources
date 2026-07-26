library(lubridate)
library(survival)
library(flexsurv)
library(survminer)
library(dplyr)
library(tidyverse)
library(ggpubr)
datfolder <- "../data"
resultsfolder <- "results"
figfolder <- "figures"
load(file.path(datfolder,"SurvivalData.RData"))
survivalData$death.age <- survivalData$last.age
survivalData$death.age[survivalData$event==0] <- NA  
survivalData$followup.years <- survivalData$last.age - survivalData$first.age

tot_all <- aggregate(event~1,data=survivalData,length)
F_all <- aggregate(event~1,
                     data=subset(survivalData,Sex=="Female"),length)
M_all <- aggregate(event~1,
                     data=subset(survivalData,Sex=="Male"),length)
deaths_all<-aggregate(event~1,data=survivalData,sum)
crude_deaths_all<-data.frame(Breed_Class="Total",
                             N=tot_all[[1]],
                             F=F_all[[1]],
                             M=M_all[[1]],
                             Deaths=deaths_all[[1]])
ages_all.mean<-data.frame(Breed_Class="Total",
                          aggregate(first.age~1,data=survivalData,mean,na.rm=T),
                          aggregate(followup.years~1,data=survivalData,mean,na.rm=T),
                          aggregate(death.age~1,data=survivalData,mean,na.rm=T))
crude_deaths_all$mort.rate <- crude_deaths_all$Deaths/(
  crude_deaths_all$N*ages_all.mean$followup.years)

ages_all.quant<-as.data.frame(cbind(data.frame(Breed_Class="Total"),
                                 aggregate(first.age~1,data=survivalData,quantile,na.rm=T)[1],
                      aggregate(followup.years~1,data=survivalData,quantile,na.rm=T)[1],
                      aggregate(death.age~1,data=survivalData,quantile,na.rm=T)[1]
                      ))
dap.desc.stat.all <- left_join(crude_deaths_all,
                                 left_join(
                                   rename(ages_all.mean,
                                          first.age.mean=first.age,
                                          followup.years.mean=followup.years,
                                          death.age.mean=death.age),
                                   ages_all.quant)
)


tot_breed <- aggregate(event~Breed_Class,data=survivalData,length)
F_breed <- aggregate(event~Breed_Class,
                     data=subset(survivalData,Sex=="Female"),length)
M_breed <- aggregate(event~Breed_Class,
                     data=subset(survivalData,Sex=="Male"),length)
deaths_breed<-aggregate(event~Breed_Class,data=survivalData,sum)
crude_deaths_breed<-left_join(left_join(rename(tot_breed,N=event),
                                        left_join(rename(F_breed,F=event),
                                                  rename(M_breed,M=event))),
                              rename(deaths_breed,Deaths=event))
ages_breed.mean<-left_join(
  left_join(aggregate(first.age~Breed_Class,
                      data=survivalData,mean,na.rm=T),
            aggregate(followup.years~Breed_Class,
                      data=survivalData,mean,na.rm=T)
  ),aggregate(death.age~Breed_Class,
              data=survivalData,mean,na.rm=T)
)
crude_deaths_breed$mort.rate <- crude_deaths_breed$Deaths/(
  crude_deaths_breed$N*ages_breed.mean$followup.years)

ages_breed.quant<-left_join(
  left_join(aggregate(first.age~Breed_Class,
                      data=survivalData,quantile,na.rm=T),
            aggregate(followup.years~Breed_Class,
                      data=survivalData,quantile,na.rm=T)
  ),aggregate(death.age~Breed_Class,
              data=survivalData,quantile,na.rm=T)
)

dap.desc.stat.breed <- left_join(crude_deaths_breed,
                                 left_join(
                                   rename(ages_breed.mean,
                                          first.age.mean=first.age,
                                          followup.years.mean=followup.years,
                                          death.age.mean=death.age),
                                   ages_breed.quant)
)

# Deaths / N different by Breed_Class?
dap.desc.stat.breed.deaths.test <- prop.test(dap.desc.stat.breed$Deaths,
                                                dap.desc.stat.breed$N)
print(dap.desc.stat.breed.deaths.test$p.value)

# Deaths / follow-up in days different by Breed_Class?
dap.desc.stat.breed.crudemort.test <- prop.test(dap.desc.stat.breed$Deaths,
                                                dap.desc.stat.breed$N*
                                                  dap.desc.stat.breed$followup.years.mean*365)
print(dap.desc.stat.breed.crudemort.test$p.value)

# Entry age different by Breed_Class?
dap.desc.stat.breed.entry.age.test<-summary(lm(first.age~Breed_Class,data=survivalData))
print(dap.desc.stat.breed.entry.age.test$coefficients)

write.csv(rbind(dap.desc.stat.all,dap.desc.stat.breed),
          file.path(resultsfolder,"DAP_desc_stat.csv"),row.names = FALSE)

### Strata

tot_strata<-aggregate(event~Size_Class_at_HLES*Breed_Class*Sex,data=survivalData,length)
deaths_strata<-aggregate(event~Size_Class_at_HLES*Breed_Class*Sex,data=survivalData,sum)
crude_deaths<-left_join(rename(tot_strata,N=event),
                        rename(deaths_strata,Deaths=event))

ages_strata.mean<-left_join(
  left_join(aggregate(first.age~Size_Class_at_HLES*Breed_Class*Sex,
                      data=survivalData,mean,na.rm=T),
            aggregate(followup.years~Size_Class_at_HLES*Breed_Class*Sex,
                      data=survivalData,mean,na.rm=T)
  ),aggregate(death.age~Size_Class_at_HLES*Breed_Class*Sex,
              data=survivalData,mean,na.rm=T)
)
crude_deaths$mort.rate <- crude_deaths$Deaths/(
  crude_deaths$N*ages_strata.mean$followup.years)

ages_strata.quant<-left_join(
  left_join(aggregate(first.age~Size_Class_at_HLES*Breed_Class*Sex,
                      data=survivalData,quantile,na.rm=T),
            aggregate(followup.years~Size_Class_at_HLES*Breed_Class*Sex,
                      data=survivalData,quantile,na.rm=T)
  ),aggregate(death.age~Size_Class_at_HLES*Breed_Class*Sex,
              data=survivalData,quantile,na.rm=T)
)

dap.desc.stat.strata <- left_join(crude_deaths,
                                  left_join(
                                    rename(ages_strata.mean,
                                           first.age.mean=first.age,
                                           followup.years.mean=followup.years,
                                           death.age.mean=death.age),
                                    ages_strata.quant)
)
write.csv(dap.desc.stat.strata,file.path(resultsfolder,"DAP_desc_stat_strata.csv"),row.names = FALSE)

write.csv(as.data.frame(summary(aov(first.age~Size_Class_at_HLES*Breed_Class*Sex,
                                    data=survivalData))[[1]]),
          file=file.path(resultsfolder,"Supp.AOV_strata.first.age.csv"))
write.csv(as.data.frame(summary(aov(followup.years~Size_Class_at_HLES*Breed_Class*Sex,
                                    data=survivalData))[[1]]),
          file=file.path(resultsfolder,"Supp.AOV_strata.followup.years.csv"))

p1<-
  ggplot(survivalData)+
  geom_boxplot(aes(x=first.age,y=Size_Class_at_HLES),outlier.size = 0.5)+
  scale_y_discrete(limits=rev)+facet_grid(Sex~Breed_Class)+
  theme_bw()+xlab("Age (yr) at entry")
p2<-
  ggplot(survivalData)+
  geom_boxplot(aes(x=followup.years,y=Size_Class_at_HLES),outlier.size = 0.5)+
  scale_y_discrete(limits=rev)+facet_grid(Sex~Breed_Class)+
  theme_bw()+xlab("Years of follow-up")
p3<-
  ggplot(survivalData)+
  geom_boxplot(aes(x=death.age,y=Size_Class_at_HLES),outlier.size = 0.5)+
  scale_y_discrete(limits=rev)+facet_grid(Sex~Breed_Class)+
  theme_bw()+xlab("Age (yr) of death")
pall <- ggarrange(p1,p2,p3,ncol=1)
print(pall)
ggsave(file.path(figfolder,"Supp.Fig.ages and followup by strata.pdf"),
       pall,scale=1.5)
