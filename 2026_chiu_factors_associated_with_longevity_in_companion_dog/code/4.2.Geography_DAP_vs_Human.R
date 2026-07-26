library(tidyverse)
library(survival)
library(usdata)
library(wCorr)
library(emmeans)
datfolder <- "../data"
resultsfolder <- "results"
figfolder <- "figures"

load(file.path(resultsfolder,"Geoeffect-Cox-results.Rdata"))
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')


###### Reference-free dog estimates (all 51 states, incl. WA) ######
## Use effect contrasts (method = "eff") = deviations from the GRAND MEAN.
## This gives EVERY state (incl. WA) a genuine estimate, SE, CI, and p-value
## -- no level pinned to 0. Plain emmeans on a coxph returns log-HRs relative
## to the reference, which zeroes WA; effect contrasts fix that. [1]
em    <- emmeans::emmeans(cox.geo.adj.st, specs = ~ state, weights = "cells")
emc   <- emmeans::contrast(em, method = "eff", infer = c(TRUE, TRUE))
em.df <- as.data.frame(emc)

## Confirm the level-string format before trusting the join.
## Effect contrasts label rows like "WA effect"; emmeans usually returns the
## BARE level, so strip the " effect" suffix and any "state" prefix. [1]
print(head(em.df))

clean_lev <- function(x) {
  x <- as.character(x)
  x <- gsub(" effect$", "", x)   # effect-contrast suffix
  x <- gsub("^state", "", x)     # in case a term prefix survives
  trimws(x)
}

hr.st <- data.frame(
  state.abbr = clean_lev(em.df$contrast),
  logHR      = em.df$estimate,   # deviation from count-weighted overall mean
  Dog.SE     = em.df$SE,         # SE of the effect contrast (non-zero for WA)
  stringsAsFactors = FALSE)
hr.st$Dog.HR <- exp(hr.st$logHR)

## Sanity checks: all 51 states, WA now non-degenerate, no dup/empty labels
stopifnot(nrow(hr.st) == 51)
stopifnot(!any(hr.st$state.abbr == "" | is.na(hr.st$state.abbr)))
stopifnot(!any(duplicated(hr.st$state.abbr)))
stopifnot(all(hr.st$Dog.SE > 0))          # <-- WA must no longer be 0 SE
stopifnot(all(is.finite(1 / hr.st$Dog.SE^2)))  # no infinite weights

# count weights per state (to center the human axis consistently)
n_by_state <- as.data.frame(table(gsub("^state", "", as.character(survivalData$state))))
names(n_by_state) <- c("state.abbr", "n_dogs")
n_by_state$state.abbr <- as.character(n_by_state$state.abbr)

###### Human mortality -> reference-free (mean-centered on log scale) ######
mr.st <- subset(read.csv(file.path(datfolder, "HDPulse_data_export.csv"), skip = 4),
                (FIPS > 0) & (!is.na(FIPS)))
mr.st$state.abbr <- state2abbr(mr.st$State)
names(mr.st)[3]  <- "MR"
mr.st <- subset(mr.st, state.abbr %in% hr.st$state.abbr)
mr.st <- left_join(mr.st, n_by_state, by = "state.abbr")
stopifnot(!any(is.na(mr.st$n_dogs)))

mr_logmean      <- weighted.mean(log(mr.st$MR), w = mr.st$Average.Annual.Count)
mr.st$logMR.dev <- log(mr.st$MR) - mr_logmean
mr.st$MR.dev    <- exp(mr.st$logMR.dev)   # ratio scale, centered at 1

###### Join + checks ######
hr.st <- left_join(hr.st,
                   mr.st[, c("state.abbr", "logMR.dev", "MR.dev")],
                   by = "state.abbr")
stopifnot(!any(is.na(hr.st$logMR.dev)))
stopifnot(nrow(hr.st) == 51)

###### Weighted regression + correlations (same variable throughout) ######
lm.res <- lm(logHR ~ logMR.dev, weights = 1/Dog.SE^2, data = hr.st)
summary(cooks.distance(lm.res)); plot(lm.res, which = 4); plot(lm.res, which = 5)
# Result: Cook's distance all < 0.5 

lm.sum      <- summary(lm.res)
hr.mr.pval  <- signif(lm.sum$coefficients["logMR.dev", 4], 2)
hr.mr.slope <- signif(lm.sum$coefficients["logMR.dev", 1], 2)

hr.mr.r   <- signif(weightedCorr(hr.st$logHR, hr.st$logMR.dev,
                                 weights = 1/hr.st$Dog.SE^2, method = "Pearson"), 2)
hr.mr.rho <- signif(weightedCorr(hr.st$logHR, hr.st$logMR.dev,
                                 weights = 1/hr.st$Dog.SE^2, method = "Spearman"), 2)

###### Plot (both axes reference-free, centered at 1) ######
plt.hr.mr.state <- ggplot(hr.st, aes(x = MR.dev, y = Dog.HR)) +
  geom_hline(yintercept = 1) + geom_vline(xintercept = 1) +
  geom_point(aes(size = 1/Dog.SE^2), shape = 21, stroke = 1.5, color = "grey33") +
  geom_text(aes(label = state.abbr), hjust = 1.3, vjust = -0.5, size = 3) +
  scale_x_log10(breaks = seq(0.6, 1.6, 0.2)) +
  scale_y_log10(breaks = seq(0.6, 1.6, 0.2)) +
  scale_size_continuous(range = c(0.01, 6)) +
  theme_bw() + annotation_logticks(side = "bl") +
  geom_smooth(aes(x = MR.dev, y = Dog.HR, weight = 1/Dog.SE^2), method = "lm") +
  annotate("text", x = 0.75, y = 1.4,
           label = list(bquote(italic(r) == .(hr.mr.r) * "," ~ italic(rho) == .(hr.mr.rho))),
           hjust = 0, vjust = 1.5, parse = FALSE) +
  annotate("text", x = 0.75, y = 1.4,
           label = list(bquote(italic(p) == .(hr.mr.pval))), hjust = 0, vjust = 3) +
  annotate("text", x = 0.75, y = 1.4,
           label = list(bquote(slope == .(hr.mr.slope))), hjust = 0, vjust = 4.5) +
  xlab("Human age-adjusted mortality ratio\n(relative to overall mean, 2019–2023)") +
  ylab("Canine adjusted mortality hazard ratio\n(relative to overall mean)") +
  theme(panel.grid.minor = element_blank(), legend.position = "none")
print(plt.hr.mr.state)

ggsave(file.path(figfolder,"Fig.2.Geoeffect_State_DogvsHumanMortRate.pdf"),
       plt.hr.mr.state,height=4,width=4,scale=1.33)

###### Sensitivity Analysis - Use only 2023 mortality rates

mr23.st <- subset(read.csv(file.path(datfolder,"HDPulse_data_export-2023.csv"),skip=4),
                (FIPS > 0 ) & (!is.na(FIPS)))
mr23.st$state.abbr <- state2abbr(mr23.st$State)
names(mr23.st)[3] <- "MR"
mr23.st <- subset(mr23.st,state.abbr %in% hr.st$state.abbr)
mr23.st <- left_join(mr23.st, n_by_state, by = "state.abbr")
stopifnot(!any(is.na(mr23.st$n_dogs)))
mr23_logmean      <- weighted.mean(log(mr23.st$MR), w = mr23.st$Average.Annual.Count)
mr23.st$logMR.dev <- log(mr23.st$MR) - mr23_logmean
mr23.st$MR.dev    <- exp(mr23.st$logMR.dev)   # ratio scale, centered at 1

###### Join + checks ######
hr23.st <- left_join(hr.st[, c("state.abbr", "logHR", "Dog.SE","Dog.HR")],
                   mr23.st[, c("state.abbr", "logMR.dev", "MR.dev")],
                   by = "state.abbr")
stopifnot(!any(is.na(hr23.st$logMR.dev)))
stopifnot(nrow(hr23.st) == 51)

###### Weighted regression + correlations (same variable throughout) ######
lm23.res <- lm(logHR ~ logMR.dev, weights = 1/Dog.SE^2, data = hr23.st)
summary(cooks.distance(lm23.res)); plot(lm23.res, which = 4); plot(lm23.res, which = 5)
# Result: Cook's distance all < 0.5 

lm23.sum      <- summary(lm23.res)
hr.mr23.pval  <- signif(lm23.sum$coefficients["logMR.dev", 4], 2)
hr.mr23.slope <- signif(lm23.sum$coefficients["logMR.dev", 1], 2)

hr.mr23.r   <- signif(weightedCorr(hr23.st$logHR, hr23.st$logMR.dev,
                                 weights = 1/hr.st$Dog.SE^2, method = "Pearson"), 2)
hr.mr23.rho <- signif(weightedCorr(hr23.st$logHR, hr23.st$logMR.dev,
                                 weights = 1/hr.st$Dog.SE^2, method = "Spearman"), 2)

###### Plot (both axes reference-free, centered at 1) ######
plt.hr.mr23.state <- ggplot(hr23.st, aes(x = MR.dev, y = Dog.HR)) +
  geom_hline(yintercept = 1) + geom_vline(xintercept = 1) +
  geom_point(aes(size = 1/Dog.SE^2), shape = 21, stroke = 1.5, color = "grey33") +
  geom_text(aes(label = state.abbr), hjust = 1.3, vjust = -0.5, size = 3) +
  scale_x_log10(breaks = seq(0.6, 1.6, 0.2)) +
  scale_y_log10(breaks = seq(0.6, 1.6, 0.2)) +
  scale_size_continuous(range = c(0.01, 6)) +
  theme_bw() + annotation_logticks(side = "bl") +
  geom_smooth(aes(x = MR.dev, y = Dog.HR, weight = 1/Dog.SE^2), method = "lm") +
  annotate("text", x = 0.75, y = 1.4,
           label = list(bquote(italic(r) == .(hr.mr23.r) * "," ~ italic(rho) == .(hr.mr23.rho))),
           hjust = 0, vjust = 1.5, parse = FALSE) +
  annotate("text", x = 0.75, y = 1.4,
           label = list(bquote(italic(p) == .(hr.mr23.pval))), hjust = 0, vjust = 3) +
  annotate("text", x = 0.75, y = 1.4,
           label = list(bquote(slope == .(hr.mr23.slope))), hjust = 0, vjust = 4.5) +
  xlab("Human age-adjusted mortality ratio\n(relative to overall mean, 2023)") +
  ylab("Canine adjusted mortality hazard ratio\n(relative to overall mean)") +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

print(plt.hr.mr23.state)

ggsave(file.path(figfolder,"SuppFig.Geoeffect_State_DogvsHuman2023MortRate.pdf"),
       plt.hr.mr23.state,height=4,width=4,scale=1.25)

###### Correlation between hazard ratios by states and US life expectancy

le.st <- read.csv(file.path(datfolder,"U.S._State_Life_Expectancy_by_Sex__2021.csv"))
le.st$state.abbr <- state2abbr(le.st$State)
le.st <- subset(le.st,Sex=="Total")
hr.le.st <- left_join(hr.st[, c("state.abbr", "logHR", "Dog.SE","Dog.HR")],
                     le.st,
                     by = "state.abbr")

###### Weighted regression + correlations (same variable throughout) ######
lm.le.res <- lm(Dog.HR ~ LE,weights=1/Dog.SE^2,data=hr.le.st)
summary(cooks.distance(lm.le.res)); plot(lm.le.res, which = 4); plot(lm.le.res, which = 5)
# Result: Cook's distance all < 0.5 

lm.le.sum      <- summary(lm.le.res)
hr.le.pval  <- signif(lm.le.sum$coefficients["LE", 4], 2)
hr.le.slope <- signif(lm.le.sum$coefficients["LE", 1], 2)
hr.le.r   <- signif(weightedCorr(hr.le.st$Dog.HR, hr.le.st$LE,
                                   weights = 1/hr.st$Dog.SE^2, method = "Pearson"), 2)
hr.le.rho <- signif(weightedCorr(hr.le.st$Dog.HR, hr.le.st$LE,
                                   weights = 1/hr.st$Dog.SE^2, method = "Spearman"), 2)

plt.hr.le.state <-
  ggplot(hr.le.st,aes(x=LE,y=Dog.HR))+
  geom_point(aes(size=1/Dog.SE^2),shape=21,stroke=1.5,color="grey33")+
  geom_text(aes(label=state.abbr),hjust=1.3,vjust=-0.5,size=3)+
  scale_x_continuous(breaks=seq(65,85))+
  theme_bw()+coord_transform(y="log10",ylim=c(0.5,2.1))+
  annotation_logticks(side="l",scaled=FALSE)+
  geom_smooth(aes(weight=1/Dog.SE^2),method="lm")+
  geom_abline(slope=-1/76.4,intercept=2,color="red",linetype="dashed")+
  annotate("text",x=71,y=2,label="Linear fit to state-level HRs",hjust=0,vjust=0)+
  annotate("text", x = 71, y = 2,
           label = list(bquote(italic(r) == .(hr.le.r) * "," ~ italic(rho) == .(hr.le.rho))),
           hjust = 0, vjust = 1.5, parse = FALSE) +
  annotate("text", x = 71, y = 2,
           label = list(bquote(italic(p) == .(hr.le.pval))), hjust = 0, vjust = 3) +
  annotate("text", x = 71, y = 2,
           label = list(bquote(slope == .(hr.le.slope))), hjust = 0, vjust = 4.5) +
  xlab("Human Life Expectancy (yrs)\nby State in 2021")+
  ylab("Demographically Adjusted DAP Mortality Hazard Ratio (HR)")+
  theme(panel.grid.minor = element_blank(),legend.position = "none")
print(plt.hr.le.state)

hr.le.slope.rel <- signif(76.4*first(lm.le.sum$coefficients["LE",]),2)

#### Try using life expectancy as a continuous variable

rownames(le.st) <- le.st$state.abbr
survivalData$state.LE <- le.st[as.character(survivalData$state),"LE"]
cox.geo.adj.state.LE <- coxph(surv ~ state.LE +
                                strata(Size_Class_at_HLES,Breed_Class,Sex),
                              data=survivalData)
cox.geo.adj.state.LE.zph <- cox.zph(cox.geo.adj.state.LE)
survivalData$state.LE.5 <- survivalData$state.LE/5
cox.geo.adj.state.LE.5 <- coxph(surv ~ state.LE.5 +
                                  strata(Size_Class_at_HLES,Breed_Class,Sex),
                                data=survivalData)
cox.geo.adj.state.LE.5.sum <- summary(cox.geo.adj.state.LE.5)
exp(c(coef(cox.geo.adj.state.LE.5),confint(cox.geo.adj.state.LE.5)))

cox.geo.adj.state.LE.5.coef <- 
  cbind(data.frame(Variable="State-Level\nHuman Life Expectancy",
                   logrankp=cox.geo.adj.state.LE.5.sum$sctest["pvalue"]),
        left_join(as.data.frame(cox.geo.adj.state.LE.5.sum$coefficients),
                  as.data.frame(cox.geo.adj.state.LE.5.sum$conf.int)),
        as.data.frame(t(cox.geo.adj.state.LE.zph$table[1,]))
  )
names(cox.geo.adj.state.LE.5.coef)[ncol(cox.geo.adj.state.LE.5.coef)-(2:0)] <- 
  paste0(names(cox.geo.adj.state.LE.5.coef)[ncol(cox.geo.adj.state.LE.5.coef)-(2:0)],".zph")
write.csv(cox.geo.adj.state.LE.5.coef,
          file=file.path(resultsfolder,
                         "Geoeffect-Cox-Human.LE.csv"),
          row.names = FALSE)

hr.le.5.hr <- signif(cox.geo.adj.state.LE.5.coef$`exp(coef)`,2)
hr.le.5.hr.lcl <- signif(cox.geo.adj.state.LE.5.coef$`lower .95`,2)
hr.le.5.hr.ucl <- signif(cox.geo.adj.state.LE.5.coef$`upper .95`,2)
hr.le.5.pval <- as.numeric(signif(cox.geo.adj.state.LE.5.coef$logrankp,2))
plt.hr.le.5 <- 
  ggplot(cox.geo.adj.state.LE.5.coef)+
  geom_hline(yintercept=1,linetype="dotted")+
  geom_point(aes(y=`exp(coef)`,x=Variable),shape=15,size=3,color="grey50")+
  geom_errorbar(aes(ymin=`lower .95`,ymax=`upper .95`,x=Variable),width=0)+
  annotate("text",x=0.5,y=2,label="Cox model for state",
           hjust=0,vjust=0)+
  annotate("text",x=0.5,y=2,label="human life expectancy",
           hjust=0,vjust=1.5)+
  annotate("text",x=0.5,y=2,label=
             list(bquote(HR == .(hr.le.5.hr)~"["*.(hr.le.5.hr.lcl) - .(hr.le.5.hr.ucl)*"]")),
           hjust=0,vjust=3)+
  annotate("text",x=0.5,y=2,label=
             list(bquote(italic(p) == .(hr.le.5.pval))),
           hjust=0,vjust=4.5)+
  ylab("")+
  xlab("per 5 yr increase in\nHuman Life Expectancy ")+
  scale_x_discrete(labels="")+
  theme_bw()+coord_trans(y="log10",ylim=c(0.4,2.1))+
  annotation_logticks(side="l",scaled=FALSE)+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x=element_blank())
print(plt.hr.le.5)
plt.hr.le.state.comb<-ggarrange(plt.hr.le.state,plt.hr.le.5,ncol=2,widths = c(2,1),labels = "auto")

ggsave(file.path(figfolder,"SuppFig.Geoeffect_State_DogvsHumanLE.pdf"),
       plt.hr.le.state.comb,height=4,width=6,scale=1.25)
