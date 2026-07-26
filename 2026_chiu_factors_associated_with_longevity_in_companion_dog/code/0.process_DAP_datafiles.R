library(lubridate)
library(dplyr)
datfolder <- "../data"
survivalfile <- "SurvivalAnalysisCuratedDogs_thru_2024_ran_2025_12_16.csv" # "survivalData_2024_12_15a.csv"

survivalData <- read.csv(file.path(
  datfolder,survivalfile)) %>% 
  rename(dog_id=study_id) %>% left_join(
    read.csv(file.path(datfolder,"DAP_2024_DogOverview_v1.0.csv"))) 
datecols <- which((grepl("date",names(survivalData),ignore.case = TRUE) |
                     grepl("dob",names(survivalData),ignore.case = TRUE) |
                     grepl("dod",names(survivalData),ignore.case = TRUE)) &
                    !(grepl("method",names(survivalData),ignore.case = TRUE)) &
                    !(grepl("source",names(survivalData),ignore.case = TRUE)) &
                    !(grepl("certainty",names(survivalData),ignore.case = TRUE)))
survivalData <- survivalData %>%
  mutate(across(all_of(datecols), ~ as.Date(.,tryFormats =c("%Y-%m-%d","%m-%d-%Y", "%m/%d/%Y"))))

# Collapse male and female
survivalData$Sex <- "Female"
survivalData$Sex[grepl("Male",survivalData$Sex_Class_at_HLES)] <- "Male"
survivalData$Sex <- factor(survivalData$Sex,
                           levels=c("Female","Male"))

# Rename weight and breed classes
survivalData$Weight_Class_10KGBin_at_HLES <-
  paste(survivalData$Weight_Class_10KGBin_at_HLES,"kg")
survivalData$Weight_Class_10KGBin_at_HLES <- 
  factor(survivalData$Weight_Class_10KGBin_at_HLES,
         levels=sort(unique(survivalData$Weight_Class_10KGBin_at_HLES)))
survivalData$Breed_Class <- factor(survivalData$Breed_Class,
                                   levels=c("Non-AKC-Recognized or Mixed Breed",
                                            "AKC-Recognized Breed"))
survivalData$Size_Class_at_HLES <- survivalData$Breed_Size_Class_at_HLES
survivalData$Size_Class_at_HLES <- gsub("AKC-recognized purebred ","",
                                        survivalData$Size_Class_at_HLES)
survivalData$Size_Class_at_HLES <- gsub("non-AKC and mixed breed ","",
                                        survivalData$Size_Class_at_HLES)
survivalData$Size_Class_at_HLES <- factor(survivalData$Size_Class_at_HLES,
                                          levels=
                                            rev(c("Giant dogs",
                                                  "Large dogs",
                                                  "Standard dogs",
                                                  "Medium dogs",
                                                  "Toy and Small dogs")
                                            ))

survivalData$Breed_Size_Class_at_HLES <- factor(survivalData$Breed_Size_Class_at_HLES,
                                                levels=
                                                  c("Giant AKC-recognized purebred dogs",
                                                    "Giant non-AKC and mixed breed dogs",
                                                    "Large AKC-recognized purebred dogs",
                                                    "Large non-AKC and mixed breed dogs",
                                                    "Standard AKC-recognized purebred dogs",
                                                    "Standard non-AKC and mixed breed dogs",
                                                    "Medium AKC-recognized purebred dogs",
                                                    "Medium non-AKC and mixed breed dogs",
                                                    "Toy and Small AKC-recognized purebred dogs",
                                                    "Toy and Small non-AKC and mixed breed dogs"
                                                  ))

if ("Estimated_DOB" %in% names(survivalData) & 
    (!("st_estimated_dob" %in% names(survivalData)))) {
  survivalData <- rename(survivalData,st_estimated_dob=Estimated_DOB)
}
# "estimated_baseline_age" is a rounded version of first.age
survivalData$first.age <- time_length(interval(survivalData$st_estimated_dob,survivalData$dap_pack_date),"years")
# "estimated_survival_age" is a rounded version of last.age
survivalData$last.age <- time_length(interval(survivalData$st_estimated_dob,survivalData$survival_date),"years")
survivalData$event <- as.numeric(survivalData$survival_status=="Dead")

survivalData <- subset(survivalData,last.age > first.age)

save(survivalData,file=file.path(datfolder,"SurvivalData.RData"))
write.csv(survivalData,file.path(datfolder,"SurvivalData.csv"),row.names = FALSE)


