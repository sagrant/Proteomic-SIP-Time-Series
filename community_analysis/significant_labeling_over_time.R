### Load libraries 
library(dplyr) 

proportionsData <- read.csv("labeled_proportions.csv")

cecum = proportionsData %>% filter(Sample_Type == "C")
stool = proportionsData %>% filter(Sample_Type == "S")


both.GLM.fit <- glm(Proportion_Labeled ~ Sample_Type,
                 data = proportionsData,
                 family = quasibinomial,
                 weights = Total_PSMs)
summary(both.GLM.fit)

cecum.GLM.fit <- glm(Proportion_Labeled ~ Time_Point,
                 data = cecum,
                 family = quasibinomial,
                 weights = Total_PSMs)
summary(cecum.GLM.fit)

stool.GLM.fit <- glm(Proportion_Labeled ~ Time_Point,
                 data = stool,
                 family = quasibinomial,
                 weights = Total_PSMs)
summary(stool.GLM.fit)
