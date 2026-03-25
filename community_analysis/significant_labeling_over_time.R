### Load libraries 
library(dplyr) 

proportionsData <- read.csv("labeled_proportions.csv")

### Subset data by sample type (cecum vs. stool)
cecum = proportionsData %>% filter(Sample_Type == "C")
stool = proportionsData %>% filter(Sample_Type == "S")

### Are there significant differences between sample types?
both.GLM.fit <- glm(Proportion_Labeled ~ Sample_Type,
                 data = proportionsData,
                 family = quasibinomial,
                 weights = Total_PSMs)
summary(both.GLM.fit)

### Are there significant differences in labeling over time in cecum samples?
cecum.GLM.fit <- glm(Proportion_Labeled ~ Time_Point,
                 data = cecum,
                 family = quasibinomial,
                 weights = Total_PSMs)
summary(cecum.GLM.fit)

### Are there significant differences in labeling over time in stool samples?
stool.GLM.fit <- glm(Proportion_Labeled ~ Time_Point,
                 data = stool,
                 family = quasibinomial,
                 weights = Total_PSMs)
summary(stool.GLM.fit)
