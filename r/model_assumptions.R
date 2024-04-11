# --------------------------------------------------------------------
# Packages
library(readxl)
library(lme4)
library(lmerTest)
library(Hmisc)
library(tidyr)
# Packages to plot model
library(sjPlot)
theme_set(theme_sjplot())
library(sjmisc)
library(ggplot2)
# package for post hoc tests
library(emmeans)
library(corrplot)
library(xtable)
library(ggpubr)
library(car)

# get data
source("~/Documents/Projects/Degeneration/code/cerebellar_degeneration/r/data_connectivity.R")


# ---- Check model assumptions ----
# Homoscedasticity: Errors are the same across all values of the independent variables
# Check by plotting residuals against predicted values
plot(modelRIM.Interaction)
res <- residuals(modelRIM.Interaction)
boxplot(res ~ base_connectivity$group)
boxplot(res ~ base_connectivity$interregional)
boxplot(res ~ base_connectivity$group * base_connectivity$interregional)

# Normality of error term: check histogram and Q-Qplot of residuals.
ggdensity(res)
qqnorm(res)
qqline(res) 

# Normality of random effect: Get the estimate of random effect (random intercept of subject), and check normality
random_effects <- ranef(modelRIM.Interaction)
random_effects_vector <- as.numeric(random_effects$subject$`(Intercept)`)
ggdensity(random_effects_vector)
qqnorm(random_effects_vector)
qqline(random_effects_vector)
