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
library(corrplot)
library(xtable)
library(ggpubr)
library(car)
library(dplyr)
library(effectsize)

# get data
source("~/Documents/Projects/Degeneration/code/cerebellar_degeneration/r/data_connectivity.R")

table_path <- "~/Documents/Projects/Degeneration/derivatives/results/tables/"



# ---------------------- Baseline connectivity differs between groups --------------------------
modelRIM.Interaction      = lmer(connectivity ~  group*region_type + (1 | subject), data=base_connectivity, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)

# ---------------------- Change in left PMD - right Cerebellar connectivity depends on vision ---------[specific]-----------------
# Region of interest
region <- "pmd-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)


# ---------------------- Change in right PMD - right Cerebellar connectivity depends on vision ---------[specific]-----------------
# Region of interest
region <- "pmd-right_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)
# ---------------------- Connectivity change in left PPC - right cerebellum motor depends on vision --------------[specific]-------
region <- "ppc-left_cereb-M3-right"
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)
# ---------------------- Connectivity change in right PPC - left cerebellum motor ---------------------
region <- "ppc-right_cereb-M3-left"
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)

# ---------------------- Connectivity increases beween M1 & Cerebellum (hand area) in all participants - expectedly contralateral to the moving hand --------------------------
# Regions of interest
region <- "m1-hand-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)

region <- "m1-hand-left_m1-hand-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)

region <- "cereb-M3-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)
# ---------------------- Connectivity change in right dlPFC - right D2 depends on feedback and group -----------[specific]-------
region <- "dlpfc-right_cereb-D2-right"
modelRIM.Interaction      = lmer(connectivity ~  group*session*feedback*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)

region <- "dlpfc-right_cereb-D2-left"
modelRIM.Interaction      = lmer(connectivity ~  group*session*feedback*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

std_model <- standardize_parameters(modelRIM.Interaction)
print(std_model)