# Packages
library(readxl)
library(tidyr)
library(dplyr)
library(DescTools)
library(ggpubr)

# ---------------------- Import data and settings --------------------------
# connectivity
data_path='~/Documents/Projects/Degeneration/derivatives/'
savepath = paste(data_path, 'figures/', sep="")
file= '../code/cerebellar_degeneration/data/connectivity.tsv'
data <- data.frame(read.csv(file.path(data_path, file), sep = "\t"))

data$subject <- as.factor(data$subject)
data$session <- factor(data$session, levels = c("ses-pre", "ses-post"))
data$group <- as.factor(data$group)
data$condition <- as.factor(data$condition)
data$sex <- as.factor(data$sex)
data$regions <- as.factor(data$regions)

unlist(lapply(data, is.numeric)) # Check that the correct variables are numeric (all apart from subject, session, group)
head(data)

# ---------------------- Variable lists --------------------------
regions_contra = c(
  'ppc-right_cereb-M3-left',
  'ppc-left_cereb-M3-right',
  'pmd-left_cereb-M3-right',
  'pmd-right_cereb-M3-left',
  'm1-hand-left_cereb-M3-right',
  'm1-hand-right_cereb-M3-left',
  'dlpfc-right_cereb-D2-left')


regions_ipsi = c(
  'ppc-right_cereb-M3-right',
  'ppc-left_cereb-M3-left',
  'pmd-right_cereb-M3-right',
  'pmd-left_cereb-M3-left',
  'm1-hand-right_cereb-M3-right',
  'm1-hand-left_cereb-M3-left',
  'dlpfc-right_cereb-D2-right')

regions_cortico = c(
  'ppc-left_ppc-right',
  'pmd-left_pmd-right',
  'm1-hand-left_m1-hand-right')

regions_cereb = c(
  'cereb-D2-left_cereb-D2-right',
  'cereb-M3-left_cereb-M3-right')

ppc_regions = c(
  'ppc-right_cereb-M3-left',
  'ppc-left_cereb-M3-right',
  'ppc-right_cereb-M3-right',
  'ppc-left_cereb-M3-left',
  'ppc-left_ppc-right')

pmd_regions = c(
  'pmd-left_cereb-M3-right',
  'pmd-right_cereb-M3-left',
  'pmd-left_pmd-right')

m1_regions = c(
  'm1-hand-left_cereb-M3-right',
  'm1-hand-right_cereb-M3-left',
  'm1-hand-left_m1-hand-right')

dlpfc_regions = c(
  'dlpfc-right_cereb-D2-left',
  'dlpfc-right_cereb-D2-right')

cereb_regions = c(
  'cereb-D2-left_cereb-D2-right',
  'cereb-M3-left_cereb-M3-right')


demographic_vars = c('age', 'sex', 'group', 'condition', 'vision', 'feedback')

# ---- Create laterality (ipsi or contra) index ----
data$laterality <- NA
data[data$regions %in% regions_ipsi,]$laterality <- c("ipsi")
data[data$regions %in% c(regions_contra, regions_cortico, regions_cereb),]$laterality <- c("contra")
data$laterality <- as.factor(data$laterality)

# ---- Create cortico-cortico / cerebello-cerebello / cortico-cerebello index ----
data$interregional <- NA
data[data$regions %in% regions_cortico,]$interregional <- c("cortico-cortico")
data[data$regions %in% c(regions_contra, regions_ipsi),]$interregional <- c("cortico-cerebello")
data[data$regions %in% regions_cereb,]$interregional <- c("cerebello-cerebello")
data$interregional <- factor(data$interregional, levels = c("cortico-cortico", "cerebello-cerebello", "cortico-cerebello"))

# ---- Create dataframes for connectivity at baseline (pre) and connectivity change ----
dropvars <- names(diff) %in% demographic_vars

# ---- Create dataframe for mean subject cortico-cortico, cortico-cerebello and cerebello-cerebello connectivity ----
# Select connectivity ROIs
connectivity_pre_select <- connectivity_pre[connectivity_pre$regions %in% c(pmd_regions, ppc_regions, m1_regions, cereb_regions),]

# Group by subject and region_type factors
base_connectivity <- aggregate(x=connectivity_pre_select$connectivity,
          by=list(connectivity_pre_select$subject,connectivity_pre_select$group,connectivity_pre_select$interregional),
          FUN=mean) 

names(base_connectivity) <- c("subject", "group", "region_type", "connectivity")

# ---- Apply Fisher's Z-transform to each region's connectivity values ----
data_non_norm <- data
base_connectivity_non_norm <- base_connectivity
for (region in levels(data$regions)) {
  data[data$regions == region,]$connectivity <- FisherZ(data[data$regions == region,]$connectivity)
}

base_connectivity$connectivity <- FisherZ(base_connectivity$connectivity)
