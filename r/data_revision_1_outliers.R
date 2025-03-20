# Packages
library(readxl)
library(tidyr)
library(dplyr)
library(DescTools)
library(ggpubr)

# ---------------------- Import data_no_outliers and settings --------------------------
# connectivity
data_path='~/Documents/Projects/Degeneration/derivatives/'
savepath = paste(data_path, 'figures/', sep="")
file= '../code/cerebellar_degeneration/data/connectivity_no_outliers.tsv'
data_no_outliers <- data.frame(read.csv(file.path(data_path, file), sep = "\t"))


data_no_outliers$subject <- as.factor(data_no_outliers$subject)
data_no_outliers$session <- factor(data_no_outliers$session, levels = c("ses-pre", "ses-post"))
data_no_outliers$group <- as.factor(data_no_outliers$group)
data_no_outliers$condition <- as.factor(data_no_outliers$condition)
data_no_outliers$sex <- as.factor(data_no_outliers$sex)
data_no_outliers$regions <- as.factor(data_no_outliers$regions)

unlist(lapply(data_no_outliers, is.numeric)) # Check that the correct variables are numeric (all apart from subject, session, group)
head(data_no_outliers)




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
data_no_outliers$laterality <- NA
data_no_outliers[data_no_outliers$regions %in% regions_ipsi,]$laterality <- c("ipsi")
data_no_outliers[data_no_outliers$regions %in% c(regions_contra, regions_cortico, regions_cereb),]$laterality <- c("contra")
data_no_outliers$laterality <- as.factor(data_no_outliers$laterality)

# ---- Create cortico-cortico / cerebello-cerebello / cortico-cerebello index ----
data_no_outliers$interregional <- NA
data_no_outliers[data_no_outliers$regions %in% regions_cortico,]$interregional <- c("cortico-cortico")
data_no_outliers[data_no_outliers$regions %in% c(regions_contra, regions_ipsi),]$interregional <- c("cortico-cerebello")
data_no_outliers[data_no_outliers$regions %in% regions_cereb,]$interregional <- c("cerebello-cerebello")
data_no_outliers$interregional <- factor(data_no_outliers$interregional, levels = c("cortico-cortico", "cerebello-cerebello", "cortico-cerebello"))

# ---- Create dataframes for connectivity at baseline (pre) and connectivity change ----
dropvars <- names(diff) %in% demographic_vars

# ---- Apply Fisher's Z-transform to each region's connectivity values ----
data_non_norm <- data_no_outliers
base_connectivity_non_norm <- base_connectivity
for (region in levels(data_no_outliers$regions)) {
  data_no_outliers[data_no_outliers$regions == region,]$connectivity <- FisherZ(data_no_outliers[data_no_outliers$regions == region,]$connectivity)
}

