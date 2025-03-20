# Packages
library(readxl)
library(tidyr)
library(dplyr)
library(DescTools)
library(ggpubr)

# ---------------------- Import data_control and settings --------------------------
# connectivity
data_path='~/Documents/Projects/Degeneration/derivatives/'
savepath = paste(data_path, 'figures/', sep="")
file= '../code/cerebellar_degeneration/data/connectivity_control.tsv'
data_control <- data.frame(read.csv(file.path(data_path, file), sep = "\t"))


data_control$subject <- as.factor(data_control$subject)
data_control$session <- factor(data_control$session, levels = c("ses-pre", "ses-post"))
data_control$group <- as.factor(data_control$group)
data_control$condition <- as.factor(data_control$condition)
data_control$sex <- as.factor(data_control$sex)
data_control$regions <- as.factor(data_control$regions)

unlist(lapply(data_control, is.numeric)) # Check that the correct variables are numeric (all apart from subject, session, group)
head(data_control)


# ---------------------- Variable lists --------------------------
demographic_vars = c('age', 'sex', 'group', 'condition', 'vision', 'feedback')



# ---- Create cortico-cortico / cerebello-cerebello / cortico-cerebello index ----
data_control$interregional <- c("cortico-cerebello")

# ---- Create dataframes for connectivity at baseline (pre) and connectivity change ----
dropvars <- names(diff) %in% demographic_vars

# ---- Apply Fisher's Z-transform to each region's connectivity values ----
data_non_norm <- data_control
for (region in levels(data_control$regions)) {
  data_control[data_control$regions == region,]$connectivity <- FisherZ(data_control[data_control$regions == region,]$connectivity)
}

