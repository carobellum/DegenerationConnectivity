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
library(dplyr)

# get data
source("~/Documents/Projects/Degeneration/code/cerebellar_degeneration/r/data_connectivity.R")

table_path <- "~/Documents/Projects/Degeneration/derivatives/results/tables/"
# ---------------------- Baseline connectivity differs between groups --------------------------
modelRIM.Interaction      = lmer(connectivity ~  group*region_type + (1 | subject), data=base_connectivity, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)
# Save table as latex table
latex_anova <- xtable(aov, caption = "Linear mixed effects model for connectivity at baseline")
file_path <- paste(table_path, "baseline.tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)

# model prediction plot
plot_model(modelRIM.Interaction, type = "pred", terms = c("region_type", "group"))

# ---- Post hoc tests : Test for patient-control differences in baseline connectivity
region_types <- c("cortico-cortico", "cerebello-cerebello", "cortico-cerebello")
p_values <- c()
results <- list()
for (i in seq_along(region_types)) {
  inter <- region_types[[i]]
  # Subset the data
  x <- base_connectivity[base_connectivity$region_type == inter & base_connectivity$group == 'p', ]$connectivity
  y <- base_connectivity[base_connectivity$region_type == inter & base_connectivity$group == 'c', ]$connectivity
  # Perform the t-test
  var.test(x,y)
  result <- t.test(x, y, paired = FALSE, var.equal = TRUE)
  results[[i]] <- result
  # Add the result to the list
  p_values <- c(p_values, result$p.value)
}
# Correct the post hoc tests for multiple comparisons
corrected_p_values <- p.adjust(p_values, method="BY")
# Print post hoc tests with multiple-comparisons-corrected p-values
for (i in seq_along(region_types)) {
  inter <- region_types[[i]]
  corrected_p_value <- corrected_p_values[i]
  
  if (corrected_p_value < 0.05) {
    cat("Significant region_type difference found for:", paste(inter, collapse=", "), "\n")
    cat("Bonferroni corrected p-value:", corrected_p_value, "\n")
    print(results[[i]])
  }
}
# ---- Post hoc tests : Test for region_type differences in baseline connectivity
region_types <- list(c("cortico-cortico", "cerebello-cerebello"),
                      c("cortico-cortico", "cortico-cerebello"),
                      c("cerebello-cerebello", "cortico-cerebello"))

p_values <- c()
results <- list()
for (i in seq_along(region_types)) {
  inter <- region_types[[i]]
  # Subset the data
  x <- base_connectivity[base_connectivity$region_type == inter[1], ]$connectivity
  y <- base_connectivity[base_connectivity$region_type == inter[2], ]$connectivity
  # Perform the t-test
  var.test(x,y)
  result <- t.test(x, y, paired = FALSE, var.equal = TRUE)
  results[[i]] <- result
  # Add the result to the list
  p_values <- c(p_values, result$p.value)
}

# Correct the post hoc tests for multiple comparisons
corrected_p_values <- p.adjust(p_values, method="BY")
# Print post hoc tests with multiple-comparisons-corrected p-values
for (i in seq_along(region_types)) {
  inter <- region_types[[i]]
  corrected_p_value <- corrected_p_values[i]
  
  if (corrected_p_value < 0.05) {
    cat("Significant region_type difference found for:", paste(inter, collapse=", "), "\n")
    cat("Bonferroni corrected p-value:", corrected_p_value, "\n")
    print(results[[i]])
  }
}




# ---------------------- Change in left PMD - right Cerebellar connectivity depends on vision ---------[specific]-----------------
# Region of interest
region <- "pmd-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)
print(aov)
# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)


# Mean and sd: per session
data_region_nonnorm_pre <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-pre'),]
data_region_nonnorm_post <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-post'),]
cat(sprintf("pre: %.2f ± %.2f\n", mean(data_region_nonnorm_pre$connectivity), sd(data_region_nonnorm_pre$connectivity)))
cat(sprintf("post: %.2f ± %.2f\n", mean(data_region_nonnorm_post$connectivity), sd(data_region_nonnorm_post$connectivity)))

# ---- Post hoc tests : Test for vision differences in connectivity change
vision <- c("True", "False")
p_values <- c()
results <- list()
for (i in seq_along(vision)) {
  vis <- vision[[i]]
  # Subset the data
  x <- data_region[data_region$vision == vis & data_region$session == 'ses-pre' & data$subject!= "sub-57", ]$connectivity
  y <- data_region[data_region$vision == vis & data_region$session == 'ses-post' & data$subject!= "sub-57", ]$connectivity
  # Perform the t-test
  result <- t.test(x, y, paired = TRUE)
  results[[i]] <- result
  # Add the result to the list
  p_values <- c(p_values, result$p.value)
}
# Correct the post hoc tests for multiple comparisons
corrected_p_values <- p.adjust(p_values, method="BH")
# Print post hoc tests with multiple-comparisons-corrected p-values
for (i in seq_along(vision)) {
  vis <- vision[[i]]
  corrected_p_value <- corrected_p_values[i]
  
  if (corrected_p_value < 0.05) {
    cat("Significant connectivity change found for vision:", paste(vis, collapse=", "), "\n")
    cat("Bonferroni corrected p-value:", corrected_p_value, "\n")
    print(results[[i]])
  }
}
# ---------------------- Change in right PMD - right Cerebellar connectivity depends on vision ---------[specific]-----------------

# Region of interest
region <- "pmd-right_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)
print(aov)
# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)



# Mean and sd: groupwise
data_region_nonnorm_p <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$group == 'p'),]
data_region_nonnorm_c <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$group == 'c'),]
cat(sprintf("patient: %.2f ± %.2f\n", mean(data_region_nonnorm_p$connectivity), sd(data_region_nonnorm_p$connectivity)))
cat(sprintf("control: %.2f ± %.2f\n", mean(data_region_nonnorm_c$connectivity), sd(data_region_nonnorm_c$connectivity)))


# ---------------------- Change in right PMD - left Cerebellar connectivity depends on vision ---------[specific]-----------------
# Region of interest
region <- "pmd-right_cereb-M3-left"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)
print(aov)
# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)


# Mean and sd: groupwise
data_region_nonnorm_p <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$group == 'p'),]
data_region_nonnorm_c <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$group == 'c'),]
cat(sprintf("patient: %.2f ± %.2f\n", mean(data_region_nonnorm_p$connectivity), sd(data_region_nonnorm_p$connectivity)))
cat(sprintf("control: %.2f ± %.2f\n", mean(data_region_nonnorm_c$connectivity), sd(data_region_nonnorm_c$connectivity)))

# Mean and sd: per session
data_region_nonnorm_pre <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-pre'),]
data_region_nonnorm_post <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-post'),]
cat(sprintf("pre: %.2f ± %.2f\n", mean(data_region_nonnorm_pre$connectivity), sd(data_region_nonnorm_pre$connectivity)))
cat(sprintf("post: %.2f ± %.2f\n", mean(data_region_nonnorm_post$connectivity), sd(data_region_nonnorm_post$connectivity)))

# model prediction plot
plot_model(modelRIM.Interaction, type = "pred", terms = c("session", "group", "vision"))


# ---- Post hoc tests : Test for vision differences in ipsilateral PMd connectivity change depending on group

vision <- c("True", "False")
groups <- c("c", "p")
p_values <- c()
results <- list()
n <- 0
for (i in seq_along(vision)) {
  for (g in seq_along(groups)) {
    vis <- vision[[i]]
    group <- groups[[g]]
    n=n+1
    print(n)
    # Subset the data
    x <- data_region[data_region$vision == vis & data_region$session == 'ses-pre' & data_region$group == group & data$subject!= "sub-57", ]$connectivity
    y <- data_region[data_region$vision == vis & data_region$session == 'ses-post' & data_region$group == group & data$subject!= "sub-57", ]$connectivity
    # Perform the t-test
    result <- t.test(x, y, paired = TRUE)
    results[[n]] <- result
    # Add the result to the list
    p_values <- c(p_values, result$p.value)
    
  }
}
# Correct the post hoc tests for multiple comparisons
corrected_p_values <- p.adjust(p_values, method="BH")

# Print post hoc tests with multiple-comparisons-corrected p-values
n <- 0
for (i in seq_along(vision)) {
  for (g in seq_along(groups)) {
    
    vis <- vision[[i]]
    group <- groups[[g]]
    n=n+1
    print(n)
    corrected_p_value <- corrected_p_values[n]
    
    if (corrected_p_value < 0.05) {
      cat("Significant connectivity change found for:", paste(vis, group, collapse=", "), "\n")
      cat("Bonferroni corrected p-value:", corrected_p_value, "\n")
      result <- results[[n]]
      print(result)
    }
  }
}


# Control analyses PMd 
# Region of interest
region <- "pmd-left_pmd-right"
# region <- "cereb-M3-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)

aov <- anova(modelRIM.Interaction)
print(aov)

feedback <- c("True", "False")
vision <- c("True", "False")
groups <- c("c", "p")
p_values <- c()
results <- list()
n <- 0
for (f in seq_along(feedback)) {
  for (i in seq_along(vision)) {
    for (g in seq_along(groups)) {
      fed <- feedback[[f]]
      vis <- vision[[i]]
      group <- groups[[g]]
      n=n+1
      print(n)
      # Subset the data
      x <- data_region[data_region$vision == vis & data_region$feedback == fed & data_region$session == 'ses-pre' & data_region$group == group & data$subject!= "sub-57", ]$connectivity
      y <- data_region[data_region$vision == vis & data_region$feedback == fed & data_region$session == 'ses-post' & data_region$group == group & data$subject!= "sub-57", ]$connectivity
      # Perform the t-test
      result <- t.test(x, y, paired = TRUE)
      results[[n]] <- result
      # Add the result to the list
      p_values <- c(p_values, result$p.value)
      
    }
  }
}
# Correct the post hoc tests for multiple comparisons
corrected_p_values <- p.adjust(p_values, method="BH")


# ---------------------- Connectivity change in left PPC - right cerebellum motor depends on vision --------------[specific]-------
region <- "ppc-left_cereb-M3-right"
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)
# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)

# model prediction plot
plot_model(modelRIM.Interaction, type = "pred", terms = c("session", "vision"))

# Mean and sd: groupwise
data_region_nonnorm_p <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$group == 'p'),]
data_region_nonnorm_c <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$group == 'c'),]
cat(sprintf("patient: %.2f ± %.2f\n", mean(data_region_nonnorm_p$connectivity), sd(data_region_nonnorm_p$connectivity)))
cat(sprintf("control: %.2f ± %.2f\n", mean(data_region_nonnorm_c$connectivity), sd(data_region_nonnorm_c$connectivity)))

# Mean and sd: per session
data_region_nonnorm_pre <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-pre'),]
data_region_nonnorm_post <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-post'),]
cat(sprintf("pre: %.2f ± %.2f\n", mean(data_region_nonnorm_pre$connectivity), sd(data_region_nonnorm_pre$connectivity)))
cat(sprintf("post: %.2f ± %.2f\n", mean(data_region_nonnorm_post$connectivity), sd(data_region_nonnorm_post$connectivity)))

# ---- Post hoc tests: Test for vision differences in contralateral PPC connectivity
# Test for connectivity change from pre to post depending on vision in all participants
vision <- c("True", "False")
results_list <- list()
p_values <- c()
for (vis in vision) {
    # Subset the data
    x <- data[data$regions == region & data$vision == vis & data$session == "ses-pre" & data$subject!= "sub-57", ]$connectivity
    y <- data[data$regions == region & data$vision == vis & data$session == "ses-post" & data$subject!= "sub-57", ]$connectivity
    # Perform the t-test
    var.test(x,y)
    result <- t.test(x, y, paired = TRUE, var.equal = TRUE)
    # Add the result to the list
    results[[vis]] <- result
    # Add the result to the list
    p_values <- c(p_values, result$p.value)
    print(paste(vis, sep=","))
}
# View the results corrected for multiple comparisons
p.adjust(p_values, method="BH")
# --> In the conditions with vision, connectivity changes significantly from pre to post


# ---------------------- Connectivity change in right PPC - left cerebellum motor ---------------------
region <- "ppc-right_cereb-M3-left"
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)
# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)

# ---- Post hoc tests: Test for vision differences in ipsilateral PPC connectivity
vision <- c("True", "False")
results_list <- list()
p_values <- c()
for (vis in vision) {
  # Subset the data
  x <- data[data$regions == region & data$vision == vis & data$session == "ses-pre" & data$subject!= "sub-57", ]$connectivity
  y <- data[data$regions == region & data$vision == vis & data$session == "ses-post" & data$subject!= "sub-57", ]$connectivity
  # Perform the t-test
  var.test(x,y)
  result <- t.test(x, y, paired = TRUE, var.equal = TRUE)
  # Add the result to the list
  results[[vis]] <- result
  # Add the result to the list
  p_values <- c(p_values, result$p.value)
  print(paste(vis, sep=","))
}
# View the results corrected for multiple comparisons
p.adjust(p_values, method="BH")
# --> Not significant in post hoc tests


# Control regions
region <- "ppc-left_ppc-right"
# region <- "cereb-M3-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)

aov <- anova(modelRIM.Interaction)
print(aov)


# ---------------------- Connectivity increases beween M1 & Cerebellum (hand area) in all participants - expectedly contralateral to the moving hand --------------------------
# Regions of interest
region <- "m1-hand-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)
# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)

# model prediction plot
plot_model(modelRIM.Interaction, type = "pred", terms = c("session", "group"))

# Mean and sd: per session
data_region_nonnorm_pre <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-pre'),]
data_region_nonnorm_post <- data_non_norm[(data_non_norm$regions == region) & (data_non_norm$session == 'ses-post'),]
cat(sprintf("pre: %.2f ± %.2f\n", mean(data_region_nonnorm_pre$connectivity), sd(data_region_nonnorm_pre$connectivity)))
cat(sprintf("post: %.2f ± %.2f\n", mean(data_region_nonnorm_post$connectivity), sd(data_region_nonnorm_post$connectivity)))

# Control analysis: Ipsilateral M1
region <- "m1-hand-right_cereb-M3-left"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)

# Control analysis: M1-M1
region <- "m1-hand-left_m1-hand-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

# Control analysis: Cereb-Cereb
region <- "cereb-M3-left_cereb-M3-right"
data_region <- data[data$regions == region,]
modelRIM.Interaction      = lmer(connectivity ~  group*session*vision*feedback + (1 | subject), data=data_region, REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)
# ---------------------- Connectivity change in right dlPFC - right D2 depends on feedback and group -----------[specific]-------
region <- "dlpfc-right_cereb-D2-right"
modelRIM.Interaction      = lmer(connectivity ~  group*session*feedback*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)
# Save table as latex table
split_regions <- strsplit(region, "_")
latex_anova <- xtable(aov, caption = paste("Linear mixed effects model for connectivity between ", split_regions[[1]][1], " and ", split_regions[[1]][2], sep=""))
file_path <- paste(table_path, region, ".tex", sep="")
print(latex_anova, include.rownames = TRUE, file = file_path)

# model prediction plot
plot_model(modelRIM.Interaction, type = "pred", terms = c("session", "group","feedback"))


# Control analysis:
region <- "dlpfc-right_cereb-D2-left"
modelRIM.Interaction      = lmer(connectivity ~  group*session*feedback*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)

# Control analysis:
region <- "cereb-D2-left_cereb-D2-right"
modelRIM.Interaction      = lmer(connectivity ~  group*session*feedback*vision + (1 | subject), data=data[data$regions == region,], REML = FALSE)
aov <- anova(modelRIM.Interaction)  
print(aov)
# ---- Post hoc tests : Test for feedback effect on dlpfc connectivity change depending on group

feedback <- c("True", "False")
groups <- c("c", "p")
p_values <- c()
results <- list()
n <- 0
for (i in seq_along(feedback)) {
  for (g in seq_along(groups)) {
    fed <- feedback[[i]]
    group <- groups[[g]]
    n=n+1
    # Subset the data
    x <- data_region[data_region$feedback == fed & data_region$session == 'ses-pre' & data_region$group == group & data$subject!= "sub-57", ]$connectivity
    y <- data_region[data_region$feedback == fed & data_region$session == 'ses-post' & data_region$group == group & data$subject!= "sub-57", ]$connectivity
    # Perform the t-test
    result <- t.test(x, y, paired = TRUE)
    results[[n]] <- result
    # Add the result to the list
    p_values <- c(p_values, result$p.value)
    
    # Print result
    print(fed)
    print(group)
    print(result)
    
    
  }
}
# Correct the post hoc tests for multiple comparisons
corrected_p_values <- p.adjust(p_values, method="BH")

# Print post hoc tests with multiple-comparisons-corrected p-values
n <- 0
for (i in seq_along(feedback)) {
  for (g in seq_along(groups)) {
    
    fed <- feedback[[i]]
    group <- groups[[g]]
    n=n+1
    print(n)
    corrected_p_value <- corrected_p_values[n]
    
    if (corrected_p_value < 0.05) {
      cat("Significant connectivity change found for:", paste(fed, group, collapse=", "), "\n")
      cat("Bonferroni corrected p-value:", corrected_p_value, "\n")
      result <- results[[n]]
      print(result)
    }
  }
}

