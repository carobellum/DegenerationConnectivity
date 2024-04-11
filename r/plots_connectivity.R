# +++ Raincloud plots +++

# install.packages("readr")
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
# library(reshape2)

source("~/code/R/geom_flat_violin.R")

# Define raincloud theme
raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

# get data
source("~/Documents/Projects/Degeneration/code/cerebellar_degeneration/r/data_connectivity.R")

extension <- ".png"

# ---------------------- Cortico-cerebellar connectivity differs more than cortico-cortico connectivity in patients --------------------------


custom_colours <- c("#153969", "#BB271A")

#  ----  Plot baseline group difference ordered by region_typeity with jittered points  ---- 

ggplot(data = base_connectivity_non_norm, aes(y = connectivity, x = region_type, fill = group, colour = group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, width=0.5) +
  geom_point(aes(y = connectivity, x = region_type, color = group), position = position_jitter(width = .15), size = 1, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 2) +
  theme_bw() +
  xlab("") +
  ylab("Connectivity") +
  scale_colour_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  raincloud_theme

# Save plot
ggsave(paste(savepath, "raincloud_connectivity_region_typexgroup_all",extension, sep=""), width = 5, height = 5)


#  ----  Plot baseline group difference ordered by region_typeity without jittered points  ---- 
ggplot(data = base_connectivity_non_norm, aes(y = connectivity, x = region_type, fill = group, colour = group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, width=0.5) +
  geom_boxplot(width = .3, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 2) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  theme_bw() +
  xlab("") +
  ylab("Connectivity") +
  scale_colour_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  raincloud_theme

# Save plot
ggsave(paste(savepath, "raincloud_connectivity_region_typexgroup_all_nopoints", extension, sep=""), width = 5, height = 5)




# ---- Plot baseline group difference with points ---- 
ggplot(data = base_connectivity_non_norm, aes(y = connectivity, x = group, fill = group, colour = group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, width=0.5) +
  geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 2) +
  geom_point(aes(y = connectivity, x = group, color = group), position = position_jitter(width = .05), size = 1, alpha = 0.8) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  theme_bw() +
  scale_colour_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  xlab("") +
  ylab("Connectivity") +
  raincloud_theme

# Save plot
ggsave(paste(savepath, "raincloud_connectivity_group_all", extension, sep=""), width = 5, height = 5)
  



