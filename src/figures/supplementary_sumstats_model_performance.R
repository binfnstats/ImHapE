#!/bin/env Rscript

# ==========================================================
# Supplementary: Tajima/Wu distributions & model performance
# ==========================================================

SUMSTATS_DIR = '~/Research/sars-cov-2/data/sumstats'
VALIDATION_DIR = '~/Research/sars-cov-2/data/validation'
SIMVALIDATION_DIR = '~/Research/sars-cov-2/data/simvalidation'
FIGURE_OUTPUT = '~/Research/sars-cov-2/figures/'

library(caret)
library(data.table)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggridges)
library(EnvStats)
library(tidyr)

setwd(SUMSTATS_DIR)
dt <- data.table()
for (file in list.files()) {
  add <- fread(file)
  dt <- rbind(dt, add)
}
dt <- dt %>% filter(!(id %in% dt[dt$value == 1e-4 & dt$sumstat %like% 'Mutation rate',]$id))
dt <- pivot_wider(data = dt, names_from=sumstat, values_from=value)
dt <- pivot_longer(data = dt, cols = c("Tajima's D", "Fay and Wu's H"))
dt$mode <- ifelse(dt$`Fitness (1 + s)` == 1, 'Neutral', 'Positive')
dt$param <- paste(dt$`Mutation rate (per site/gen)`, dt$window, sep = ' | ')

# Plot summary statistics across mutation rates, window sizes, and prob beneficial
# --------------------------------------------------------------------------------
A <- ggplot(dt, aes(x = mode, y = value)) + 
  #geom_density(aes(fill = mode, color = as.factor(`P(Beneficial)`))) +
  #geom_boxplot(aes(fill = as.factor(`P(Beneficial)`), color = as.factor(`P(Beneficial)`)), outlier.size = 0.6) +
  geom_violin(aes(fill = as.factor(`P(Beneficial)`), color = as.factor(`P(Beneficial)`)), scale = 'width') +
  facet_wrap(name ~ param, scale = 'free_y', ncol = 3) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 7, color = 'black'),
        axis.title = element_text(size = 7.5),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 7.5),
        strip.text = element_text(size = 8),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = 'bottom') +
  labs(x = 'Evolutionary Mode', y = 'Statistic Value', color = 'P(Mutation will be beneficial)', fill = 'P(Mutation will be beneficial)') +
  scale_fill_npg() +
  scale_color_npg() +
  guides(colour = guide_legend(nrow = 1))

ggsave(A, filename = paste0(FIGURE_OUTPUT, 'supplementary_TajWu_eval.png'), units = 'in', width = 7.07, height = 6.61)