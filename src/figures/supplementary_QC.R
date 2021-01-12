#!/bin/env Rscript

# Supplementary figure: evaluate empirical 
# ========================================
COG_META = '~/Research/sars-cov-2/data/cog_meta_withweeks.csv'
GISAID_META =  '~/Research/sars-cov-2/data/aligned_meta.csv'
COG_SITES_MASKED = '~/Research/sars-cov-2/data/sites_masked_nonstandard_COG.csv'
GISAID_SITES_MASKED = '~/Research/sars-cov-2/data/sites_masked_nonstandard_GISAID.csv'
COG_MUTATION = '~/Research/sars-cov-2/data/mutation_counts_COG_wuhan.csv'
COG_MUTATION_ENGLAND = '~/Research/sars-cov-2/data/mutation_counts_COG_england.csv'
GISAID_MUTATION =  '~/Research/sars-cov-2/data/mutation_counts_GISAID_wuhan.csv'

FIGURE_OUTPUT = '~/Research/sars-cov-2/figures/'

library(data.table)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggridges)
library(EnvStats)

# Proportion masked per month
# ---------------------------
cog <- fread(COG_META) %>% filter(month != 2)
gisaid <- fread(GISAID_META) %>% filter(!month %in% c(8, 12))

PA <- ggplot(cog, aes(x = prop_masked*100, y = as.factor(month))) +
  geom_density_ridges(alpha = 0.8, fill = pal_npg()(2)[2]) + 
  geom_vline(xintercept = 1.8, color = 'red') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  xlim(c(NA, 3.5)) +
  labs(x = 'Proportion of genome masked (%)', y = 'Month (2020)')

PB <- ggplot(gisaid, aes(x = prop_masked*100, y = as.factor(month))) +
  geom_density_ridges(alpha = 0.8, fill = pal_npg()(2)[2]) + 
  geom_vline(xintercept = 0.1, color = 'red') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  xlim(c(NA, 1.5)) +
  labs(x = 'Proportion of genome masked (%)', y = 'Month (2020)')

# Masking per position
# --------------------
cog_masked <- fread(COG_SITES_MASKED)
cog_masked$p_masked <- cog_masked$n_masked / 130140
gisaid_masked <- fread(GISAID_SITES_MASKED)
gisaid_masked$p_masked <- gisaid_masked$n_masked / 49320

library(ggrepel)
MA <- ggplot(cog_masked, aes(x = position, y = p_masked*100)) + 
  geom_point(size = 0.2) +
  geom_line() +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  ylim(c(0,100)) +
  labs(x = 'Genomic position (bp)', y = 'Proportion of samples masked (%)') +
  annotate(geom = 'text', label = 'COG United Kingdom Sequences (n = 130140)', x = -Inf, y = Inf, hjust = -0.2, vjust = 1.65, size = 3.5)

MB <- ggplot(gisaid_masked, aes(x = position, y = p_masked*100)) + 
  geom_point(size = 0.2) +
  geom_line() +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  ylim(c(0,100)) +
  labs(x = 'Genomic position (bp)', y = 'Proportion of samples masked (%)') +
  annotate(geom = 'text', label = 'Global GISAID Sequences (n = 49320)', x = -Inf, y = Inf, hjust = -0.25, vjust = 1.65, size = 3.5)

# Mutations per sample
# --------------------
cog_mutation <- fread(COG_MUTATION)
cog_mutation$group <- 'COG England Sequences (Mutations called with Wuhan reference)'
cog_mutation_england <- fread(COG_MUTATION_ENGLAND)
cog_mutation_england$group <- 'COG England Sequences (Mutations called with England-LOND-D51C5-2020 reference)'
gisaid_mutation <- fread(GISAID_MUTATION)
gisaid_mutation$group <- 'GISAID Sequences (Mutations called with Wuhan reference)'
all_mutation <- rbind(cog_mutation, cog_mutation_england, gisaid_mutation)

CA <- ggplot(all_mutation, aes(x = mut_counts)) + 
  geom_density(aes(fill = group), alpha = 0.7) +
  scale_fill_npg() +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = c(0, 1), 
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  scale_x_continuous(trans = 'log10', limits = c(1, 175)) +
  ylim(c(0, 5)) +
  labs(x = 'Mutations per virus', y= 'Density', fill = 'Group')

output <- ggarrange(ggarrange(PA, MA, PB, MB, ncol = 2, nrow = 2, labels = c('a', 'b', 'c', 'd')), CA, ncol = 1, heights = c(2, 1), labels = c(NA, 'e'))

# Save output
# -----------
ggsave(output, filename = paste0(FIGURE_OUTPUT, 'supplementary_QC.png'), width = 9.18, height = 9.09, units = 'in')
