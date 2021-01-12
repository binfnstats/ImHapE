#!/bin/env Rscript

# ===============
# Figure 2 panels
# ===============

SLIDING_DIRECTORY = '~/Research/sars-cov-2/data/sliding_output/'
FIGURE_OUTPUT = '~/Research/sars-cov-2/figures/'

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)

setwd(SLIDING_DIRECTORY)

# GISAID data
# -----------
gisaid.counts <- count(fread('../aligned_meta.csv') %>% filter(region %in% c('Asia', 'Europe', 'North America') & prop_masked <= 0.001), month, region) %>% `colnames<-`(c('month', 'region', 'n'))
gisaid <- data.table()
for (g in list.files()[list.files() %like% 'GISAID']) {
  add <- fread(g)
  add$file <- g
  gisaid <- rbind(gisaid, add)
}

gisaid$month <- as.numeric(sapply(strsplit(gisaid$file, '_'), '[', 8))
gisaid$region <- sapply(strsplit(gisaid$file, '_'), '[', 7)
gisaid[gisaid$region == 'NorthAmerica',]$region <- 'North America'

# Remove Asia Month 7, not enough samples
gisaid.agg <- data.table(aggregate(gisaid$rnn_regression, by = list(gisaid$month, gisaid$region), FUN = function(x) c(mean(x), mean(x) - qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x)), mean(x) + qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))))) %>% `colnames<-`(c('month', 'region', 'mean', 'lower', 'upper'))
A1_ <- ggplot(gisaid %>% filter(!paste0(region, month) %in% paste0('Asia', 7)), aes(x = as.numeric(month), y = rnn_regression)) + geom_point(aes(color = region), show.legend = FALSE) +
  geom_errorbar(data = gisaid.agg %>% filter(!paste0(region, month) %in% paste0('Asia', 7)), aes(ymin = lower, y = mean, ymax = upper)) +
  geom_line(data = gisaid.agg %>% filter(!paste0(region, month) %in% paste0('Asia', 7)), aes(x = as.numeric(month), y = mean)) +
  labs(x = 'Month (2020)', y = 'Fitness (1 + s)') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  facet_wrap(~region, nrow = 1) +
  scale_color_npg() +
  scale_y_continuous(breaks = seq(1, 1.6, 0.1))

gisaid.agg.class <- data.table(aggregate(gisaid$rnn_classification, by = list(gisaid$month, gisaid$region), FUN = function(x) c(mean(x), mean(x) - qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x)), mean(x) + qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))))) %>% `colnames<-`(c('month', 'region', 'mean', 'lower', 'upper'))
A2_ <- ggplot(gisaid %>% filter(!paste0(region, month) %in% paste0('Asia', 7)), aes(x = as.numeric(month), y = rnn_classification)) + geom_point(aes(color = region), show.legend = FALSE) +
  geom_crossbar(data = gisaid.agg.class %>% filter(!paste0(region, month) %in% paste0('Asia', 7)), aes(ymin = lower, y = mean, ymax = upper)) +
  labs(x = 'Month (2020)', y = 'Fitness (1 + s)') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  facet_wrap(~region, nrow = 1) +
  scale_color_npg() +
  labs(y = 'P(Selection)') +
  geom_hline(yintercept = c(0.5), linetype = 'dashed') +
  scale_y_continuous(limits = c(0.35, 1))#+
#scale_y_continuous(breaks = seq(1, 1.6, 0.1))

A3_ <- ggplot(gisaid.counts %>% filter(month %in% c(3,4,5,6,7)), aes(x = as.numeric(month), y = n)) + 
  geom_col(aes(fill = ifelse(n < 200, 'green', 'blue')), show.legend = FALSE) +
  labs(x = 'Month (2020)', y = 'Samples') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm')) +
  scale_fill_manual(values = c('grey50', pal_npg()(1)[1])) +
  facet_wrap(~region, nrow = 1)



# COG UK data
# -----------
cog.counts <- count(fread('../cog_meta_withweeks.csv') %>% filter(prop_masked <= 0.018), epi_week) %>% `colnames<-`(c('week', 'n'))
cog_uk <- data.table()
for (g in list.files()[list.files() %like% 'cog_england']) {
  add <- fread(g)
  add$file <- g
  cog_uk <- rbind(cog_uk, add)
}
cog_uk$week <- as.numeric(sapply(strsplit(cog_uk$file, '_'), '[', 7))
cog_uk <- cog_uk %>% filter(week >= 16) # LONDON reference genome from week 16, only use weeks above 17
cog_uk.agg <- data.table(aggregate(cog_uk$rnn_regression, by = list(cog_uk$week), FUN = function(x) c(mean(x), mean(x) - qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x)), mean(x) + qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))))) %>% `colnames<-`(c('week', 'mean', 'lower', 'upper'))

B1_ <- ggplot(cog_uk %>% filter(week %in% filter(cog.counts, n > 200)$week), aes(x = as.numeric(week), y = rnn_regression)) + 
  geom_point(size = 0.5) +
  geom_errorbar(data = cog_uk.agg %>% filter(week %in% filter(cog.counts, n > 200)$week), aes(ymin = lower, y = mean, ymax = upper)) +
  geom_line(data = cog_uk.agg %>% filter(week %in% filter(cog.counts, n > 200)$week), aes(x = week, y = mean)) +
  labs(x = 'Week (2020)', y = 'Fitness (1 + s)') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm'),
        title = element_text(size = 8)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = 16, linetype = 'dashed') +
  scale_y_continuous(breaks = seq(1, 1.6, 0.1)) +
  scale_x_continuous(breaks = c(16, seq(20, 50, 5)), expand = c(0.01,0.01), limits = c(NA, 51)) +
  ggtitle('COG UK')

cog_uk.agg.class <- data.table(aggregate(cog_uk$rnn_classification, by = list(cog_uk$week), FUN = function(x) c(mean(x), mean(x) - qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x)), mean(x) + qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))))) %>% `colnames<-`(c('week', 'mean', 'lower', 'upper'))
B2_ <- ggplot(cog_uk %>% filter(week %in% filter(cog.counts, n > 200)$week), aes(x = as.numeric(week), y = rnn_classification)) + 
  geom_point(size = 0.5) +
  geom_line(data = cog_uk.agg.class %>% filter(week %in% filter(cog.counts, n > 200)$week), aes(x = week, y = mean)) +
  geom_crossbar(data = cog_uk.agg.class %>% filter(week %in% filter(cog.counts, n > 200)$week), aes(x = week, ymin = lower, ymax = upper, y = mean), fill = pal_npg()(6)[6], alpha = 0.8) +
  labs(x = 'Week (2020)', y = 'P(Selection)') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm')) +
  geom_hline(yintercept = c(0.5,1), linetype = 'dashed') +
  scale_x_continuous(breaks = c(16, seq(20, 50, 5)), expand = c(0.01,0.01), limits = c(NA, 51)) +
  scale_y_continuous(limits = c(0.35, 1))

lineage.agg <- count(cog.meta, epi_week, pangolin_lineage)
lineage.agg <- left_join(lineage.agg, count(cog.meta, epi_week), by = 'epi_week')
lineage.agg$p <- lineage.agg$n.x / lineage.agg$n.y
B3_ <- ggplot(lineage.agg %>% filter(epi_week >= 16), aes(x = epi_week, y = p*100)) + 
  geom_col(aes(fill = pangolin_lineage, color = ifelse(pangolin_lineage == 'B.1.1.7' | pangolin_lineage == 'B.1.177', 'red', 'black')), show.legend = FALSE) +
  scale_fill_viridis_d(option = 3, direction = -1) +
  scale_color_manual(values = c(NA, 'black')) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm')) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(16, seq(20, 50, 5)), expand = c(0.01,0.01), limits = c(NA, 51)) +
  labs(x = 'Week (2020', y = 'Pangolin\nlineage (%)')

B4_ <- ggplot(cog.counts %>% filter(week >= 16), aes(x = as.numeric(week), y = n)) + 
  geom_col(aes(fill = ifelse(n < 200, 'green', 'blue')), show.legend = FALSE) +
  labs(x = 'Week (2020)', y = 'Samples') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm')) +
  scale_x_continuous(breaks = c(16, seq(20, 50, 5)), expand = c(0.01,0.01), limits = c(NA, 51)) +
  scale_fill_manual(values = c('grey50', pal_npg()(1)[1]))



# Build panels
# ------------
A <- ggarrange(A1_, A2_, A3_, align = 'v', heights = c(1.2, 1, 1.1), ncol = 1)

B <- ggarrange(B1_, B2_, B3_, B4_, align = 'v', heights = c(1.2, 1, 1.1, 1.1), ncol = 1)

#ggsave(A, filename = paste0(FIGURE_OUTPUT, 'supplementary_gisaid_empirical.png'), units = 'in', width = 7.25, height = 4.7)
#ggsave(B, filename = paste0(FIGURE_OUTPUT, 'figure2_raw.pdf'), units = 'in', width = 7.25, height = 4.93)

out <- ggarrange(A, B, widths = c(1, 1.5), nrow = 1)
ggsave(out, filename = paste0(FIGURE_OUTPUT, 'figure2_v2_raw.pdf'), units = 'in', width = 9.67, height = 5.24)
