# SARS-CoV-2 Figures
# ------------------
# -----
# Extended Data Figure 2: Relationship between CNN and RNN estimates and VAF
# -----

DATA_DIR = 'SET_WD_TO_DIRECTORY_STORING_DATA'
SIM_DIR = 'SET_WD_TO_DIRECTORY_STORING_SIMULATED_GENOME_DATA'

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggsci)
library(ggpubr)
library(EnvStats)
library(ggrepel)

# Read in data
# =============================
setwd(SIM_DIR)
dt.all <- data.table()
ids <- unique(sapply(strsplit(sapply(strsplit(list.files(), '_', 4), '[', 6), '.', 2), '[', 1))
for (i in ids) {
  dt.pred <- fread(list.files()[list.files() %like% i & list.files() %like% 'predictions'])
  dt.loci <- fread(list.files()[list.files() %like% i & list.files() %like% 'loci']) %>% filter(VAF >= 0.01)
  
  # Remove neutral sites within beneficial mutation windows
  dt.remove <- dt.loci %>% filter(mode == 'Beneficial') %>% select(loci, mode) %>% `colnames<-`(c('ben.loci', 'mode'))
  dt.loci$keep <- 0
  for (i in 1:length(dt.loci$loci)) {
    keep <- 0
    for (j in dt.remove$ben.loci) {
      if (dt.loci$mode[i] == 'Beneficial') {
        keep <- 1
      } else if (abs(dt.loci$loci[i] - j) < 4000) {
        keep <- 0
        break
      } else {
        keep <- 1
      }
    }
    dt.loci$keep[i] <- keep
  }
  
  # Add values to loci table
  dt.loci <- dt.loci %>% filter(keep == 1)
  dt.loci$mean <- sapply(dt.loci$loci, FUN = function(x) max(dt.pred[x > dt.pred$index1 & x < dt.pred$index2,]$mean))
  dt.loci$lower <- sapply(dt.loci$loci, FUN = function(x) max(dt.pred[x > dt.pred$index1 & x < dt.pred$index2,]$lower))
  dt.loci$upper <- sapply(dt.loci$loci, FUN = function(x) max(dt.pred[x > dt.pred$index1 & x < dt.pred$index2,]$upper))
  dt.loci$RNN <- sapply(dt.loci$loci, FUN = function(x) max(dt.pred[x > dt.pred$index1 & x < dt.pred$index2,]$RNN))
  #dt.both <- merge(dt.pred, dt.loci, by = 'index1') %>% filter(VAF >= 0.01)
  #dt.both$fitness <- strsplit(list.files()[list.files() %like% i][1], '_')[[1]][2]
  dt.all <- rbind(dt.all, dt.loci)
}

dt.sample <- dt.all[c(sample(which(dt.all$mode == 'Beneficial'), 900), sample(which(dt.all$mode == 'Neutral'), 900)),]
A <- ggplot(dt.sample, aes(x = mode, y = mean - RNN)) + 
  geom_violin(draw_quantiles = 0.5, aes(fill = mode), show.legend = FALSE) + 
  geom_jitter(alpha = 0.5, size = 1, width = 0.1) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  labs(x = 'Mutation type', y = 'CNN estimate - RNN estimate') +
  scale_fill_npg()

B <- ggplot(dt.sample, aes(x = VAF)) + 
  geom_histogram(aes(fill = mode), show.legend = FALSE) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  labs(x = 'Allele Frequency (900 loci per mutation type from\n 600 genomes from independent population simulations)', y = 'Count') +
  scale_fill_npg() +
  facet_wrap(~mode)

C <- ggplot(dt.all, aes(x = VAF, y = mean - RNN)) + 
  geom_point(aes(color = mode), alpha = 0.5, size = 1) + 
  geom_smooth(aes(color = mode), method = 'lm') + 
  scale_x_log10() +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5),
        legend.position = c(1,1),
        legend.background = element_rect(fill = rgb(0,0,0,0)),
        legend.key = element_rect(fill = rgb(0,0,0,0)),
        legend.justification = c(1.5,1.5),
        legend.key.size = unit(0.4, 'cm')) +
  labs(x = 'Allele Frequency', y = 'CNN estimate - RNN estimate', color = 'Mutation type') +
  scale_color_npg()

abplot <- ggarrange(ggarrange(B, C, nrow = 1, labels = c('a', 'b')), A, ncol = 1, labels = c(NA, 'c'))

ggsave(abplot, filename = '~/Research/sars-cov-2/figures/extended_data_figure2.png', units = 'in', width = 8.73, height = 8.78)