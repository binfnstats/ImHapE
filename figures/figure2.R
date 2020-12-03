# SARS-CoV-2 Figures
# ------------------
# -----
# Figure 2
# -----

DATA_DIR = 'DIRECTORY_STORING_DATA'
SLIDE_DIR = 'DIRECTORY_STORING_SLIDING_WINDOW_ANALYSIS_OF_EMPIRICAL_SARS-COV-2_GENOMES'
VAF_FILE = 'PATH_TO_VAF_FILE'
SAVE_PLOT_PDF = 'figure2.pdf'
SAVE_PLOT_PNG = 'figure2.png'

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggsci)
library(ggpubr)
library(EnvStats)
library(ggrepel)

# Figure: COVID sample overview
# =============================
setwd(DATA_DIR)

sample_ids <- na.omit(merge(fread('covid_ids_20268_18877.tsv'), fread('covid_hamming.csv'), by = 'filename'))
sample_ids <- sample_ids  %>% filter(!month %in% c(1, 2, 8, 12))
sample_ids$month[sample_ids$month == 3] <- 'March 2020'
sample_ids$month[sample_ids$month == 4] <- 'April 2020'
sample_ids$month[sample_ids$month == 5] <- 'May 2020'
sample_ids$month[sample_ids$month == 6] <- 'June 2020'
sample_ids$month[sample_ids$month == 7] <- 'July 2020'
sample_ids$month <- factor(sample_ids$month, levels = c('March 2020', 'April 2020', 'May 2020', 'June 2020', 'July 2020'))
sample_ids <- sample_ids %>% select(region, month, hamming)

# Simulation counts
sim_counts <- fread('simulation_snp_counts.tsv') %>% filter(!sample %like% 'RNN') # Remove redundant sims
sim_counts$region <- 'Simulated Genomes'
sim_counts <- sim_counts %>% select(region, persample)
colnames(sim_counts)[2] <- 'hamming'

A.samples <- ggplot(sample_ids, aes(x = fct_rev(region))) + 
  geom_bar(aes(fill = region), color = 'black', show.legend = FALSE) + 
  theme(panel.border = element_rect(color = "black", fill = rgb(0,0,0,0)),
        panel.background = element_rect(fill = rgb(0,0,0,0)),
        strip.background = element_blank(),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.text = element_text(size = 8.5, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8.5),
        axis.text.x = element_text(angle = 45, hjust = 0.98),
        plot.title = element_text(size = 8),
        legend.position = 'top',
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  scale_fill_npg() +
  labs(y = 'Number of samples', x = 'Global region\n') +
  facet_wrap(~as.factor(month), ncol = 1) +
  scale_y_log10()

B.hamming <- ggplot(sample_ids %>% filter(grepl('Europe|North America', region)), aes(x = fct_rev(month), y = hamming)) + 
  geom_boxplot(aes(fill = region), show.legend = FALSE) +
  #scale_y_log10() +
  facet_wrap(~region, ncol = 1) +
  scale_fill_manual(values = pal_npg()(8)[c(2,5,6,8)]) +
  theme(panel.border = element_rect(color = "black", fill = rgb(0,0,0,0)),
        panel.background = element_rect(fill = rgb(0,0,0,0)),
        strip.background = element_blank(),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.text = element_text(size = 8.5, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8.5),
        axis.title.y = element_text(color = 'white'),
        plot.title = element_text(size = 8),
        legend.position = 'top',
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  coord_flip() +
  labs(y = 'Number of mutations per sample', x = 'Month') +
  stat_n_text(size = 2.5, hjust = 0.3)

C.sims <- ggplot(sim_counts[1:1000,], aes(x = region, y = hamming)) + 
  geom_boxplot(aes(fill = region), outlier.shape = NA, show.legend = FALSE) +
  #scale_y_log10() +
  facet_wrap(~region, ncol = 1) +
  scale_fill_manual(values = pal_npg()(8)[c(2,5,6,8)]) +
  theme(panel.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 8.5, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8.5),
        axis.ticks.y = element_line(color = 'white'),
        axis.text.y = element_text(color = 'white', size = 3.8),
        axis.title.y = element_text(color = 'white'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm')) +
  coord_flip() +
  labs(y = 'Number of mutations per sample', x = 'Month') +
  stat_n_text(size = 2.5, hjust = 1.25) +
  ylim(c(-80,896))

# Width: 7.17, Height: 5.53
out <- ggarrange(A.samples, ggarrange(B.hamming, C.sims, align = 'v', ncol = 1, heights = c(4, 1)), ncol = 2, widths = c(1,2.3), labels = c('a', 'b'))  
#ggsave(out, file = '../figures/fig1_sample_info.pdf', width = 6.42, height = 7.08, units = 'in')

# Figure: COVID sample overview
# =============================
setwd(SLIDE_DIR)

slide <- data.table()
for (i in list.files()) {
  add <- fread(i)
  add$month <- strsplit(strsplit(i, '[.]')[[1]][1], '[_]')[[1]][5]
  add$month[add$month == 3] <- 'March 2020'
  add$month[add$month == 4] <- 'April 2020'
  add$month[add$month == 5] <- 'May 2020'
  add$month[add$month == 6] <- 'June 2020'
  add$month[add$month == 7] <- 'July 2020'
  add$month[add$month == 0] <- 'All months combined'
  add$region <- strsplit(strsplit(i, '[.]')[[1]][1], '[_]')[[1]][6]
  add$region[add$region == 'EastAsiaPacific'] <- 'East Asia & Pacific'
  add$region[add$region == 'NorthAmerica'] <- 'North America'
  add$region[add$region == 'EuropeCentralAsia'] <- 'Europe & Central Asia'
  slide <- rbind(slide, add)
}

# Get VAF data
# ------------
vafs <- data.table()
for (i in list.files(VAF_FILE)) {
  vaf <- fread(paste(VAF_FILE, i, sep=''))
  colnames(vaf) <- 'VAF'
  vaf$region <- strsplit(strsplit(i, '_', 3)[[1]][3], '.', 2)[[1]][1]
  vaf$index <- seq(1, nrow(vaf), 1)
  vafs <- rbind(vafs, vaf)
}
vafs$region[vafs$region == 'EastAsiaPacific'] <- 'East Asia & Pacific'
vafs$region[vafs$region == 'NorthAmerica'] <- 'North America'
vafs$region[vafs$region == 'EuropeCentralAsia'] <- 'Europe & Central Asia'
vafs <- vafs %>% filter(grepl('Europe|North', region))

# Get genome features
# -------------------
genome_features <- fread("../covid_gene_annotation.tsv")  %>% filter(Category %in% c('proteins')) %>% select(Start, End, Label, Category)
#genome_features <- genome_features[sample(1:nrow(genome_features), size = nrow(genome_features), replace = FALSE),]
genome_features <- distinct(genome_features, Label, .keep_all = TRUE)
genome_features$id <- rev(seq(1, nrow(genome_features)))
gf_unmelt <- genome_features
genome_features <- melt(genome_features, id.vars=c("Label", 'Category', 'id'))
genome_features <- genome_features[order(genome_features$Label),]

plotTemporal <- function(data, latest_month, vaf, location) {
  data$month <- factor(data$month, levels = c('March 2020', 'April 2020', 'May 2020', 'June 2020', 'July 2020'))
  pltA <- ggplot(data %>% filter(month != latest_month & region == location), aes(x = (index1+index2)/2, y = rnn_mean)) + 
    geom_ribbon(aes(ymin = rnn_lower, ymax = rnn_upper), fill = "grey50", alpha = 0.8) +
    geom_line(color = pal_npg()(4)[4], show.legend = FALSE, alpha = 1) +
    geom_point(aes(color = ifelse(rnn_lower > 0.5, 'black', rgb(0,0,0,0)), fill = ifelse(rnn_lower > 0.5, 'black', rgb(0,0,0,0))),shape = 21, size = 1, alpha = 0.3, show.legend = FALSE) +
    scale_x_continuous(breaks = c(seq(0, 24000, 4000), 29903)) +#expand = expansion(mult = c(0.08, 0.05))) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = function(x) sprintf("%.2f", x)) +
    theme(panel.border = element_rect(color = "black", fill = rgb(0,0,0,0)),
          panel.background = element_rect(fill = rgb(0,0,0,0)),
          strip.background = element_blank(),
          panel.grid = element_line(color = rgb(0,0,0,0)),
          strip.text = element_text(size = 8.5, color = "black"),
          axis.text = element_text(size = 7, color = "black"),
          axis.title = element_text(size = 8.5, hjust = 0),
          #axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 8),
          legend.key.size = unit(0.8, 'lines'),
          legend.position = 'top',
          legend.key = element_rect(fill = 'white'),
          legend.text = element_text(size = 7.5),
          legend.title = element_text(size = 8)) +
    labs(y = 'Probability Positive Selection (CNN)', x = 'Genomic Position (bp)') +
    geom_hline(yintercept = 0.5, linetype = 'dotted') +
    scale_color_manual(values = c(NA, pal_npg()(4))) +
    scale_fill_manual(values = c(NA, rgb(0,0,0,0.1))) +
    facet_wrap(~month, ncol = 1)
  
  pltB <- ggplot(data %>% filter(region == location & month == 'July 2020'), aes(x = (index1+index2)/2, y = rnn_mean)) + 
    geom_ribbon(aes(ymin = rnn_lower, ymax = rnn_upper), fill = "grey50", alpha = 0.8) +
    geom_line(color = pal_npg()(4)[4], show.legend = FALSE, alpha = 0.8) +
    geom_col(data = vaf %>% filter(region == location), aes(x = index, y = VAF), size = 1, color = pal_npg()(5)[2], fill = pal_npg()(5)[2], alpha = 0.7) +
    geom_point(aes(color = ifelse(rnn_lower > 0.5, 'black', rgb(0,0,0,0)), fill = ifelse(rnn_lower > 0.5, 'black', rgb(0,0,0,0))),shape = 21, size = 1, alpha = 0.3, show.legend = FALSE) +
    scale_x_continuous(breaks = c(seq(0, 24000, 4000), 29903), expand = c(0.001,0.001)) +#expand = expansion(mult = c(0.08, 0.05))) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.5, 1)) +
    theme(panel.border = element_rect(color = "black", fill = rgb(0,0,0,0)),
          panel.background = element_rect(fill = rgb(0,0,0,0)),
          strip.background = element_blank(),
          panel.grid = element_line(color = rgb(0,0,0,0)),
          strip.text = element_text(size = 8.5, color = "black"),
          axis.text = element_text(size = 7, color = "black"),
          axis.title = element_text(size = 8.5),
          axis.title.y = element_text(color = 'white'),
          plot.title = element_text(size = 8),
          legend.key.size = unit(0.8, 'lines'),
          legend.position = 'top',
          legend.key = element_rect(fill = 'white'),
          legend.text = element_text(size = 7.5),
          legend.title = element_text(size = 8)) +
    labs(y = 'Probability Positive Selection (CNN)', x = 'Genomic Position (bp)') +
    geom_hline(yintercept = 0.5, linetype = 'dotted') +
    scale_color_manual(values = c(NA, pal_npg()(4))) +
    scale_fill_manual(values = c(NA, rgb(0,0,0,0.1))) +
    geom_line(data = genome_features, aes(x = value, y = -0.3), color = rep(c(rev(pal_npg()(10)), 'grey'), 2), show.legend = FALSE, size = 5) +
    geom_text_repel(data = aggregate(genome_features$value, by = list(Label = genome_features$Label), FUN = mean), aes(x = x, y = -0.3, label=Label), segment.alpha = 0.5, min.segment.length = 0, force = 15, size=2.5) +
    facet_wrap(~month, ncol = 1)
  plt <- ggarrange(pltA, pltB, ncol = 1, heights = c(1.5, 1)) + theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
  
  return(plt)

}

NAR <- plotTemporal(data = slide %>% filter(month != 'All months combined'), vaf = vafs, latest_month = 'July 2020', location = 'North America')
ECAR <- plotTemporal(data = slide %>% filter(month != 'All months combined'), vaf = vafs, latest_month = 'July 2020', location = 'Europe & Central Asia')

out_plot <- ggarrange(out, ggarrange(NAR, ECAR, ncol = 1, labels = c('c', 'd')), ncol = 2)

ggsave(out_plot, filename = SAVE_PLOT_PDF, units = 'in', width = 12.61, height = 8.78)
ggsave(out_plot, filename = SAVE_PLOT_PNG, units = 'in', width = 12.61, height = 8.78)
