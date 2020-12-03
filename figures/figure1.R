# SARS-CoV-2 Figures
# ------------------

# -----
# Figure 1: CNN and RNN performance and an example simulated genome
# -----

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggsci)
library(EnvStats)

setwd('SET_WD_TO_DIRECTORY_STORING_DATA')

# Plot CNN and RNN performance
# -----------------------------
library(caret)

strSample <- function(cnn) {
  inds <- c()
  for (i in unique(cnn$modes)[unique(cnn$modes) != 'neutral']) {
    ind.samp <- sample(which(cnn$modes == i), 20)
    inds <- c(inds, ind.samp)
  }
  inds <- c(inds, which(cnn$modes == 'neutral'))
  return(inds)
}

cnn <- fread('validation_2500_none_2020-10-29.csv')
cnn.samp <- cnn[strSample(cnn)]
cnn.samp$actual <- ifelse(cnn.samp$modes == 'neutral', 'Neutral', 'Positive selection')
cnn.samp$class <- ifelse(cnn.samp$class == FALSE, 'Neutral', 'Positive selection')

cnn.performance <- confusionMatrix(data = as.factor(cnn.samp$actual), reference = as.factor(cnn.samp$class), positive = 'Positive selection')$table
cnn.recall <- data.table(prop.table(cnn.performance, margin = 1))
colnames(cnn.recall)[3] <- 'Recall'
cnn.counts <- data.table(cnn.performance)
cnn.performance <- merge(cnn.recall, cnn.counts, by = c('Prediction','Reference'))
cnn.cm <- ggplot(cnn.performance, aes(x = Prediction, y = Reference, fill = Recall)) + 
  geom_tile(show.legend = FALSE) +
  geom_text(aes(label = N, color = ifelse(Recall > 0.5, 'white', 'black')), size = 2.8, vjust = -1, show.legend = FALSE) +
  geom_text(aes(label = paste0('(', round(Recall*100, 1), '%)'), color = ifelse(Recall > 0.5, 'white', 'black')), size = 2.8, vjust = 1, show.legend = FALSE) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = '#E5EEFD', high = '#364973') +
  scale_color_manual(values = c('black', 'white')) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm'))

rnn <- fread('RNN_performance.csv')
rnn$actual <- ifelse(rnn$actual == 0, 'No beneficial\n mutation in window', 'Beneficial mutation\n in window')
rnn$pred_binary <- ifelse(rnn$pred_binary == FALSE, 'No beneficial\n mutation in window', 'Beneficial mutation\n in window')
rnn.performance <- confusionMatrix(data = as.factor(rnn$actual), reference = as.factor(rnn$pred_binary), positive = 'Beneficial mutation\n in window')$table
rnn.recall <- data.table(prop.table(rnn.performance, margin = 1))
colnames(rnn.recall)[3] <- 'Recall'
rnn.counts <- data.table(rnn.performance)
rnn.performance <- merge(rnn.recall, rnn.counts, by = c('Prediction','Reference'))

rnn.cm <- ggplot(rnn.performance, aes(x = Prediction, y = Reference, fill = Recall)) + 
  geom_tile(show.legend = FALSE) +
  geom_text(aes(label = N, color = ifelse(Recall > 0.5, 'white', 'black')), size = 2.8, vjust = -1, show.legend = FALSE) +
  geom_text(aes(label = paste0('(', round(Recall*100, 1), '%)'), color = ifelse(Recall > 0.5, 'white', 'black')), size = 2.8, vjust = 1, show.legend = FALSE) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = '#E5EEFD', high = '#364973') +
  scale_color_manual(values = c('black', 'white')) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 8.5, color = 'black'),
        axis.title = element_text(size = 9),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0.1)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 6.5),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 7.5),
        legend.key.size = unit(0.4, 'cm'))

ggsave(cnn.cm, file = '../figures/fig1_cnn_cm.pdf', units = 'in', width = 3.08, height = 2.1)
ggsave(rnn.cm, file = '../figures/fig1_rnn_cm.pdf', units = 'in', width = 3.46, height = 2.2)

# Example simulated genomes
# -------------------------
plotSimulatedGenome <- function(id, directory) {
  # Get data
  dir = getwd()
  setwd(directory)
  loci <- fread(list.files()[list.files() %like% id & list.files() %like% 'loci']) %>% filter(VAF > 0.001)
  loci$modes <- factor(loci$modes, levels = c('Beneficial', 'Neutral'))
  prob <- fread(list.files()[list.files() %like% id & list.files() %like% 'predictions'])
  fitness <- strsplit(list.files()[list.files() %like% id][1], '_')[[1]][4]
  window <- strsplit(list.files()[list.files() %like% id][1], '_')[[1]][2]
  prob$sig <- ifelse(prob$lower > 0.5, 'yes', 'no')
  loci$mode <- factor(loci$mode, levels = c('Neutral', 'Beneficial'))
  
  # Plot data
  plt <- ggplot(prob, aes(x = (index1+index2)/2, y = mean)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey50", alpha = 0.7) +
    geom_line(color = pal_npg()(3)[3], size = 0.2, alpha = 1) +
    geom_line(aes(x = (index1+index2)/2, y = RNN), color = pal_npg()(4)[4], size = 1, alpha = 1) +
    geom_point(aes(alpha = sig), size = 0.8, show.legend = FALSE) +
    scale_alpha_manual(values = c(0, 1)) +
    geom_hline(yintercept = 0.5, color = 'black', linetype = 'dashed') +
    theme(panel.border = element_rect(color = "black", fill = rgb(0,0,0,0)),
          panel.background = element_rect(fill = rgb(0,0,0,0)),
          strip.background = element_blank(),
          panel.grid = element_line(color = rgb(0,0,0,0)),
          plot.title = element_text(size = 8),
          legend.key.size = unit(0.6, 'lines'),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          strip.text = element_text(size = 8, color = "black"),
          axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 8.5)) +
    ggtitle(paste('Genome size = 29903 bp, Sliding window =', window, 'bp, Step size = 100 bp, Fitness (1 + s) = 2')) +
    labs(x = 'Genomic position (bp)', y = 'P(Selection)', fill = 'Mutation type', color = 'Mutation type') +
    scale_x_continuous(breaks = c(seq(0, 24000, 6000), 29903)) +
    scale_y_continuous(expand = expansion(c(0.01, 0.01)), limits = c(0, 1), breaks = c(0, 0.250, 0.500, 0.750, 1)) +
    geom_col(data = loci, aes(x = loci, y = VAF, color = mode, fill = mode), position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_fill_manual(values = pal_npg()(2)[c(2,1)]) +
    scale_color_manual(values = pal_npg()(2)[c(2,1)])
  setwd(dir)
  return(plt)
  
}

example <- plotSimulatedGenome(id = '9493491', directory = 'model_simgenome/')
ggsave(example, file = '../figures/fig1_example_genome.pdf', units = 'in', width = 7.52, height = 2.58)