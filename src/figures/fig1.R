#!/bin/env Rscript

# ===============
# Figure 1 panels
# ===============

FIGURE_OUTPUT = '~/Research/sars-cov-2/figures/'

setwd("~/Research/sars-cov-2/data")
set.seed(123367)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(caret)
library(splitstackshape)

# C. RNN classification performance
# ---------------------------------
rnn_class <- fread('rnn_classification_predictions.csv')

cmat <- data.table()
for (i in 1:100) {
  print(i)
  #dt <- stratified(rnn_class, group = 'fitness', size = 100)
  dt.melt <- melt(rnn_class, id.vars = 'fitness')
  for (i in unique(dt.melt$variable)) {
    #for (fit in unique(dt.melt$fitness)[unique(dt.melt$fitness) != 1]) {
    add <- dt.melt %>% filter(variable == i)#%>% filter((variable == i & fitness == fit) | (fitness == 1 & variable == i))
    add$actual <- as.factor(ifelse(add$fitness == 1, 'Neutral', 'Positive'))
    add <- stratified(add, group = 'actual', size = 100)
    add$prediction <- as.factor(ifelse(add$value <= 0.5, 'Neutral', 'Positive'))
    add <- confusionMatrix(data = add$prediction, reference = add$actual, positive = 'Positive')$byClass
    add <- data.table(t(add))
    add$model <- i
    #add$fitness <- fit
    cmat <- rbind(cmat, add)
    #}
  }
}

cmat_filter <- cmat %>% filter(F1 > 1/1.5) # Only keep models that did not overfit
cmat_filter$sorting <- lapply(sapply(strsplit(cmat_filter$model, '_'), '[', -1), FUN = function(x) paste(x, collapse = '_'))
cmat_filter <- filter(cmat_filter, model %in% distinct(cmat_filter[order(cmat_filter$F1, decreasing = TRUE),], sorting, .keep_all = TRUE)$model)

cmat_filter[cmat_filter$model == '0.0_wuH_CNN']$model <- 'CNN + FWuH'
cmat_filter[cmat_filter$model == '0.0_tajD_wuH']$model <- 'TajD + FWuH'
cmat_filter[cmat_filter$model == '0.0_multi']$model <- 'CNN + TajD + FWuH'
cmat_filter[cmat_filter$model == '0.0_CNN']$model <- 'CNN'
cmat_filter[cmat_filter$model == '0.001_tajD_CNN']$model <- 'CNN + TajD'
cmat_filter[cmat_filter$model == '0.0_tajD']$model <- 'TajD'

C_ <- ggplot(cmat_filter, aes(x = reorder(model, F1), y = F1*100)) + 
  geom_boxplot(outlier.shape = NA, fill = pal_npg()(6)[6]) + 
  geom_jitter(size = 0.4, width = 0.25) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 7.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = c(1, 1.05),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, 'cm')) +
  coord_flip() +
  scale_y_continuous(breaks = c(75, 80, 85, 90, 95, 100), limits = c(NA, 101)) +
  labs(x = 'RNN features', y = 'F1 Score (%)')

# D. RNN classification performance
# ---------------------------------
rnn_reg <- fread('rnn_regression_predictions.csv')
rnn_reg <- rnn_reg %>% select(all_of(unique((cmat %>% filter(F1 > 1/1.5))$model)), 'fitness')

# Ensemble the predictions
trainIndex <- createDataPartition(rnn_reg$fitness, p = .2, list = FALSE, times = 1)
train <- rnn_reg[trainIndex,]
test <- rnn_reg[-trainIndex,]
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 1)
fits <- train(fitness ~ ., data = train, method = "gbm", trControl = fitControl, verbose = FALSE)

preds <- predict(fits, newdata = test)
postResample(pred = preds, obs = test$fitness)
test$pred <- preds
test <- stratified(test, group = 'fitness', size = 175)
mpe <- (100/nrow(test)) * sum(((test$fitness - test$pred) / test$fitness))
D_ <- ggplot(test, aes(x = as.factor(fitness), y = pred)) + 
  geom_jitter(width = 0.1, size = 0.4) +
  geom_boxplot(alpha = 0.5, width = 0.5, fill = pal_npg()(6)[6], outlier.shape = NA) +
  geom_smooth(method = 'lm') +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 7.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = c(1, 1.05),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, 'cm')) +
  labs(x = 'Actual Fitness (1 + s)', y = 'Predicted Fitness (1 + s)') +
  scale_y_continuous(breaks = seq(1, 2, 0.1), limits = c(NA, 2)) +
  annotate(geom = 'text', label = paste0('Mean percentage error (MPE) = ',  as.character(round(mpe,2)), '%'), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.8, size = 3.25)

# E. Number of samples
# -----------------
gisaid_samples <- fread('aligned_meta.csv')
gisaid_samples$group <- 'GISAID'
cog_samples <- fread('cog_meta_withweeks.csv') %>% select(-epi_week)
cog_samples$group <- 'COG UK'
meta_samples <- rbind(gisaid_samples, cog_samples)

E_ <- ggplot(meta_samples, aes(x = as.factor(month))) + 
  geom_bar(aes(fill = group), color = 'black', width = 0.7) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 7.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = c(1, 1.04),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, 'cm')) +
  scale_fill_npg() +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = 'Month (2020)', y = 'Number of samples', fill = 'Samples')

# F. Mutation counts
# -------------------
mut_cog <- fread('mutation_counts_COG_england.csv')
mut_cog$group <- 'COG UK        \n(England reference)'
mut_gisaid <- fread('mutation_counts_GISAID_wuhan.csv')
mut_gisaid$group <- 'GISAID        \n(Wuhan reference)'
mut_sim <- fread('mutation_counts_simgenome.csv')
mut_sim$group <- 'Simulated genomes\n(1000 x 8000 sims)'
mut_sim$mutrate <- mut_sim$mutrate/29903
mut_emp <- rbind(mut_cog, mut_gisaid)

F_ <- ggplot(mut_emp, aes(x = group, y = mut_counts)) + 
  geom_violin(fill = pal_npg()(6)[6], draw_quantiles = 0.5) +
  geom_violin(data = mut_sim, aes(x = group, y = mut_counts, fill = as.factor(mutrate)), scale = 'width') +
  coord_flip() +
  scale_y_log10() +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 7.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm')) +
  scale_fill_npg() +
  labs(x = '', y = 'Number of mutations per genome', fill = ' Mutation rate\n(bp/generation)')

out <- ggarrange(ggarrange(C_, D_, nrow = 1), ggarrange(E_, F_, widths = c(1, 1.2), nrow = 1), nrow = 2)
ggsave(out, filename = paste0(FIGURE_OUTPUT, 'fig1_cdef.pdf'), units = 'cm', width = 18.17, height = 11.68)

# SUPPLEMENTARY MODEL SELECTION FIGURES
# =====================================

# Classification model selection
cmat2 <- cmat
cmat2$regularization <- sapply(strsplit(cmat2$model, '[_]'), '[', 1)
cmat2$features <- lapply(sapply(strsplit(cmat2$model, '[_]'), '[', -1), FUN = function(x) paste(x, collapse = ' + '))
cmat2$features <- gsub('tajD', 'TajD', cmat2$features)
cmat2$features <- gsub('wuH', 'FWuH', cmat2$features)
cmat2$features <- gsub('multi', 'CNN + TajD + FWuH', cmat2$features)
cmat2 <- data.table(aggregate(cmat2$F1, by = list(regularization = cmat2$regularization, features = cmat2$features), FUN = mean))
s1 <- ggplot(cmat2, aes(x = regularization, y = features, fill = x)) + 
  geom_tile() + 
  scale_fill_viridis_c(option = 'magma') +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = 'L1/L2 regularization', y = 'RNN features for classification of positive selection', fill = 'F1 Score') +
  geom_text(aes(label = ifelse(x > 1/1.5, round(x, 2), ''))) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm'))

# Regression model selection
mod.select <- data.table()
for (i in colnames(test)[colnames(test) != 'fitness']) {
  mpe <- (100/nrow(test)) * sum(((test$fitness - as.vector(select(test, all_of(i)))) / test$fitness))
  mae <- as.vector(postResample(as.vector(select(test, all_of(i))), test$fitness)[3])
  mod.select <- rbind(mod.select, data.table(features = i, MPE = mpe, MAE = mae))
}
mod.select[mod.select$features == 'pred',]$features <- 'RNN ensemble'
mod.select <- melt(mod.select, id.vars = c('features'))
s2 <- ggplot(mod.select %>% filter(variable == 'MPE'), aes(x = variable, y = reorder(features, value, decreasing = TRUE), fill = value)) + 
  geom_tile() + 
  scale_fill_viridis_c(option = 'magma') +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = 'Performance metric', y = 'RNN features for quantifying fitness', fill = 'Mean % error') +
  geom_text(aes(label = round(value, 2))) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm'))

s3 <- ggplot(mod.select %>% filter(variable == 'MAE'), aes(x = variable, y = reorder(features, -value), fill = value)) + 
  geom_tile() + 
  scale_fill_viridis_c(option = 'magma') +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = 'Performance metric', y = 'RNN features for quantifying fitness', fill = 'Mean absolute error') +
  geom_text(aes(label = round(value, 2), color = ifelse(value < 0.22, 'white', 'black')), show.legend = FALSE) +
  scale_color_manual(values = c('black', 'white')) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 9, hjust = 0.5, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        panel.border = element_rect(color = 'black', fill = rgb(0,0,0,0)),
        panel.grid = element_line(color = rgb(0,0,0,0)),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 8.5),
        strip.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.key.size = unit(0.8, 'cm'))

out <- ggarrange(s1, ggarrange(s2, s3, nrow = 1), labels = c('a', 'b'), ncol = 1)
ggsave(out, filename = paste0(FIGURE_OUTPUT, 'supplementary_model_selection.png'), width = 10.14, height = 8.72)
       