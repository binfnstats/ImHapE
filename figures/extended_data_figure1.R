# SARS-CoV-2 Figures
# ------------------
# -----
# Extended Data Figure 1: Relationship between CNN and RNN estimates and VAF
# -----
setwd('SET_WORKING_DIRECTORY_TO_FOLDER_CONTAINING_MODEL_VALIDATION_FOLDER')
OUTPUT <- 'PATH_TO_FIGURE_DIRECTORY/'

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(caret)

# MODEL HISTORY: Examine model history across different window sizes and sorting methods
# ============================================================================================
dt2 <- data.table()
for (i in list.files('./model_validation/')) {
  add <- fread(paste('./model_validation/', i, sep = '')) %>% dplyr::select(-V1)
  add$`Window Size` <- strsplit(i, '_')[[1]][2]
  dt2 <- rbind(dt2, add)
}
dt2$`Window Size` <- ifelse(dt2$`Window Size` == 1000, 'Window: 1000bp (Row/Column Sorted)', ifelse(dt2$`Window Size` == 2500, 'Window: 2500bp (No Sorting)', 'Window: 5000bp (No Sorting)'))
dt2$modes <- ifelse(dt2$modes == 'neutral', 'N', dt2$modes)
dt2$modes <- factor(dt2$modes, levels = c('N', names(table(dt2$modes))[1:10]))
dt2$class <- ifelse(dt2$class == 1, 'Positive', 'Neutral')

B1 <- ggplot(dt2, aes(x = modes)) + 
  geom_bar(aes(fill = class), position = 'fill') +
  facet_wrap(~`Window Size`, nrow = 1) +
  theme(panel.border = element_rect(color = "black", fill = rgb(0,0,0,0)),
        panel.background = element_rect(fill = rgb(0,0,0,0)),
        strip.background = element_blank(),
        panel.grid = element_line(color = rgb(0,0,0,0.05)),
        strip.text = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8.5),
        plot.title = element_text(size = 8.5),
        legend.key.size = unit(0.8, 'lines'),
        legend.position = 'top',
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = 'Fitness (1 + s)', y = 'Proporition(%)', fill = 'Class') +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_npg()

B2 <- ggplot(dt2, aes(x = modes, y = scores)) + 
  #geom_jitter(width = 0.05, alpha = 0.3, size = 0.8, show.legend = FALSE) +
  geom_violin(alpha = 0.5, scale = 'width', draw_quantiles = 0.5) +
  facet_wrap(~`Window Size`, nrow = 1) +
  theme(panel.border = element_rect(color = "black", fill = rgb(0,0,0,0)),
        panel.background = element_rect(fill = rgb(0,0,0,0)),
        strip.background = element_blank(),
        panel.grid = element_line(color = rgb(0,0,0,0.05)),
        #strip.text = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8.5),
        plot.title = element_text(size = 8.5),
        legend.key.size = unit(0.8, 'lines'),
        legend.position = 'bottom',
        strip.text = element_blank(),
        legend.key = element_rect(fill = 'white'),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8)) +
  labs(x = 'Fitness (1 + s)', y = 'P(Selection)')

B <- ggarrange(B1, B2, common.legend = TRUE, legend = 'bottom', ncol = 1)

# Save plots
ggsave(B, file = paste(OUTPUT, 'extended_data_figure1.png', sep = ''), width = 8.73, height = 5.35, units = 'in')
