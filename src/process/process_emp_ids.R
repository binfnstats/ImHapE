#!/bin/env Rscript

# Generate identifier tables for running empirical analysis 
# =========================================================
library(data.table)
library(dplyr)

OUTPUT_DIR = '/OUTPUT_DIR'

# Build a list of all ids to be analyzed across timepoints
# ---------------------------------------------------------------------
cog <- fread('cog_meta_withweeks.csv') %>% filter(prop_masked <= 0.018)
for (week in unique(cog$epi_week)) {
	print(week)
	sv <- cog %>% filter(epi_week == week) %>% select(filename, epi_week)
	outname <- paste0(OUTPUT_DIR, 'COG_', week, '_.csv')
	write.table(sv, file = outname, quote = FALSE, row.names = FALSE, sep =',')
}

gisaid <- fread('aligned_meta.csv') %>% filter(prop_masked <= 0.001)
for (reg in c('Asia', 'Europe', 'North America')) {
	print(reg)
	for (mt in c(3, 4, 5, 6, 7)) {
		print(mt)
		sv <- gisaid %>% filter(region == reg & month == mt) %>% select(filename, region, month)
		print(head(sv))
		outname <- paste0(OUTPUT_DIR, 'GISAID_', gsub(' ', '', reg), '_', mt, '_.csv')
		write.table(sv, file = outname, quote = FALSE, row.names = FALSE, sep =',')
	}
}