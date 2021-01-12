#!/usr/bin/env python3
''' Description: Evaluating quality of the sequences from GISAID and COG '''

cog_sequences = '/COG_PATH/cog_sequences/'
cog_haplotypes_wuhan = '/COG_WUHAN_REFERENCE/haplotypes_cog_wuhan/'
cog_haplotypes_england = '/COG_ENGLAND_REFERENCE/haplotypes_cog_england/'
gisaid_sequences = '/GISAID_WUHAN_REFERENCE/gisaid_sequences/'
gisaid_haplotypes = '/GISAID_HAPLOTYPES/haplotypes_gisaid_wuhan/'

import numpy as np 
import pandas as pd 
import os

def examineSites(target_dir):
	masked = [0 for i in range(29903)]
	nonstandard = [0 for i in range(29903)]
	os.chdir(target_dir)
	for i in os.listdir():
		seq = [k.rstrip() for k in open(i, 'r')][0]
		if len(seq) == 29903:
			for j in range(len(seq)):
				if seq[j] == 'N':
					masked[j] += 1
				if seq[j] not in ['A', 'G', 'C', 'T']:
					nonstandard[j] += 1
	data = pd.DataFrame({'n_masked':masked, 'n_nonstandard':nonstandard, 'position':[k for k in range(1, 29904)]})
	return (data)

def mutationCounts(target_dir):
	os.chdir(target_dir)
	mut_counts = []
	for i in os.listdir():
		print(i)
		seq = [k.rstrip() for k in open(i, 'r')][0]
		if len(seq) == 29903:
			seq = np.array([int(k) for k in seq])
			nmuts = np.sum(seq)
			mut_counts = mut_counts + [nmuts]
	data = pd.DataFrame({'filename':os.listdir(), 'mut_counts':mut_counts})
	return (data)

# Write files
# -----------
meta_output_directory = '/OUTPUT_PATH/data/'

cog_sites = examineSites(cog_sequences)
cog_sites.to_csv(meta_output_directory + 'sites_masked_nonstandard_COG.csv')

gisaid_sites = examineSites(gisaid_sequences)
gisaid_sites.to_csv(meta_output_directory + 'sites_masked_nonstandard_GISAID.csv')

cog_mc_wuhan = mutationCounts(cog_haplotypes_wuhan)
cog_mc_england = mutationCounts(cog_haplotypes_england)
gisaid_mc_wuhan = mutationCounts(gisaid_haplotypes)

# Write dfs
# ---------
cog_mc_wuhan.to_csv(meta_output_directory + 'mutation_counts_COG_wuhan.csv')
cog_mc_england.to_csv(meta_output_directory + 'mutation_counts_COG_england.csv')
gisaid_mc_wuhan.to_csv(meta_output_directory + 'mutation_counts_GISAID_wuhan.csv')
