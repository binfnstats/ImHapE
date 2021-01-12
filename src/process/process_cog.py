#!/usr/bin/env python3
''' Description: Process the aligned COG sequences and link back to GISAID meta data'''

import numpy as np
import pandas as pd

meta_file = '/COG_PATH/cog_metadata_withweeks.csv'
alignment_fasta = '/COG_PATH/cog_alignment.fasta'
haplotype_output_directory = '/COG_PATH/cog_sequences/'
meta_output_directory = '/OUTPUT_PATH/data/'

def processSequences(alignment_fasta, output_directory):
	aligned = [seq.rstrip() for seq in open(alignment_fasta,'r')]
	identifiers = [aligned[i].split('>')[1] for i in range(0, len(aligned), 2)]
	sequences = [aligned[i] for i in range(1, len(aligned), 2)]
	filename = [i.replace('/', '-') + '.txt' for i in identifiers]
	masked_sites = []
	prop_masked = []
	for i in range(len(sequences)):
		seq, fname = sequences[i], filename[i]
		masked = np.sum([i == 'N' for i in seq])
		prop = masked / 29903
		masked_sites, prop_masked = masked_sites + [masked], prop_masked + [prop]
		with open(output_directory + fname, "w") as text:
			text.write(seq)
	data = pd.DataFrame({'strain': identifiers, 'filename': filename, 'masked_sites': masked_sites, 'prop_masked': prop_masked})
	return (data)

def combineMeta(data, meta_file, delim = '\t'):
	meta = pd.read_csv(meta_file, delimiter = delim)
	meta = meta.loc[:,['gisaid_epi_isl', 'strain', 'date', 'country', 'division', 'region', 'GISAID_clade', 'Nextstrain_clade', 'pangolin_lineage']]
	combined = data.merge(meta, on = 'strain')
	combined = combined.loc[combined['date'] != '2020',:] # Only get samples with exact sequencing dates
	combined['month'] = [i.split('-')[1] for i in combined['date']]
	return (combined)

df = processSequences(alignment_fasta = alignment_fasta, output_directory = haplotype_output_directory)
combined = combineMeta(data = df, meta_file = meta_file)
df.to_csv
combined.to_csv(meta_output_directory + 'cog_meta.csv')