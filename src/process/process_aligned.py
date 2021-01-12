#!/usr/bin/env python3
''' Description: Takes a directory of EMBOSS stretcher aligned sequences and generates a pseudo-VCF file and an identifier file'''

import pandas as pd
import numpy as np
import os

meta_file = '/META_PATH/metadata_2020-12-18_09-25.tsv'
alignment_dir = '/ALIGNMENT_PATH/alignments/'
haplotype_output_directory = '/HAPLOTYPE_OUTPUT_PATH/gisaid_sequences/'
meta_output_directory = '/META_OUTPUT/data/'

def processStretched(alignment_dir, output_directory):
	os.chdir(alignment_dir)
	masked_sites, prop_masked, filenames, gisaid_epi_isl = [], [], [], []
	for file in os.listdir():
		print(file)
		try:
			alignment = [seq.strip() for seq in open(file, 'r')]
			end = alignment.index('>EMBOSS_002')
			alignment = ''.join(alignment[end+1:])
			print(len(alignment))
			if len(alignment) == 29903:
				masked = np.sum([i == 'N' for i in alignment])
				prop = masked/29903
				masked_sites, prop_masked, filenames, gisaid_epi_isl = masked_sites + [masked], prop_masked + [prop], filenames + [file], gisaid_epi_isl + [file.split('+')[1]]
				with open(output_directory + file, "w") as text:
					text.write(alignment)
		except:
			pass
	data = pd.DataFrame({'filename': filenames, 'masked_sites': masked_sites, 'prop_masked':prop_masked, 'gisaid_epi_isl':gisaid_epi_isl})
	return(data)

def combineMeta(data, meta_file, delim = '\t'):
	meta = pd.read_csv(meta_file, delimiter = delim)
	meta = meta.loc[:,['gisaid_epi_isl', 'strain', 'date', 'country', 'division', 'region', 'GISAID_clade', 'Nextstrain_clade', 'pangolin_lineage']]
	combined = data.merge(meta, on = 'gisaid_epi_isl')
	combined = combined.loc[combined['date'] != '2020',:] # Only get samples with exact sequencing dates
	combined['month'] = [i.split('-')[1] for i in combined['date']]
	return (combined)


df = processStretched(alignment_dir = alignment_dir, output_directory = haplotype_output_directory)
combined = combineMeta(data = df, meta_file = meta_file)
combined.to_csv(meta_output_directory + 'aligned_meta.csv')