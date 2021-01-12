import numpy as np 
import os
import pandas as pd

EMPIRICAL_ANALYSIS = '/LISTS_OF_EMPIRICAL_SAMPLES_PER_TIMEPOINT/empirical_analysis'
GISAID_WUHAN = '/GISAID_WUHAN_REFERENCE/haplotypes_gisaid_wuhan/'
COG_ENGLAND = '/COG_ENGLAND_REFERENCE/haplotypes_cog_england/'

# COG VAFs over time
# ------------------
cog_vafs = pd.DataFrame({'position':[i+1 for i in range(29903)]})
os.chdir(EMPIRICAL_ANALYSIS)
for i in [j for j in os.listdir() if 'COG' in j]:
	print(i)
	ids = pd.read_csv(i)
	week = i.split('_')[1]
	n = len(ids['filename'].tolist())
	counts = np.sum(np.array([[int(r) for r in list([k.rstrip() for k in open(COG_ENGLAND + file, 'r')][0])] for file in ids['filename'].tolist()]), axis = 0)
	vafs = counts / n
	cog_vafs[week] = vafs

cog_vafs.to_csv('cog_england_vafs_by_week.csv', index = False)

gisaid_vafs = pd.DataFrame({'position':[i+1 for i in range(29903)]})
for i in [j for j in os.listdir() if 'GISAID' in j]:
	print(i)
	ids = pd.read_csv(i)
	region = i.split('_')[1]
	month = i.split('_')[2]
	both = region + '_' + month
	n = len(ids['filename'].tolist())
	counts = np.sum(np.array([[int(r) for r in list([k.rstrip() for k in open(GISAID_WUHAN + file, 'r')][0])] for file in ids['filename'].tolist()]), axis = 0)
	vafs = counts / n
	gisaid_vafs[both] = vafs

gisaid_vafs.to_csv('gisaid_wuhan_vafs_by_week.csv', index = False)