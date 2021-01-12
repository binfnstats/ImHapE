#!/bin/env python3
''' Description: This script takes CNN estimates from sliding window analysis of full 29903 bp simulated genomes and trains an LSTM to predict if the window actually contains a beneficial mutation '''

OUTPUT = '/OUTPUT_PATH/mutation_counts_simgenome.csv'

import os
import numpy as np
import pandas as pd 

filename, mut_counts = [], []
for j in [i for i in os.listdir() if 'count' in i]:
	a = np.load(j)
	filename, mut_counts = filename + [j for k in range(1000)], mut_counts + a.tolist()

c = pd.DataFrame({'filename':filename, 'mut_counts':mut_counts})
c['fitness'] = [i.split('_')[3] for i in c['filename']]
c['mutrate'] = [i.split('_')[5] for i in c['filename']]
c.to_csv(OUTPUT)