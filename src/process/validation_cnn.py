#!/usr/bin/env python3
''' Description: This script computes summary statistics and CNN performance across the simulated haplotypes '''

import tensorflow as tf
import pandas as pd
import numpy as np
import datetime
import argparse
import os
import datetime
from itertools import combinations

parser = argparse.ArgumentParser(description = 'Sumstats input and output.')
parser.add_argument('-dir', '--directory', help='Haplotype directory')
parser.add_argument('-out', '--output', help='Output path file name')
parser.add_argument('-mod', '--model', help='CNN model')
parser.add_argument('-s', '--sort', default='none', help='Sorting')

args = parser.parse_args()
DIRECTORY = args.directory
OUTPUT_PATH = args.output
CNN_MODEL = args.model
SORTING = args.sort

ImHapE = tf.keras.models.load_model(CNN_MODEL)

# Compute summary statistic metrics
# ---------------------------------
def mutateStats(haplotypes):
	average_mutations_per_virus = np.mean(np.sum(haplotypes, axis = 0))
	average_mutations_per_site = np.mean(np.sum(haplotypes, axis = 1))
	mutations_per_alignment = np.sum(haplotypes)	
	return [average_mutations_per_virus, average_mutations_per_site, mutations_per_alignment]

def popgenStats(haplotypes, n = 200):
	counts = np.array([i for i in np.sum(haplotypes, axis = 0) if i < 200 and i > 0])
	# Watterson's estimator
	S = len(counts)
	an = np.sum([1/i for i in range(1, n)])
	theta_w = S / an
	# Tajima's pi
	vafs = counts / n
	d = np.sum(vafs * (1 - vafs)) * (n^2)
	theta_pi = d / len([j for j in combinations([i for i in range(n)], 2)])
	# Variance for Tajima's D
	a1 = np.sum([1/i for i in range(1, n)])
	a2 = np.sum([1/(i*i) for i in range(1, n)])
	b1 = (n + 1) / (3 * (n - 1))
	b2 = (2 * ((n*n) + n + 3)) / (9 * n * (n - 1))
	c1 = b1 - (1 / a1)
	c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1*a1))
	e1 = c1 / a1
	e2 = c2 / ((a1*a1) + a2)
	seD = np.sqrt((e1 *S) + ((e2 * S) * (S - 1)))
	TajimasD = (theta_pi - theta_w) / seD
	# Fay and Wu's H
	theta_l = np.sum(vafs) * (n / (n - 1))
	bn = np.sum([1/(i*i) for i in range(1, n)])
	bn1 = np.sum([1/(i*i) for i in range(1, n+1)])
	theta2 = (S * (S - 1)) / ((an*an) + bn)
	seH = np.sqrt((((n - 2) / (6*(n - 1))) * theta_w) + (theta2 * (((18 * (n*n) * (3 * n + 2) * bn1) - (88 * (n*n*n) + 9 * (n*n) - 13*n + 6))/(9 * n * ((n - 1)*(n-1))))))
	WuH = 2 * ((theta_pi - theta_l) / seH)
	return [TajimasD, WuH]

# Organize and save output
# ------------------------
os.chdir(DIRECTORY)

average_mutations_per_virus = [] 
average_mutations_per_site = []
mutations_per_alignment = []
tajimasD = []
faywuH = []
mutrate = []
probben = []
fitness = []
ids = []
window = []
scores = []
count = 1
for file in os.listdir():
	print(count)
	count += 1
	haplotypes = np.load(file)[0]
	if SORTING == 'row' or SORTING == 'rowcol':
		haplotypes = np.array(haplotypes[np.argsort(haplotypes.sum(axis=1))[::-1],:].tolist())
	if SORTING == 'col' or SORTING == 'rowcol':
		haplotypes = np.array(haplotypes[:,np.argsort(haplotypes.sum(axis=0))[::-1]].tolist())
	descriptive, sumstats = mutateStats(haplotypes), popgenStats(haplotypes)
	ids = ids + [file]
	# Get descriptive data
	average_mutations_per_virus = average_mutations_per_virus + [descriptive[0]]
	average_mutations_per_site = average_mutations_per_site + [descriptive[1]]
	mutations_per_alignment = mutations_per_alignment + [descriptive[2]]
	window = window + [int(file.split('_')[0])]
	# Get simulation parameters
	mutrate = mutrate + [float(file.split('_')[1])/float(file.split('_')[0])]
	probben = probben + [float(file.split('_')[2])]
	fitness = fitness + [float(file.split('_')[3])]
	# Get summary statistics
	tajimasD = tajimasD + [sumstats[0]]
	faywuH = faywuH + [sumstats[1]]
	# Get CNN prediction
	haplotypes = haplotypes.reshape((1,) + haplotypes.shape + (1,))
	scores = scores + ImHapE.predict(haplotypes)[0].tolist()

df = pd.DataFrame({'id':ids, 'Window size (bp)':window, 'Average mutations per haplotype':average_mutations_per_virus, 'Average mutations per site':average_mutations_per_site, 'Mutations per alignment':mutations_per_alignment, 'Tajima\'s D':tajimasD, 'Fay and Wu\'s H':faywuH, 'CNN':scores, 'Mutation rate (per site/gen)':mutrate, 'P(Beneficial)':probben, 'Fitness (1 + s)':fitness})

df.to_csv(OUTPUT_PATH + '.csv', index = False)

