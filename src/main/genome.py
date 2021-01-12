#!/usr/bin/env python3
''' Description: Runs sliding window analysis on simulated genomes '''

# Import modules 
from simulation import simulateViralEvolution
import model as mod
import numpy as np
import argparse
import random
import time
from virus import Virus, Population
from copy import deepcopy,copy
import math
import random
import pandas as pd
import tensorflow as tf
from itertools import combinations

# Arguments
parser = argparse.ArgumentParser(description = 'Simulation parameters')
parser.add_argument('-m', '--model', help='Path to the trained CNN.')
parser.add_argument('-w', '--fitness', help='Fitness of beneficial mutations.')
parser.add_argument('-u', '--mutrate', help='Fitness of beneficial mutations.')
parser.add_argument('-as', '--alignment_size', default=1000, help = 'Alignment size that CNN was trained on.')
parser.add_argument('-out', '--out', help='Output directory. Requires / at the end.')
parser.add_argument('-p', '--pb', default = 0.1, help='Probability of mutation being beneficial.')
parser.add_argument('-s', '--sort', default = 0.1, help='Sort the haplotypes.')

# Set parameters
args = parser.parse_args()
STEP_SIZE = 50
BUFFER = 50
NUMBER_SUBSAMPLES = 10
GENOME_LENGTH = 29903
ALIGNMENT_SIZE = int(args.alignment_size)
MODEL = str(args.model)
OUTPUT_DIRECTORY = str(args.out)
FITNESS = float(args.fitness)
PROB_BEN = float(args.pb)
MUTRATE = float(args.mutrate)
SORTING = args.sort

# Function
# --------
def popgenStats(haplotypes, n = 200):
	counts = np.array([i for i in np.sum(haplotypes, axis = 0) if i < 200 and i > 0])
	# Watterson's estimator
	S = len(counts)
	if S > 0:
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
	else:
		TajimasD, WuH = 0, 0
	return [TajimasD, WuH]

# Run simulation at fitness of 1.1
# --------------------------------
print('[1] Simulating data.....')
randid = str(random.randint(1, 1e7)) # Assign random id to for saving simulations
sim = simulateViralEvolution(r = 2.05, w = FITNESS, x = 1, probBen = PROB_BEN, mutRate = MUTRATE, initSize = 110, genomeSize = GENOME_LENGTH, gens = 250, maxPopSize = 1e5)

# Sample haplotypes and compute frequencies for positive loci
# -----------------------------------------------------------
haplotypes, modes = mod.sampleData(sims = sim, size = 1000, reps = 1, sort_row = False, sort_col = False)
positiveLoci = sim.positiveLoci
neutralLoci = sim.neutralLoci
positiveFreq = np.sum(np.array(haplotypes[0][:,positiveLoci]), axis = 0)/1000
neutralFreq = np.sum(np.array(haplotypes[0][:,neutralLoci]), axis = 0)/1000

counts_per_virus = np.sum(np.array(haplotypes)[0], axis = 1)
np.save(OUTPUT_DIRECTORY + str(SORTING) + '_' + str(ALIGNMENT_SIZE) + '_' + 'simgenome_' + str(FITNESS) + '_' + str(PROB_BEN) + '_' + str(MUTRATE) + '_counts_' + str(randid), counts_per_virus)


# Build list of indices for sliding window
steps = [(i*STEP_SIZE)+BUFFER for i in range(0, math.floor(GENOME_LENGTH/STEP_SIZE))]
steps = [[steps[i], steps[i+1]+(ALIGNMENT_SIZE - STEP_SIZE)] for i in range(len(steps)-1)]
steps = steps + [[steps[-1][0] + GENOME_LENGTH - steps[-1][1] - BUFFER, steps[-1][1] + GENOME_LENGTH - steps[-1][1] - BUFFER]]
steps = [i for i in steps if i[1] <= GENOME_LENGTH]

# Load model
hapCNN = tf.keras.models.load_model(MODEL)

print('[2] Performing slide window analysis.....')
# Run CNN over subsamples across steps
index1, index2, upper, mean, lower, tajimasD, fayandwuH = [], [], [], [], [], [], []
for i1, i2 in steps:
	print('[----] On to step: ' + str(i2))
	scores = []
	store = []
	for i in range(NUMBER_SUBSAMPLES):
		print('Subsample ' + str(i))
		row_index = random.sample([i for i in range(haplotypes[0].shape[0])], k = 200)
		haplotypeSample = haplotypes[0][row_index,i1:i2]
		print(' - Number mutations in alignment: ' + str(np.sum(haplotypeSample)))
		if SORTING == 'row' or SORTING == 'rowcol':
			haplotypeSample = np.array(haplotypeSample[np.argsort(haplotypeSample.sum(axis=1))[::-1],:].tolist())
		if SORTING == 'col' or SORTING == 'rowcol':
			haplotypeSample = np.array(haplotypeSample[:,np.argsort(haplotypeSample.sum(axis=0))[::-1]].tolist())
		store.append(haplotypeSample)
	store = np.array(store)
	tajimasD = tajimasD + [np.mean([popgenStats(k)[0] for k in store])]
	fayandwuH = fayandwuH + [np.mean([popgenStats(k)[1] for k in store])]
	store = store.reshape(store.shape + (1,))
	scores = hapCNN.predict(store)
	upper = upper + [np.quantile(scores, 0.975)]
	mean = mean + [np.mean(scores)]
	lower = lower + [np.quantile(scores, 0.025)]
	index1 = index1 + [i1]
	index2 = index2 + [i2]

output = pd.DataFrame({'index1': index1, 'index2':index2, 'upper':upper, 'mean':mean, 'lower':lower, 'Tajima\'s D':tajimasD, 'Fay and Wu\'s H':fayandwuH})
output.to_csv(OUTPUT_DIRECTORY + str(SORTING) + '_' + str(ALIGNMENT_SIZE) + '_' + 'simgenome_' + str(FITNESS) + '_' + str(PROB_BEN) + '_' + str(MUTRATE) + '_predictions_' + str(randid) + '.csv', index = False)

output2 = pd.DataFrame({'loci': positiveLoci + neutralLoci, 'VAF': positiveFreq.tolist() + neutralFreq.tolist(), 'mode': ['Beneficial' for i in range(len(positiveLoci))] + ['Neutral' for i in range(len(neutralLoci))]})
output2.to_csv(OUTPUT_DIRECTORY + str(SORTING) + '_' + str(ALIGNMENT_SIZE) + '_' + 'simgenome_' + str(FITNESS) + '_' + str(PROB_BEN) + '_' + str(MUTRATE) + '_loci_' + str(randid) + '.txt')