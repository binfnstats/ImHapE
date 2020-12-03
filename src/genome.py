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

# Arguments
parser = argparse.ArgumentParser(description = 'Simulation parameters')
parser.add_argument('-m', '--model', help='Path to the trained CNN.')
parser.add_argument('-w', '--fitness', help='Fitness of beneficial mutations.')
parser.add_argument('-as', '--alignment_size', default=1000, help = 'Alignment size that CNN was trained on.')
parser.add_argument('-out', '--out', help='Output directory. Requires / at the end.')
parser.add_argument('-s', '--sort', help='One of col, row, row_col, or none.')
parser.add_argument('-p', '--pb', default = 0.1, help='Probability of mutation being beneficial.')

# Set parameters
args = parser.parse_args()
STEP_SIZE = 100
BUFFER = 50
NUMBER_SUBSAMPLES = 10
GENOME_LENGTH = 29903
ALIGNMENT_SIZE = int(args.alignment_size)
MODEL = str(args.model)
OUTPUT_DIRECTORY = str(args.out)
SORT_DATA = str(args.sort)
FITNESS = float(args.fitness)
PROB_BEN = float(args.pb)

# Run simulation at fitness of 1.1
# --------------------------------
print('[1] Simulating data.....')
randid = str(random.randint(1, 1e7)) # Assign random id to for saving simulations
sim = simulateViralEvolution(r = 2.05, w = FITNESS, x = 1, probBen = PROB_BEN, mutRate = 0.03, initSize = 110, genomeSize = 29903, gens = 250, maxPopSize = 1e5)

# Sample haplotypes and compute frequencies for positive loci
# -----------------------------------------------------------
haplotypes, modes = mod.sampleData(sims = sim, size = 1000, reps = 1, sort_row = False, sort_col = False)
positiveLoci = sim.positiveLoci
neutralLoci = sim.neutralLoci
positiveFreq = np.sum(np.array(haplotypes[0][:,positiveLoci]), axis = 0)/1000
neutralFreq = np.sum(np.array(haplotypes[0][:,neutralLoci]), axis = 0)/1000

# Build list of indices for sliding window
steps = [(i*STEP_SIZE)+BUFFER for i in range(0, math.floor(GENOME_LENGTH/STEP_SIZE))]
steps = [[steps[i], steps[i+1]+(ALIGNMENT_SIZE - STEP_SIZE)] for i in range(len(steps)-1)]
steps = steps + [[steps[-1][0] + GENOME_LENGTH - steps[-1][1] - BUFFER, steps[-1][1] + GENOME_LENGTH - steps[-1][1] - BUFFER]]
steps = [i for i in steps if i[1] <= GENOME_LENGTH]

# Load model
hapCNN = tf.keras.models.load_model(MODEL)

print('[2] Performing slide window analysis.....')
# Run CNN over subsamples across steps
index1, index2, upper, median, lower = [], [], [], [], []
for i1, i2 in steps:
	print('[----] On to step: ' + str(i2))
	scores = []
	store = []

	for i in range(NUMBER_SUBSAMPLES):
		print('Subsample ' + str(i))
		row_index = random.choices([i for i in range(haplotypes[0].shape[0])], k = 200)
		haplotypeSample = haplotypes[0][row_index,i1:i2]
		if SORT_DATA == 'col' or SORT_DATA == 'row_col':
			haplotypeSample = haplotypeSample[:,np.argsort(haplotypeSample.sum(axis=0))[::-1]]
		if SORT_DATA == 'row' or SORT_DATA == 'row_col':
			haplotypeSample = haplotypeSample[np.argsort(haplotypeSample.sum(axis=1))[::-1],:]
		print(' - Number mutations in alignment: ' + str(np.sum(haplotypeSample)))
		store.append(haplotypeSample)

	store = np.array(store)
	store = store.reshape(store.shape + (1,))
	scores = hapCNN.predict(store)
	upper = upper + [np.quantile(scores, 0.975)]
	median = median + [np.mean(scores)]
	lower = lower + [np.quantile(scores, 0.025)]
	index1 = index1 + [i1]
	index2 = index2 + [i2]

output = pd.DataFrame({'index1': index1, 'index2':index2, 'upper':upper, 'mean':median, 'lower':lower})
output.to_csv(OUTPUT_DIRECTORY + str(ALIGNMENT_SIZE) + '_' + 'simgenome_' + str(FITNESS) + '_' + 'predictions_' + str(randid) + '.csv', index = False)

output2 = pd.DataFrame({'loci': positiveLoci + neutralLoci, 'VAF': positiveFreq.tolist() + neutralFreq.tolist(), 'mode': ['Beneficial' for i in range(len(positiveLoci))] + ['Neutral' for i in range(len(neutralLoci))]})
output2.to_csv(OUTPUT_DIRECTORY + str(ALIGNMENT_SIZE) + '_' + 'simgenome_' + str(FITNESS) + '_' + 'loci_' + str(randid) + '.txt')