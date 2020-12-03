#!/usr/bin/env python3
''' Description: Scans aligned haplotypes for signatures of selection. Subsamples from the population to build confidence intervals.'''

import pandas as pd
import numpy as np
import tensorflow as tf
import math
import random
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'Parameters for running CNN-based analysis')
parser.add_argument('-gl', '--genome_length', help='Size of organism genome (e.g. COVID ~ 29903')
parser.add_argument('-ss', '--step_size', help='Sliding window step size.')
parser.add_argument('-as', '--alignment_size', help='Length of haplotype alignments CNN was trained on.')
parser.add_argument('-bs', '--buffer_size', help='Starting and ending buffer size')
parser.add_argument('-ns', '--number_subsamples', help='Number of subsamples to perform for each sliding window.')
parser.add_argument('-g', '--group', type = str, help='The name of the column in geography data frame.')
parser.add_argument('-s', '--samples', type = str, help='The specific value in the group column to filter geography data frame.')
parser.add_argument('-mon', '--months', default=0, type = int, help='The specific month to evaluate. Default of 0 equals all months.')
parser.add_argument('-geo', '--geography', help='Geography data frame path.')
parser.add_argument('-cnn', '--cnn_model', help='The tensorflow CNN path.')
parser.add_argument('-rnn', '--rnn_model', help='The tensorflow CNN path.')
parser.add_argument('-d', '--haplotype_directory', help='The directory containing the binary encoded haplotypes.')
parser.add_argument('-out', '--output', help='The output directory.')
parser.add_argument('-back', '--backwards', type = str, default = 'False', help='Runs the sliding window in the opposite direction. Combining forward and backward slides refines predictions.')

args = parser.parse_args()
GENOME_LENGTH = int(args.genome_length)
STEP_SIZE = int(args.step_size)
ALIGNMENT_SIZE = int(args.alignment_size)
BUFFER = int(args.buffer_size)
NUMBER_SUBSAMPLES = int(args.number_subsamples)
GROUP = args.group
SAMPLES = args.samples
MONTH = int(args.months)
GEOGRAPHY = str(args.geography)
CNN_MODEL = str(args.cnn_model)
RNN_MODEL = str(args.rnn_model)
DIRECTORY = str(args.haplotype_directory)
OUTPUT_DIRECTORY = str(args.output)
BACKWARDS = False if args.backwards == 'False' else True

# Example
#python3 /.mounts/labs/awadallalab/private/touellette/covid/src/slide.py -gl 29903 -ss 5 -as 1000 -bs 50 -ns 10 -g region -s South Asia -mon 0 -geo /.mounts/labs/awadallalab/private/touellette/covid/data/covid_ids.csv -m /.mounts/labs/awadallalab/private/touellette/covid/models_final/positive_ALL_2020-10-08-20-02_cnn.tf -d /.mounts/labs/awadallalab/private/touellette/covid/haplotypes/ -out /.mounts/labs/awadallalab/private/touellette/covid/sliding_output/

# Arguments
# GENOME_LENGTH = 29903
# STEP_SIZE = 100
# BUFFER = 50
# ALIGNMENT_SIZE = 1000
# NUMBER_SUBSAMPLES = 5
# GROUP = 'country'
# SAMPLES = 'Canada'
# MONTH = 4
# GEOGRAPHY = '/.mounts/labs/awadallalab/private/touellette/covid/data/covid_ids.csv'
# MODEL = '/.mounts/labs/awadallalab/private/touellette/covid/models_final/positive_ALL_2020-10-08-20-02_cnn.tf'
# DIRECTORY = '/.mounts/labs/awadallalab/private/touellette/covid/haplotypes/'
# RNN_MODEL = '/.mounts/labs/awadallalab/private/touellette/sars-cov-2/models/RNN_2500.tf'
# CNN_MODEL = '/.mounts/labs/awadallalab/private/touellette/covid/models_final/positive_ALL_2020-10-08-20-02_cnn.tf'

# Generating predictions along the covid genome (29903bp)
# -------------------------------------------------------
# i. Skip/buffer first and last bases to avoid biases due to sequencing error/etc.
# ii. Randomly sample 200 haplotypes from specified group
# iii. Make prediction
# iv. Repeat ii - iii, N times to get a subsampled confidence interval
# v. Slide k base pairs, and repeat ii - iv 

# Run script
# ----------

# Get file ids for specific group
geo = pd.read_csv(GEOGRAPHY)
geo = geo[geo[GROUP] == SAMPLES]
if MONTH != 0:
	geo = geo[geo['month'] == MONTH]

# Build list of indices for sliding window
steps = [(i*STEP_SIZE)+BUFFER for i in range(0, math.floor(GENOME_LENGTH/STEP_SIZE))]
steps = [[steps[i], steps[i+1]+(ALIGNMENT_SIZE - STEP_SIZE)] for i in range(len(steps)-1)]
steps = steps + [[steps[-1][0] + GENOME_LENGTH - steps[-1][1] - BUFFER, steps[-1][1] + GENOME_LENGTH - steps[-1][1] - BUFFER]]
steps = [i for i in steps if i[1] <= GENOME_LENGTH]

if BACKWARDS == True:
	steps = steps[::-1]

# Load model
hapCNN = tf.keras.models.load_model(CNN_MODEL)
hapRNN = tf.keras.models.load_model(RNN_MODEL)

# Run CNN over subsamples across steps
index1, index2, upper, mean, lower = [], [], [], [], []
rnn_array = []
iteration = 0
for i1, i2 in steps:
	iteration += 1
	print('[----] On to step: ' + str(i2))
	scores = []
	store = []
	for i in range(NUMBER_SUBSAMPLES):
	        print('Subsample ' + str(i))
	        filenames = random.choices(list(geo['filename']), k = 200)
	        haplotypes = np.array([list(''.join(open(DIRECTORY + i, 'r').readlines()).rstrip())[i1:i2] for i in filenames]).astype(int)
	        #haplotypes = haplotypes[:,np.argsort(haplotypes.sum(axis=0))[::-1]]
	        #haplotypes = haplotypes.reshape((1,) + haplotypes.shape + (1,))
	        print(' - Number mutations in alignment: ' + str(np.sum(haplotypes)))
	        store.append(haplotypes)
	        #score = hapCNN.predict(haplotypes)
	        #scores = scores + [score] 
	store = np.array(store)
	store = store.reshape(store.shape + (1,))
	scores = hapCNN.predict(store)
	upper = upper + [np.quantile(scores, 0.975)]
	mean = mean + [np.mean(scores)]
	lower = lower + [np.quantile(scores, 0.025)]
	index1 = index1 + [i1]
	index2 = index2 + [i2]
	if iteration == 1:
		rnn_array = scores
	else:
		rnn_array = np.hstack((rnn_array, scores))

# Compute RNN predictions
rnn_array2 = rnn_array.reshape(rnn_array.shape + (1,))
rnn_out = hapRNN.predict(rnn_array2)
rnn_mean = np.mean(rnn_out, axis = 0).T[0]
rnn_upper = np.quantile(rnn_out, 0.975, axis = 0).T[0]
rnn_lower = np.quantile(rnn_out, 0.025, axis = 0).T[0]
output = pd.DataFrame({'index1': index1, 'index2':index2, 'upper':upper, 'mean':mean, 'lower':lower, 'rnn_upper':rnn_upper, 'rnn_mean':rnn_mean, 'rnn_lower':rnn_lower})

output.to_csv(OUTPUT_DIRECTORY + str(ALIGNMENT_SIZE) + '_predictions_' + str(GROUP) + '_' + str(MONTH) + '_' + ''.join(SAMPLES.split()) + '.csv', index = False)