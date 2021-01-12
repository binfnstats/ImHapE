#!/usr/bin/env python3
''' Description: Scans aligned haplotypes for signatures of selection. Subsamples from the population to build confidence intervals.'''

import pandas as pd
import numpy as np
import tensorflow as tf
import math
import random
import argparse
from itertools import combinations

# Arguments
parser = argparse.ArgumentParser(description = 'Parameters for running CNN-based analysis')
parser.add_argument('-d', '--haplotype_directory', help='The directory containing the binary encoded haplotypes.')
parser.add_argument('-t', '--target', help='Target file.')
parser.add_argument('-out', '--output', help='The output directory.')

args = parser.parse_args()
DIRECTORY = str(args.haplotype_directory)
OUTPUT_DIRECTORY = str(args.output)
TARGET_FILE = str(args.target)

STEP_SIZE = 50
BUFFER = 50
NUMBER_SUBSAMPLES = 25
GENOME_LENGTH = 29903
ALIGNMENT_SIZE = 2500
BUFFER = 50
CNN_MODEL = 'CNN_MODEL example ==> newest_2500_5e-6_row_2020-12-28-07-07_cnn.tf'
RNN_CLASSIFICATION_MODEL = 'RNN_CLASSIFCATION_MODE example ===> classification_0.0_wuH_CNN_rnn.tf'
RNN_REGRESSION_MODELS = 'PATH_TO_FOLDER_OF_CONTAINING_ALL_REGRESSION_MODELS_USED_IN_ENSEMBLE_INCLUDING_RNN_ENSEMBLE_JOBLIB'

# Pop gen stats function
# ----------------------
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
	if np.isnan([TajimasD]) == True and np.isnan([WuH]) == True:
		return [-3, 0]
	else:
		return [TajimasD, WuH]

def ensemblePredict(tajD, wuH, CNN, rnn_models = ['0.0_tajD', '0.001_tajD_CNN', '0.002_CNN', '0.0_wuH_CNN', '0.003_CNN', '0.0_CNN', '0.0_multi', '0.0_tajD_CNN', '0.0_tajD_wuH', '0.001_tajD'], model_directory = RNN_REGRESSION_MODELS):
	from joblib import load
	ensembled = load(RNN_REGRESSION_MODELS + 'rnn_ensemble.joblib')
	regression = ['regression_' + i + '_rnn.tf' for i in rnn_models]
	df = pd.DataFrame({})
	for mod in regression:
		print(mod)
		vars()[mod] = tf.keras.models.load_model(RNN_REGRESSION_MODELS + mod)
		identifier = mod.split('regression_')[1].split('_rnn.tf')[0]
		fts = []
		if 'CNN' in mod or 'multi' in mod:
			fts.append(CNN)
		if 'wuH' in mod or 'multi' in mod:
			fts.append(wuH)
		if 'tajD' in mod or 'multi' in mod:
			fts.append(tajD)
		if 'multi' in mod:
			fts = np.transpose(fts)
		elif len(mod.split('_')) > 4:
			fts = np.transpose(fts)
		print('past features')
		if len(mod.split('_')) > 4 or 'multi' in mod:
			fts = np.array(fts)
			fts = fts.reshape((1,) + fts.shape)
			pred = vars()[mod].predict(np.array(fts))
		else:
			fts = np.array(fts)
			fts = fts.reshape(fts.shape + (1,))
			pred = vars()[mod].predict(fts)
		df[identifier] = pred[0]
	preds = ensembled.predict(df)
	return preds[0]

# Run script
# ----------

# Get file ids for specific group
geo = pd.read_csv(TARGET_FILE)

# Build list of indices for sliding window
steps = [(i*STEP_SIZE)+BUFFER for i in range(0, math.floor(GENOME_LENGTH/STEP_SIZE))]
steps = [[steps[i], steps[i+1]+(ALIGNMENT_SIZE - STEP_SIZE)] for i in range(len(steps)-1)]
steps = steps + [[steps[-1][0] + GENOME_LENGTH - steps[-1][1] - BUFFER, steps[-1][1] + GENOME_LENGTH - steps[-1][1] - BUFFER]]
steps = [i for i in steps if i[1] <= GENOME_LENGTH]

# Load model
hapCNN = tf.keras.models.load_model(CNN_MODEL)
hapRNN_classification = tf.keras.models.load_model(RNN_CLASSIFICATION_MODEL)

# Run CNN over subsamples across steps
rnn_classification, rnn_regression = [], []
for j in range(NUMBER_SUBSAMPLES):
	index1, index2, upper, mean, lower, wuH, tajD  = [], [], [], [], [], [], []
	#tajD, wuH = [], []
	for i1, i2 in steps:
		print('[----] On to step: ' + str(i2))
		scores = []
		for i in range(1):
		        print('Subsample ' + str(j))
		        if len(list(geo['filename'])) < 200:
		        	filenames = random.choices(list(geo['filename']), k = 200)
		        else:
		        	filenames = random.sample(list(geo['filename']), k = 200)
		        haplotypes = np.array([list(''.join(open(DIRECTORY + i, 'r').readlines()).rstrip())[i1:i2] for i in filenames]).astype(int)
		        haplotypes[:,np.where(np.sum(haplotypes, axis = 0) == 200)[0]] = 0 # Remove fixed sites
		        haplotypes = haplotypes[np.argsort(haplotypes.sum(axis=1))[::-1],:].tolist() # Row sort
		        store = np.array(haplotypes)
		        tajD = tajD + [popgenStats(store)[0]]
		        wuH = wuH + [popgenStats(store)[1]]
		        store = store.reshape((1,) + store.shape + (1,))
		        scores = hapCNN.predict(store)
		        upper = upper + [np.quantile(scores, 0.975)]
		        mean = mean + [np.mean(scores)]
		        lower = lower + [np.quantile(scores, 0.025)]
		index1 = index1 + [i1]
		index2 = index2 + [i2]
	# Classification
	rnn_classify = np.transpose(np.array([mean, wuH]))
	rnn_classify = rnn_classify.reshape((1,) + rnn_classify.shape)
	rnn_classify = hapRNN_classification.predict(rnn_classify)[0][0]
	rnn_classification.append(rnn_classify)
	# Regressuib
	rnn_regressor = ensemblePredict(tajD = tajD, wuH = wuH, CNN = mean)
	rnn_regression.append(rnn_regressor)
	print(rnn_regressor)

# Compute RNN predictions
output = pd.DataFrame({'n_subsample': [i for i in range(NUMBER_SUBSAMPLES)], 'rnn_classification': rnn_classification, 'rnn_regression': rnn_regression})

output.to_csv(OUTPUT_DIRECTORY + 'empirical_predictions_' + TARGET_FILE.split('/')[-1].split('.csv')[0] + '.csv', index = False)