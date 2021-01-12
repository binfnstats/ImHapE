#!/usr/bin/env python3
''' Description: This script computes performance of RNN across independent validation simulations'''

import tensorflow as tf
import pandas as pd
import numpy as np
import datetime
import argparse
import os
import re
import datetime
from itertools import combinations

MODELS_DIRECTORY = '/MODELS_PATH/simmodels'
PROCESSED_DIRECTORY = '/TEST_SET_SIMULATIONS/simprocess_test_set'

os.chdir(MODELS_DIRECTORY)

rnns = []
for mod in os.listdir():
	if 'classification' in mod:
		vars()[mod] = tf.keras.models.load_model(mod)
		rnns.append(mod)

os.chdir(PROCESSED_DIRECTORY)

CNN_all = np.load('CNN_all.npy')
tajD = np.load('tajD.npy')
wuH = np.load('wuH.npy')
CNN = np.load('CNN.npy')
tajD_wuH = np.load('tajD_wuH.npy')
tajD_CNN = np.load('tajD_CNN.npy')
wuH_CNN = np.load('wuH_CNN.npy')
multi = np.load('multi.npy')
fitness = np.load('fitness.npy')
#mutrate = np.load('mutrate.npy')
#probben = np.load('probben.npy')

dt_classification = pd.DataFrame({'fitness':fitness})
#dt_regression = pd.DataFrame({'fitness':fitness})
for mod in rnns:
	print(mod)
	if 'classification' in mod:
		l1_l2 = mod.split('_')[1]		
		fts = re.split('classification_.*?_', mod)[1].split('_rnn.tf')[0]
		identifier =  l1_l2 + '_' + fts
		try:
			pred = vars()[mod].predict(np.array(vars()[fts]))
		except:
			fts = np.array(vars()[fts])
			fts = fts.reshape(fts.shape + (1,))
			pred = vars()[mod].predict(fts)
		dt_classification[identifier] = pred
	else:
		continue
	# if 'regression' in mod:
	# 	l1_l2 = mod.split('_')[1]
	# 	fts = re.split('regression_.*?_', mod)[1].split('_rnn.tf')[0]
	# 	identifier =  l1_l2 + '_' + fts
	# 	try:
	# 		pred = vars()[mod].predict(np.array(vars()[fts]))
	# 	except:
	# 		fts = np.array(vars()[fts])
	# 		fts = fts.reshape(fts.shape + (1,))
	# 		pred = vars()[mod].predict(fts)
	# 	dt_regression[identifier] = pred

dt_classification.to_csv('rnn_classification_predictions.csv', index = False)
#dt_regression.to_csv('/.mounts/labs/awadallalab/private/touellette/covid-19/analysis/rnn_regression_predictions.csv', index = False)