#!/usr/bin/env python3
''' Description: This script takes CNN estimates from sliding window analysis of full 29903 bp simulated genomes and trains an LSTM to predict if the window actually contains a beneficial mutation '''

import tensorflow as tf
import pandas as pd
import numpy as np
import os
from tensorflow.keras import datasets, layers, models, regularizers
from sklearn.metrics import confusion_matrix
import random

random.seed(123456)

RNN_PERFORMANCE_OUTPUT = 'NAME_OF_FILE_FOR_EVALUATING_RNN_PERFORMANCE.csv'
RNN_MODEL_OUTPUT = 'RNN_LSTM_TENSORFLOW_MODEL_NAME.tf'

# Extract CNN predictions of simulated genomes using 2500 bp sliding windows
# -----------------------------------------------------------------------------------
# *NOTE: Must be run from directory containing CNN estimates across simulated genomes

files = [i for i in os.listdir() if '2500' in i]
random.shuffle(files)
ids = [i.split('_')[-1].split('.')[0] for i in files]
features = []
targets = []
vafs = []
for id in ids:
	print(id)
	check = [i for i in files if id in i and '_1.0_' in i] # Remove samples without any beneficial loci
	if len(check) > 0:
		continue
	preds = pd.read_csv([i for i in files if id in i and 'predictions' in i and '2500' in i][0])
	loci = pd.read_csv([i for i in files if id in i and 'loci' in i and '2500' in i][0])
	dt_loci = loci[loci['mode'] == 'Beneficial'][loci['VAF'] > 0.1]
	preds['selection'] = 0
	inds = [0 for i in range(29903)]
	for i in range(preds.shape[0]):
		bens = dt_loci['loci'].tolist()
		for j in bens:
			preds.loc[(j > preds.index1) & (j < preds.index2), 'selection'] = 1
	features.append(preds['mean'].tolist())
	targets.append(preds['selection'].tolist()) 

# Train an LSTM on the sliding window output 80% training / 20% test
# ------------------------------------------------------------------
features2 = np.array(features)
features2 = features2.reshape(features2.shape + (1,))
targets2 = np.array(targets)
targets2 = targets2.reshape(targets2.shape + (1,))
train_x = features2[0:int(0.8*len(features2))]
train_y = targets2[0:int(0.8*len(targets2))]
test_x = features2[int(0.8*len(features2)):int(len(features2))]
test_y = targets2[int(0.8*len(targets2)):int(len(targets2))]
model = tf.keras.Sequential()
model.add(layers.LSTM(275, return_sequences=True))
model.add(tf.keras.layers.Dense(units=1))
model.add(layers.Activation('sigmoid'))
model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=['accuracy'])
model.fit(train_x, train_y, epochs=25, batch_size = 16, validation_data=(test_x, test_y))

# Evaluate predictions on test set
# --------------------------------
checking = model.predict(test_x)
checking = checking.reshape(checking.shape[0:2])
test_y2 = test_y
test_y2 = test_y2.reshape(test_y2.shape[0:2])
dt_add = pd.DataFrame()
for i in range(checking.shape[0]):
	pred1 = checking[i,:].tolist()
	actu1 = test_y2[i,:].tolist()
	add = pd.DataFrame({'pred':pred1, 'actual':actu1})
	dt_add = pd.concat([dt_add, add], ignore_index=True)

dt_add.to_csv(RNN_PERFORMANCE_OUTPUT)

# Run the RNN on 9 example simulated genomes 
# -------------------------------------------
count = 1
for id in ids[int(0.8*len(features2)):int(len(features2))]:
	count += 1
	print(count)
	check = [i for i in files if id in i and '_1.0_' in i]
	if len(check) > 0:
		continue
	preds = pd.read_csv([i for i in files if id in i and 'predictions' in i and '2500' in i][0])
	loci = pd.read_csv([i for i in files if id in i and 'loci' in i and '2500' in i][0])
	evaluate_mean = np.array(preds['mean'].tolist())
	evaluate_mean = evaluate_mean.reshape((1,) + evaluate_mean.shape + (1,))
	preds['RNN'] = model.predict(evaluate_mean)[0]
	loci.to_csv('RNN_' + [i for i in files if id in i and 'loci' in i and '2500' in i][0])
	preds.to_csv('RNN_' + [i for i in files if id in i and 'predictions' in i and '2500' in i][0])

tf.keras.models.save_model(model, filepath = RNN_MODEL_OUTPUT)




