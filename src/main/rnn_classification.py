#!/usr/bin/env python3
''' Description: This script takes CNN estimates from sliding window analysis of full 29903 bp simulated genomes and trains an LSTM to predict if the window actually contains a beneficial mutation '''

import tensorflow as tf
import pandas as pd
import numpy as np
import os
from tensorflow.keras import datasets, layers, models, regularizers
from sklearn.model_selection import train_test_split
import random
import argparse
from collections import Counter
import model as mod

# Arguments
parser = argparse.ArgumentParser(description = 'Simulation parameters')
parser.add_argument('-f', '--feature', help='Which feature to use for training.')
parser.add_argument('-r', '--regularization', help='L1/L2 regularization')

# Set parameters
args = parser.parse_args()
FT = args.feature
REG = float(args.regularization)

INPUT_DIR = 'PATH_TO_SLIDING_WINDOW_FEATURES'
OUTPUT_MODEL = 'OUTPUT_PATH_FOR_RNN_CLASSIFICATION_' + str(REG) + '_' + FT + '_rnn.tf'

# Load inputs
# -----------
fts = np.load(INPUT_DIR + FT + '.npy')
fitness = np.load(INPUT_DIR + 'fitness.npy')

# Train an LSTM on the sliding window output 80% training / 20% test
# ------------------------------------------------------------------
random.seed(645621)

if FT in ['tajD_wuH', 'tajD_CNN', 'wuH_CNN', 'multi', 'CNN_all']:
    m = mod.trainLSTM(features = fts, fitness = fitness, multi = True, regularization = REG)
else:
    m = mod.trainLSTM(features = fts, fitness = fitness, regularization = REG)

tf.keras.models.save_model(m, filepath = OUTPUT_MODEL)

