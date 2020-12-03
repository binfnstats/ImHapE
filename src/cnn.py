#!/usr/bin/env python
''' Description: This script trains the deep learning (CNN) model on haplotype data from simulations '''

import tensorflow as tf
import pandas as pd
import numpy as np
import datetime
import argparse
import model as mod

# Arguments
parser = argparse.ArgumentParser(description = 'Parameters for running CNN-based analysis')
parser.add_argument('-pos', '--positive', help='Directory storing positive simulated data')
parser.add_argument('-neu', '--neutral', help='Directory storing neutral simulated data')
parser.add_argument('-p', '--prop', help='Percentage of data assigned to test set')
parser.add_argument('-o', '--out', help='Output directory ending with /')
parser.add_argument('-n', '--num', help='Number of training samples per evolutionary class')
parser.add_argument('-w', '--fitness', type = str, default = 'none', help='Fitness of positive simulations used to maintain standardized file labels.')
parser.add_argument('-nr', '--nrows', type = int, default = 200, help='Number of aligned haplotypes.')
parser.add_argument('-nc', '--ncols', type = int, default = 1000, help='Number of base pairs in alignment.')

args = parser.parse_args()
pos_dir = str(args.positive)
neu_dir = str(args.neutral)
out = str(args.out)
p = float(args.prop)
n = int(args.num)
name = args.fitness
nrows = args.nrows
ncols = args.ncols

# Build training and test sets
print('[1] Merging haplotype data \n')
haplotypes, modes = mod.mergeData(positive_dir = pos_dir, neutral_dir = neu_dir, n = n)

# Generate training and test tensors
print('[2] Generating training and testing sets \n')
train, test = mod.trainTestData(haplotypes = haplotypes, modes = modes, p = p)

haplotypes = None
modes = None

# Run CNN
print('[3] Running CNN \n')
model, history, test_loss, test_acc = mod.trainCNN(train = train, test = test, nrows = nrows, ncols = ncols)

# Save history
instant = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M')

pd.DataFrame(history.history).to_csv(out + name + '_' + instant + '_' + 'history.csv')
tf.keras.models.save_model(model, filepath = out + instant + '_' + 'cnn.tf')