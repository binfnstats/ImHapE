#!/usr/bin/env python3
''' Description: This script takes CNN estimates from sliding window analysis of full 29903 bp simulated genomes and trains an LSTM to predict if the window actually contains a beneficial mutation '''

import tensorflow as tf
import pandas as pd
import numpy as np
import os
from tensorflow.keras import datasets, layers, models, regularizers
from sklearn.metrics import confusion_matrix
import random

OUTPUT = '/OUTPUT_PATH/simprocess/'
SET_DIR = '/SIMULATED_GENOME_PATH/row/'

# Extract CNN predictions of simulated genomes using 2500 bp sliding windows
# -----------------------------------------------------------------------------------
# *NOTE: Must be run from directory containing CNN estimates across simulated genomes
os.chdir(SET_DIR)
files = [i for i in os.listdir()]
ids = [i.split('_')[7].split('.')[0] for i in files]
ids = set([x for x in ids if ids.count(x) > 1])
random.shuffle(files)
tajD, wuH, CNN, CNN_all, tajD_wuH, tajD_CNN, wuH_CNN, multi, fitness = [], [], [], [], [], [], [], [], []
targets = []
count = 1
for id_ in ids:
    print(count)
    count += 1
    w = float([i for i in files if id_ in i and 'predictions' in i][0].split('_')[3])
    preds = pd.read_csv([i for i in files if id_ in i and 'predictions' in i][0])
    loci = pd.read_csv([i for i in files if id_ in i and 'loci' in i][0])
    dt_loci = loci[loci['mode'] == 'Beneficial'][loci['VAF'] > 0.05]
    if w > 1 and dt_loci.shape[0] == 0:
        continue # Only keep positive samples with at least one beneficial mutation above a VAF of 0.05 
    #preds['selection'] = 0
    # inds = [0 for i in range(29903)]
    # for i in range(preds.shape[0]):
    #     if dt_loci.shape[0] == 0:
    #         break
    #     bens = dt_loci['loci'].tolist()
    #     for j in bens:
    #         preds.loc[(j > preds.index1) & (j < preds.index2), 'selection'] = 1
    fitness.append(w)
    tajD.append(preds["Tajima's D"].tolist())
    wuH.append(preds["Fay and Wu's H"].tolist())
    CNN.append(preds["mean"].tolist())
    CNN_all.append(np.transpose(np.array([preds["upper"].tolist(), preds["mean"].tolist(), preds["lower"].tolist()])))
    tajD_wuH.append(np.transpose(np.array([preds["Fay and Wu's H"].tolist(),preds["Tajima's D"].tolist()])))
    tajD_CNN.append(np.transpose(np.array([preds['mean'].tolist(), preds["Tajima's D"].tolist()])))
    wuH_CNN.append(np.transpose(np.array([preds['mean'].tolist(), preds["Fay and Wu's H"].tolist()])))
    multi.append(np.transpose(np.array([preds['mean'].tolist(), preds["Fay and Wu's H"].tolist(),preds["Tajima's D"].tolist()])))
    #targets.append(preds['selection'].tolist())

# Save outputs
np.save(OUTPUT + 'CNN_all', CNN_all)
np.save(OUTPUT + 'tajD', tajD)
np.save(OUTPUT + 'wuH', wuH)
np.save(OUTPUT + 'CNN', CNN)
np.save(OUTPUT + 'tajD_wuH', tajD_wuH)
np.save(OUTPUT + 'tajD_CNN', tajD_CNN)
np.save(OUTPUT + 'wuH_CNN', wuH_CNN)
np.save(OUTPUT + 'multi', multi)
np.save(OUTPUT + 'fitness', fitness)
#np.save(OUTPUT + 'targets', targets) # Sequence target