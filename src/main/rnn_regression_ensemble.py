#!/bin/env python3

# Parameters from caret R 10-fold cross-validation
# n.trees = 100, depth = 2, shrinkage = 0.1, n.minobsinnode = 10

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_percentage_error
import pandas as pd 
import numpy as np

df = pd.read_csv('PATH_TO_RNN_REGRESSION_PREDICTIONS')

df = df.loc[:, ['0.0_tajD', '0.001_tajD_CNN', '0.002_CNN', '0.0_wuH_CNN', '0.003_CNN', '0.0_CNN', '0.0_multi', '0.0_tajD_CNN', '0.0_tajD_wuH', '0.001_tajD', 'fitness']]
X = df.loc[:,df.columns != "fitness"]
y = df.loc[:, "fitness"]

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size = 0.1)
ensembled = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=2, min_samples_leaf = 10, random_state=0).fit(train_X, train_y)

from joblib import dump
dump(ensembled, 'OUTPUT_PATH_TO_MODEL_FOLDER/rnn_ensemble.joblib')