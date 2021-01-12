#!/usr/bin/env python3
''' Description: This script contains functions relevant to sampling simulated haplotypes (sampleData), merging neutral and positive simulations (mergeData), and training CNN (trainTestData and trainCNN) '''

def sampleData(sims, size, reps, sort_row = False, sort_col = False):
    '''Generates a haplotype array (predictor) and evolutionary mode and fitness (responses)'''
    import numpy as np
    import tensorflow as tf
    store_mode = []
    store_haplotypes = []
    fitness = sims.w
    mode = str(fitness) if fitness > 1 else 'neutral'
    for i in range(reps):
        haplotypes = sims.mutations
        indices = np.random.choice(haplotypes.shape[0], size = size, replace = False)
        samples = haplotypes[indices,:]
        haplotypes = [[0 if samples[i,j] == sims.reference[j] else samples[i,j] for j in range(samples.shape[1])] for i in range(samples.shape[0])]
        haplotypes = (np.array(haplotypes) > 0).astype(int).tolist()
        if sort_row == True:
            haplotypes = np.array(haplotypes)
            haplotypes = haplotypes[np.argsort(haplotypes.sum(axis=1))[::-1],:].tolist()
        if sort_col == True:
            haplotypes = np.array(haplotypes)
            haplotypes = haplotypes[:,np.argsort(haplotypes.sum(axis=0))[::-1]].tolist()
        store_mode.append(mode)
        store_haplotypes.append(haplotypes)
    return (np.array(store_haplotypes), np.array(store_mode))

def mergeData(positive_dir, neutral_dir, n = 100, remove = 'NAN'):
    '''Merges positive and neutral simulations into one array for CNN training'''
    import numpy as np
    import random
    import os
    pos_simulations = random.sample(os.listdir(positive_dir), k=n)
    neu_simulations = random.sample(os.listdir(neutral_dir), k=n)
    positive_haplotypes = np.vstack(np.array([np.load(positive_dir + i) for i in pos_simulations]))
    neutral_haplotypes = np.vstack(np.array([np.load(neutral_dir + i) for i in neu_simulations]))
    haplotypes = np.vstack((neutral_haplotypes, positive_haplotypes))
    modes = ['neutral' for i in range(len(neutral_haplotypes))] + ['positive' for i in range(len(positive_haplotypes))]
    return (haplotypes, modes)

def trainTestData(haplotypes, modes, p, sort_row = True, sort_col = True):
    '''Splits data into training and test set for CNN training'''
    import numpy as np
    import tensorflow as tf
    from sklearn.model_selection import train_test_split
    tf.random.set_seed(123456)
    if sort_row == True:
        haplotypes = np.array(haplotypes)
        haplotypes = np.array([i[np.argsort(i.sum(axis=1))[::-1],:].tolist() for i in haplotypes])
    if sort_col == True:
        haplotypes = np.array(haplotypes)
        haplotypes = np.array([i[:,np.argsort(i.sum(axis=0))[::-1]].tolist() for i in haplotypes])
    train_x, test_x, train_y, test_y = train_test_split(haplotypes, modes, stratify = modes, test_size = p, shuffle = True)
    train_dataset, test_dataset = (train_x, train_y), (test_x, test_y)
    return (train_dataset, test_dataset)

def trainCNN(train, test, nrows = 200, ncols = 2500, size_batch = 32):
    '''Train CNN on aligned haplotypes. Setup to train on nrows (haplotypes) X ncols (base pair) alignments'''
    import tensorflow as tf
    from tensorflow.keras import datasets, layers, models, regularizers
    import numpy as np
    tf.random.set_seed(123456)
    # Load data
    train_x, train_y = train
    test_x, test_y = test
    # Some more preprocessing
    train_x = train_x.reshape(train_x.shape + (1,))
    test_x = test_x.reshape(test_x.shape + (1,))
    train_y = np.array([1 if i == 'positive' else 0 for i in train_y])
    test_y = np.array([1 if i == 'positive' else 0 for i in test_y])
    # Early stopping
    callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=3)
    # Initialize model
    model = models.Sequential()
    model.add(layers.Conv2D(32, (3, 3), activation='relu', input_shape=(nrows, ncols, 1), kernel_regularizer=regularizers.l1_l2(l1=0.008, l2=0.008)))
    model.add(layers.MaxPooling2D((2, 2)))
    model.add(layers.Conv2D(64, (3, 3), activation='relu', kernel_regularizer=regularizers.l1_l2(l1=0.006, l2=0.006)))
    model.add(layers.MaxPooling2D((2, 2)))
    model.add(layers.Conv2D(64, (3, 3), activation='relu', kernel_regularizer=regularizers.l1_l2(l1=0.002, l2=0.002)))
    model.add(layers.Flatten())
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(1, activation='sigmoid'))
    METRICS = [tf.keras.metrics.BinaryAccuracy(name='accuracy'),tf.keras.metrics.Precision(name='precision'),tf.keras.metrics.Recall(name='recall'),tf.keras.metrics.AUC(name='auc')]
    model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=METRICS)
    if size_batch > 0:
        history = model.fit(train_x, train_y, batch_size = size_batch, epochs=30, validation_data=(test_x, test_y), callbacks = [callback])
    else:
        history = model.fit(train_x, train_y, epochs=25, validation_data=(test_x, test_y), callbacks = [callback])
    return (model, history)

def trainLSTMreg(features, fitness, multi = False, nepochs = 30, regularization):
    import tensorflow as tf
    import pandas as pd
    import numpy as np
    import os
    from tensorflow.keras import datasets, layers, models, regularizers
    from sklearn.model_selection import train_test_split
    import random
    tf.random.set_seed(1234567)
    train_x, test_x, train_y, test_y = train_test_split(features, fitness, stratify = fitness, test_size = 0.1, shuffle = True)
    train_y, test_y = train_y.reshape(train_y.shape + (1,) + (1,)), test_y.reshape(test_y.shape + (1,) +(1,))
    if multi == False:
        train_x, test_x = train_x.reshape(train_x.shape + (1,)), test_x.reshape(test_x.shape + (1,))
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Bidirectional(layers.LSTM(549, kernel_regularizer=regularizers.l1_l2(l1=regularization, l2=regularization))))
    model.add(tf.keras.layers.Dense(units=1))
    METRICS = [tf.keras.metrics.MeanSquaredError(name='mean_squared_error', dtype=None), tf.keras.metrics.MeanAbsolutePercentageError(name='mean_absolute_percentage_error', dtype=None)]
    model.compile(loss='mean_squared_error', optimizer='adam', metrics = METRICS)
    callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=4)
    model.fit(train_x, train_y, epochs=nepochs, batch_size = 32, validation_data=(test_x, test_y), callbacks = callback)
    return (model)

def trainLSTM(features, fitness, multi = False, nepochs = 30, regularization):
    import tensorflow as tf
    import pandas as pd
    import numpy as np
    import os
    from tensorflow.keras import datasets, layers, models, regularizers
    from sklearn.model_selection import train_test_split
    import random
    tf.random.set_seed(1234567)
    train_x, test_x, train_y, test_y = train_test_split(features, fitness, stratify = fitness, test_size = 0.1, shuffle = True)
    train_y = np.array([0 if i in np.where(np.array(train_y) == 1)[0] else 1 for i in range(len(train_y))])
    test_y = np.array([0 if i in np.where(np.array(test_y) == 1)[0] else 1 for i in range(len(test_y))])
    train_y, test_y = train_y.reshape(train_y.shape + (1,) + (1,)), test_y.reshape(test_y.shape + (1,) +(1,))
    if multi == False:
        train_x, test_x = train_x.reshape(train_x.shape + (1,)), test_x.reshape(test_x.shape + (1,))
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Bidirectional(layers.LSTM(549, kernel_regularizer=regularizers.l1_l2(l1=regularization, l2=regularization))))
    model.add(tf.keras.layers.Dense(units=1))
    model.add(layers.Activation('sigmoid'))
    METRICS = [tf.keras.metrics.BinaryAccuracy(name='accuracy'), tf.keras.metrics.Precision(name='precision'), tf.keras.metrics.Recall(name='recall'), tf.keras.metrics.AUC(name='auc')]
    model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=METRICS)
    callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=4)
    model.fit(train_x, train_y, epochs=nepochs, batch_size = 32, validation_data=(test_x, test_y), callbacks = callback)
    return (model)

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