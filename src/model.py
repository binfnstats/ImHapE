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

def mergeData(positive_dir, neutral_dir, n = 100):
    '''Merges positive and neutral simulations into one array for CNN training'''
    import numpy as np
    import random
    import os
    pos_simulations = random.choices(os.listdir(positive_dir), k=n)
    neu_simulations = random.choices(os.listdir(neutral_dir), k=n)
    positive_haplotypes = np.vstack(np.array([np.load(positive_dir + i) for i in pos_simulations]))
    neutral_haplotypes = np.vstack(np.array([np.load(neutral_dir + i) for i in neu_simulations]))
    haplotypes = np.vstack((neutral_haplotypes, positive_haplotypes))
    modes = ['neutral' for i in range(len(neutral_haplotypes))] + ['positive' for i in range(len(positive_haplotypes))]
    return (haplotypes, modes)

def trainTestData(haplotypes, modes, p):
    '''Splits data into training and test set for CNN training'''
    import numpy as np
    import tensorflow as tf
    from sklearn.model_selection import train_test_split
    train_x, test_x, train_y, test_y = train_test_split(haplotypes, modes, stratify = modes, test_size = p, shuffle = True)
    train_dataset, test_dataset = (train_x, train_y), (test_x, test_y)
    return (train_dataset, test_dataset)

def plotTrainingData(train):
    '''Visualize some of the training data'''
    class_names = ['positive', 'neutral']
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10,10))
    train_x, train_y = train
    for i in range(25):
        plt.subplot(5,5,i+1)
        plt.xticks([])
        plt.yticks([])
        plt.grid(False)
        plt.imshow(train_x[i], cmap=plt.cm.binary)
        plt.xlabel(train_y[i])
    plt.show()

def trainCNN(train, test, nrows = 200, ncols = 2500, size_batch = 32):
    '''Train CNN on aligned haplotypes. Setup to train on 200 haplotype X 1000 base pair alignments'''
    import tensorflow as tf
    from tensorflow.keras import datasets, layers, models, regularizers
    import numpy as np
    # Load data
    train_x, train_y = train
    test_x, test_y = test
    # Some more preprocessing
    train_x = train_x.reshape(train_x.shape + (1,))
    test_x = test_x.reshape(test_x.shape + (1,))
    train_y = np.array([1 if i == 'positive' else 0 for i in train_y])
    test_y = np.array([1 if i == 'positive' else 0 for i in test_y])
    # Early stopping
    callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5)
    # Initialize model
    model = models.Sequential()
    model.add(layers.Conv2D(32, (3, 3), activation='relu', input_shape=(nrows, ncols, 1), kernel_regularizer=regularizers.l1_l2(l1=0.008, l2=0.008)))
    model.add(layers.MaxPooling2D((2, 2)))
    model.add(layers.Conv2D(64, (3, 3), activation='relu', kernel_regularizer=regularizers.l1_l2(l1=0.006, l2=0.006)))
    model.add(layers.MaxPooling2D((2, 2)))
    model.add(layers.Conv2D(64, (3, 3), activation='relu', kernel_regularizer=regularizers.l1_l2(l1=0.004, l2=0.004)))
    model.add(layers.Flatten())
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(1, activation='sigmoid'))
    # Compile and train model
    model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=['accuracy'])
    if size_batch > 0:
        history = model.fit(train_x, train_y, batch_size = size_batch, epochs=50, validation_data=(test_x, test_y), callbacks = [callback])
    else:
        history = model.fit(train_x, train_y, epochs=50, validation_data=(test_x, test_y), callbacks = [callback])
    # Evaluate performance
    test_loss, test_acc = model.evaluate(test_x, test_y, verbose=2)
    return (model, history, test_loss, test_acc)

def trainLSTM(features, targets):
    '''Train RNN on CNN sliding window output'''
    import tensorflow as tf
    from tensorflow.keras import datasets, layers, models, regularizers
    import numpy as np
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
    
    return (model, dt_add)