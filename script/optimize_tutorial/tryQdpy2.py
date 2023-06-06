from qdpy import algorithms, containers, benchmarks, plots
from qdpy.base import ParallelismManager
import math
import pickle
import os
import numpy as np
from keras.models import load_model

import tensorflow as tf
print("TensorFlow version:", tf.__version__)

from tensorflow.keras.layers import Dense
from tensorflow.keras import Model
from tensorflow.keras import layers
from tensorflow import keras
import sys
import matplotlib.pyplot as plt
import pandas as pd
import csv
import pickle
from scipy.stats import linregress
from sklearn.model_selection import train_test_split
from sklearn.model_selection import train_test_split
from datetime import datetime
import numpy as np
from sklearn.preprocessing import Normalizer


def eval_fn2(ind):
    """An example evaluation function. It takes an individual as input, and returns the pair ``(fitness, features)``, where ``fitness`` and ``features`` are sequences of scores."""

    normalization = sum((x for x in ind))
    k = 10.
    score = 1. - sum(( math.cos(k * ind[i]) * math.exp(-(ind[i]*ind[i])/2.) for i in range(len(ind)))) / float(len(ind))
    fit0 = sum((x * math.sin(abs(x) * 2. * math.pi) for x in ind)) / normalization
    fit1 = sum((x * math.cos(abs(x) * 2. * math.pi) for x in ind)) / normalization
    features = (fit0, fit1)
    return (score,), features


def runQDpy1():
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
    grid = containers.Grid(
        shape=(64,64), 
        max_items_per_bin=1, 
        fitness_domain=((0., 1.),), 
        features_domain=((0., 1.), (0., 1.)))
    algo = algorithms.RandomSearchMutPolyBounded(
        grid, 
        budget=60000, 
        batch_size=500,
        dimension=3, 
        optimisation_task="maximisation")

    # Create a logger to pretty-print everything and generate output data files
    logger = algorithms.AlgorithmLogger(algo)

    # Define evaluation function
    eval_fn = algorithms.partial(
        benchmarks.illumination_rastrigin_normalised,
        nb_features = len(grid.shape))

    # Run illumination process !
    best = algo.optimise(eval_fn)

    # Print results info
    print(algo.summary())

    # Plot the results
    plots.default_plots_grid(logger)

    print("All results are available in the '%s' pickle file." % logger.final_filename)

def runQDpy2():
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
    grid = containers.Grid(
        shape=(16,16), 
        max_items_per_bin=1, 
        fitness_domain=((-math.pi, math.pi),), 
        features_domain=((0., 1.), (0., 1.)))
    algo = algorithms.RandomSearchMutPolyBounded(
        grid, 
        budget=3000, 
        batch_size=500,
        dimension=3, 
        optimisation_task="minimisation")

    # Create a logger to pretty-print everything and generate output data files
    logger = algorithms.TQDMAlgorithmLogger(algo)

    # Run illumination process !
    with ParallelismManager("none") as pMgr:
        best = algo.optimise(eval_fn2, executor = pMgr.executor, batch_mode=False) # Disable batch_mode (steady-state mode) to ask/tell new individuals without waiting the completion of each batch

    # Print results info
    print("\n" + algo.summary())

    # Plot the results
    plots.default_plots_grid(logger)

    print("\nAll results are available in the '%s' pickle file." % logger.final_filename)
    print(f"""
To open it, you can use the following python code:
    import pickle
    # You may want to import your own packages if the pickle file contains custom objects
    with open("{logger.final_filename}", "rb") as f:
        data = pickle.load(f)
    # ``data`` is now a dictionary containing all results, including the final container, all solutions, the algorithm parameters, etc.
    grid = data['container']
    print(grid.best)
    print(grid.best.fitness)
    print(grid.best.features)
    """)

def main():
    # prepare input
    target = "../../data/dic/L1_data_4.pkl"
    
    dic = {}
    with open(target, "rb") as t:
        dic = pickle.load(t)
        # print(dic)

    model = load_model('../../saved_model/l1_ave_230530')
    x_l1_mean_fullPath = os.path.abspath('../../data/npy/x_random_l1_6_mean.npy')
    y_l1_mean_fullPath = os.path.abspath('../../data/npy/y_random_l1_6_mean.npy')

    x_l1_mean_path = tf.keras.utils.get_file('x_random_l1_6_mean.npy', 'file://'+x_l1_mean_fullPath)
    y_l1_mean_path = tf.keras.utils.get_file('y_random_l1_6_mean.npy', 'file://'+y_l1_mean_fullPath)

    x_data = np.load(x_l1_mean_path)
    y_data = np.load(y_l1_mean_path, allow_pickle=True)

    print(x_data[0])
    print(model.predict(x_data[0:1]))

    runQDpy1()
    runQDpy2()



if __name__ == '__main__':
    main()
