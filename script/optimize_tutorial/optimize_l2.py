#!/usr/bin/env python
# coding: utf-8

from qdpy import algorithms, containers, benchmarks, plots
from qdpy.base import ParallelismManager
import math
import pickle
import os
import numpy as np
from keras.models import load_model

import tensorflow as tf
print("TensorFlow version:", tf.__version__)
#tf.get_logger().setLevel('ERROR')


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

import functools


def getModel(path="../../saved_model/l1_ave_230530"):
    model = load_model(path)
    return model


def set_eval(ind, averageModel, deviationModel, scale=10.0):
    # X, Y = getXY()
    
    indexes = np.array(ind[:-1]) > ind[-1]
    strands = [1 if a else 0 for a in indexes]
    score = 2*math.atan(averageModel.predict([strands], verbose = 0)[0]/scale)/math.pi
    fit0 = ind[-1]#deviationModel.predict([strands])[0,0]
    fit1 = np.sum(indexes)
    features = (fit0, fit1)
    return (score,), features


def run_qdpy(dirpath="test"):
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
    
    #評価結果を配置するgridを作成。
    #L1に対して
    grid = containers.Grid(
        shape=(16,16), 
        max_items_per_bin=1, 
        fitness_domain=((0.0, 1.),), #評価関数が返す値の範囲 average 
        #評価関数はどんな関数でもいいが、返すデータ型は"fitness"と(特徴1,特徴2)という形 
        features_domain=((0., 1.), (1, 256))) #軸 deviation, number of strands 
    
    #配置アルゴリズムを指定。今回はエネルギーが小さいほど高評価なので、minimization。
    #??
    algo = algorithms.RandomSearchMutPolyBounded(
        grid, 
        # budget=10000, 
        budget=1000, 
        batch_size=100,
        # dimension=17, #1つのストランドセットに幾つパラメータがあるか # one bit per strand
        dimension=257, #1つのストランドセットに幾つパラメータがあるか # one bit per strand
        optimisation_task="maximization")
    
    # Create a logger to pretty-print everything and generate output data files
    #すべてをプリティプリントするロガーを作成し、出力データファイルを生成する。
    #配置されたデータはpickleファイルから全て取得可能。
    logger = algorithms.AlgorithmLogger(algo)

    # Run illumination process !
    #配置を実行する。
    averageModel = getModel('../../saved_model/l2_ave_230613')
    deviationModel = getModel('../../saved_model/l2_dev_230613')
    eval_fn = functools.partial(set_eval,averageModel=averageModel,deviationModel=deviationModel)
    best = algo.optimise(eval_fn)
    print(algo.summary())
    #print(type(algo))
    #print("wow")
    
    # Plot the results
    logger.final_filename = dirpath + "/qdpy_log_l1_230613.p"
    print(logger.final_filename)
    plots.default_plots_grid(logger)
    print("All results are available in the '%s' pickle file." % logger.final_filename)


def getXY():
    x_l2_mean_fullPath = os.path.abspath('../../data/npy/x_random_l2_6_mean.npy')
    y_l2_mean_fullPath = os.path.abspath('../../data/npy/y_random_l2_6_mean.npy')

    x_l2_mean_path = tf.keras.utils.get_file('x_random_l2_6_mean.npy', 'file://'+x_l2_mean_fullPath)
    y_l2_mean_path = tf.keras.utils.get_file('y_random_l2_6_mean.npy', 'file://'+y_l2_mean_fullPath)

    x_data = np.load(x_l2_mean_path)
    y_data = np.load(y_l2_mean_path, allow_pickle=True)
    return x_data, y_data


def main():
    run_qdpy()
    x_data, y_data = getXY()
    # 256
    print(len(x_data[0]))


    

if __name__ == "__main__":
    main()


