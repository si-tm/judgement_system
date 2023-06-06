#!/usr/bin/env python
# coding: utf-8

# In[29]:


import sys
import functools
from qdpy import algorithms, containers, benchmarks, plots
import pickle
import numpy as np
import importlib
import glob

# reference from takechan-san


# Define evaluation function
# arctan
def eval_from_svm(array):
    val = np.arctan(array)
    arrs = np.zeros(256)
    a = arrs.sum() #存在するストランド組み合わせの種類数。2~7までの値の列。
    b = val.std() #SVM100個で予測した値の標準偏差の列。予測値の信頼度を表す。
    return (np.average(val),), (a,b)


def run_qdpy(dirpath="test"):
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
    #コンテナとアルゴリズムを作成する。ここでは、進化によってGridコンテナを照らすことで、MAP-Elitesを使用しています。
    
    #評価結果を配置するgridを作成。
    grid = containers.Grid(
        shape=(64,64), 
        max_items_per_bin=1, 
        fitness_domain=((0.0, 1.),), #評価関数が返す値の範囲
        #評価関数はどんな関数でもいいが、返すデータ型は"fitness"と(特徴1,特徴2)という形
        features_domain=((1., 256.), (-2., 2.)))#軸。横strand(特徴１）数(1~7ぐらい)、縦（特徴2）はブースティングで求めたstd
    
    #配置アルゴリズムを指定。今回はエネルギーが小さいほど高評価なので、minimization。
    algo = algorithms.RandomSearchMutPolyBounded(
        grid, 
        budget=10000, 
        batch_size=500,
        dimension=256, #1つのストランドセットに幾つパラメータがあるか # one bit per strand
        optimisation_task="minimization")
    
    # Create a logger to pretty-print everything and generate output data files
    #すべてをプリティプリントするロガーを作成し、出力データファイルを生成する。
    #配置されたデータはpickleファイルから全て取得可能。
    logger = algorithms.AlgorithmLogger(algo)

    logger.final_filename = dirpath + "/qdpy_log.p"
    with open(dirpath + "/bootstrap_models.p","rb") as f:
        SVRs = pickle.load(f)
        
    # Run illumination process !
    #配置を実行する。
    #評価関数はfunctools.partialによって「svmとthが指定された、引数がarrayの新しいオブジェクト」になる。
    best = algo.optimise(functools.partial(eval_from_svm,svm=SVRs,th=0.95))
    #print(algo.summary())
    
    # Plot the results
    plots.default_plots_grid(logger)
    print("All results are available in the '%s' pickle file." % logger.final_filename)


# 2023/05/09<br>
# https://blog.amedama.jp/entry/2015/11/28/000432<br>
# functools.partialは、関数やメソッドの引数の一部をある値に固定した形で新しい呼び出し可能オブジェクトを作ることができる。partialを使って、eval_from_svmの新しいオブジェクトを作り、algo.optimiseに入力できるようにする必要がある。

def getData(target_file):
    target_file = sys.argv[1]
    data = {}
    with open(target_file, 'rb') as f:
        data = pickle.load(f)

    return data


def main():
    if len(sys.argv) != 2:
        print("usage : python3 tryQdpy.py [target pickle file]")
        exit()

    data = getData(target_file=sys.argv[1])
    run_qdpy()
    

if __name__ == "__main__":
    main()


