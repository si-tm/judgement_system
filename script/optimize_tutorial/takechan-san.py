#!/usr/bin/env python
# coding: utf-8

# In[29]:


import sys
import functools
from qdpy import algorithms, containers, benchmarks, plots
import pickle
import numpy as np
import importlib
import read_learning_result as rlr
importlib.reload(rlr)


# In[2]:


# Define evaluation function
#QDPYの評価関数を定義する
#100個のモデルで予測したエネルギー平均値が小さいほど良いとする。
def eval_from_svm(array,svm,th):
    args = np.argsort(array)
    arrs = np.zeros(256)
    for i in args[-7:]:
        if array[i] > th:
            arrs[i] = 1.0#数字のリストIndividualをを0,1の値のリストに変換する。
    val = np.array([s.predict([arrs]) for s in svm])#評価対象。今回はSVM100個で予測したエネルギー値の平均値のリスト。
    #valはブートストラップされたモデルそれぞれについて予測を行った結果のリスト。
    a = arrs.sum() #存在するストランド組み合わせの種類数。2~7までの値の列。
    b = val.std() #SVM100個で予測した値の標準偏差の列。予測値の信頼度を表す。
    return (np.average(val),), (a,b)
#返す値は予測エネルギー平均値と、a,b


#bのうち、上位N%に対応するものをoxDNAに渡すようにしたい


# In[1]:


def run_qdpy(dirpath):
    # Create container and algorithm. Here we use MAP-Elites, by illuminating a Grid container by evolution.
    #コンテナとアルゴリズムを作成する。ここでは、進化によってGridコンテナを照らすことで、MAP-Elitesを使用しています。
    
    #評価結果を配置するgridを作成。
    grid = containers.AutoScalingGrid(
        shape=(64,64), 
        max_items_per_bin=1, 
        fitness_domain=((- np.inf, 1.),), #評価関数が返す値の範囲
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

# In[31]:


def main(args):
    dirpath = "/".join(args[1:])
    run_qdpy(dirpath)


# In[ ]:


if __name__ == "__main__":
    args = sys.argv
    #main()
    sys.exit(main(sys.argv))



