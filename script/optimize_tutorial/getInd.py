import pickle
import numpy as np
from qdpy import algorithms, containers, benchmarks, plots
from qdpy.base import *

# final.pからスコアの高いindを取り出す

def Ind2strand(ind):
    indexes = np.array(ind[:-1]) > ind[-1]
    strands = [1 if a else 0 for a in indexes]
    return strands


# https://gitlab.com/leo.cazenille/qdpy/-/blob/master/qdpy/containers.py
def readFinal(path):
    with open(path, "rb") as f:
        data = pickle.load(f)
    # ``data`` is now a dictionary containing all results, including the final container, all solutions, the algorithm parameters, etc.
    grid = data['container']
    # print(grid)
    for ind in grid.best:
        if ind > 0.9:
            pass
    
    lst = Ind2strand(grid.best)
    print(lst)
    print(grid.best.fitness)
    print(grid.best.features)

def test():
    path="optimizationresults_20230620023641/final.p"
    readFinal(path)

if __name__ == '__main__':
    test()
    pass
