import numpy as np
import matplotlib.pylab as plt
import load_results as lr
from measuring_volume import run_output_bonds as rob
from measuring_volume import convexhull_volume as cv
import pickle

def get_y():
    # y_dic[dir_name] = mean_val
    y_dic = {}
    dirs = lr.load_directory()
    for d in dirs:
        print(d)
        # get bonds 
        rob.make_bonds(d)
        # convex_hull
        volume = cv.convexhull_volume_all_strands(d)
        y_dic[d] = volume
    
    print(y_dic)
    with open("../data/dic/y_1.pkl", "wb") as tf:
        pickle.dump(y_dic,tf)

def load_y():
    with open("../data/dic/y_1.pkl", "rb") as tf:
        new_dict = pickle.load(tf)

    print(new_dict)

def main():
    get_y()
    # load_y()

if __name__ == '__main__':
  main()
  