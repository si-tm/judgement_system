import numpy as np
import matplotlib.pylab as plt
import load_results as lr
import get_x as gx
import get_y as gy
from classify_seq import make_input_seq as mis

def dic2lst(dic, lst):
    param = []
    for l in lst:
        param.append(dic[l])
    return param

def make_xy():
    dirs = lr.load_directory()
    dic_x1 = gx.load_x1()
    dic_x2 = gx.load_x1()
    dic_y = gy.load_y()
    x1_lst = mis.seq_lst1()
    x2_lst = mis.seq_lst2()

    x1 = []
    x2 = []
    y = []

    for d in dirs:
        x1.append(dic2lst(dic_x1[d], x1_lst))
        x2.append(dic2lst(dic_x2[d], x2_lst))
        y.append(dic_y[d])

    np.save('../data/npy/x1.npy',x1)
    np.save('../data/npy/x2.npy',x2)
    np.save('../data/npy/y.npy',y)

def make_xy_l1():
    dirs_l1 = lr.load_l1_directory()
    dic_x_l1 = gx.load_x("../data/dic/x_l1_1.pkl")
    dic_y_l1 = gy.load_y("../data/dic/y_l1_1.pkl")
    x_l1_lst = mis.seq_lst("../input/input_seq_L1.csv")

    x_l1 = []
    y_l1 = []

    for d in dirs_l1:
        x_l1.append(dic2lst(dic_x_l1[d], x_l1_lst))
        y_l1.append(dic_y_l1[d])
    print(y_l1)
    
    np.save('../data/npy/x_l1_1.npy',x_l1)
    np.save('../data/npy/y_l1_1.npy',y_l1)

def load_npy(path):
    z = np.load(path) 
    return z

def main():
    # make_xy_l1()
    print(load_npy("../data/npy/y_l1_1.npy"))

if __name__ == '__main__':
  main()
