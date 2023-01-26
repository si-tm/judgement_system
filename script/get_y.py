import numpy as np
import matplotlib.pylab as plt
import load_results as lr
from measuring_volume import run_output_bonds as rob
from measuring_volume import convexhull_volume as cv
import pickle

# def get_y():
#     # y_dic[dir_name] = mean_val
#     y_dic = {}
#     dirs = lr.load_directory()
#     for d in dirs:
#         print(d)
#         # get bonds 
#         rob.make_bonds(d)
#         # convex_hull
#         volume = cv.convexhull_volume_all_strands(d)
#         y_dic[d] = volume
    
#     print(y_dic)
#     with open("../data/dic/y_1.pkl", "wb") as tf:
#         peickle.dump(y_dic,tf)

def get_y_l1():
    # y_dic[dir_name] = mean_val
    y_dic = {}
    dirs = lr.load_l1_directory()
    for d in dirs:
        print(d)
        # get bonds 
        rob.make_bonds(d)
        # convex_hull
        volume = cv.convexhull_volume_all_strands(d)
        y_dic[d] = volume
    
    print(y_dic)
    with open("../data/dic/y_l1_1.pkl", "wb") as tf:
        pickle.dump(y_dic,tf)

def get_y_l2():
    # y_dic[dir_name] = mean_val
    y_dic = {}
    dirs = lr.load_l2_directory()
    for d in dirs:
        print(d)
        # get bonds 
        rob.make_bonds(d)
        # convex_hull
        volume = cv.convexhull_volume_all_strands(d)
        y_dic[d] = volume
    
    print(y_dic)
    with open("../data/dic/y_l2_1.pkl", "wb") as tf:
        pickle.dump(y_dic,tf)

def get_y_l3():
    # y_dic[dir_name] = mean_val
    y_dic = {}
    dirs = lr.load_l3_directory()
    for d in dirs:
        print(d)
        # get bonds 
        rob.make_bonds(d)
        # convex_hull
        volume = cv.convexhull_volume_all_strands(d)
        y_dic[d] = volume
    
    print(y_dic)
    with open("../data/dic/y_l3_1.pkl", "wb") as tf:
        pickle.dump(y_dic,tf)

def get_y(result_path, dirs):
    y_dic = {}
    for d in dirs:
        rob.make_bonds(d)
        volume = cv.convexhull_volume_all_strands(d)
        y_dic[d] = volume
    
    with open(result_path, "wb") as tf:
        pickle.dump(y_dic,tf)
    

def load_y(path):
    with open(path, "rb") as tf:
        new_dict = pickle.load(tf)
        print(new_dict)
        return new_dict

def test():
    dirs_l1 = lr.load_random_dir("oxdna_random_3", "L1")
    y_l1_path = "../data/dic/y_random_l1_3.pkl"
    # dirs_l2 = lr.load_random_dir("oxdna_random_3","L2")
    # y_l2_path = "../data/dic/y_random_l2_2.pkl"
    # dirs_l3 = lr.load_random_dir("oxdna_random_3","L3")
    # y_l3_path = "../data/dic/y_random_l3_2.pkl"
    
    get_y(y_l1_path, dirs_l1)
    # get_y(y_l2_path, dirs_l2)
    # get_y(y_l3_path, dirs_l3)


def main():
    y_l1_path = "../data/dic/y_l1_1.pkl"
    y_l2_path = "../data/dic/y_l2_1.pkl"
    y_l3_path = "../data/dic/y_l3_1.pkl"
    # get_y_l1()
    # get_y_l2()
    # get_y_l3()
    load_y("../data/dic/y_1.pkl")
    # load_y(y_l1_path)
    # load_y(y_l2_path)
    # load_y(y_l3_path)

if __name__ == '__main__':
    # main()
    test()
  