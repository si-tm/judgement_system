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

def make_xy(dirs, target_name, seq_csv):
    dic_x = gx.load_x("../data/dic/x_" + target_name + ".pkl")
    dic_y = gy.load_y("../data/dic/y_" + target_name + ".pkl")
    x_lst = mis.seq_lst(seq_csv)

    x = []
    y = []

    print(len(dic_x))
    print(len(dic_y))

    for d in dirs:
        # print(d)
        x.append(dic2lst(dic_x[d], x_lst))
        y.append(dic_y[d])
    
    # print(x)
    # print(y)
    
    np.save("../data/npy/x_" + target_name + ".npy",x)
    np.save("../data/npy/y_" + target_name + ".npy",y)

def make_xy_ratioB(dirs, target_name, seq_csv):
    dic_x = gx.load_x("../data/dic/x_random_l1_6.pkl")
    dic_y = gy.load_y("../data/dic/L1_data_1.pkl")
    x_lst = mis.seq_lst(seq_csv)

    x = []
    y = []

    print(len(dic_x))
    print(len(dic_y))

    for d in dirs:
        # print(d)
        x.append(dic2lst(dic_x[d], x_lst))
        y.append(dic_y[d]["ratio_of_volume"])
    
    # print(x)
    # print(y)
    
    np.save("../data/npy/x_" + target_name + "_rb.npy",x)
    np.save("../data/npy/y_" + target_name + "_rb.npy",y)

def make_xy_ratioB(dirs, target_name, seq_csv):
    dic_x = gx.load_x("../data/dic/x_random_l2_6.pkl")
    dic_y = gy.load_y("../data/dic/L2_data_1.pkl")
    x_lst = mis.seq_lst(seq_csv)

    x = []
    y = []

    print(len(dic_x))
    print(len(dic_y))

    for d in dirs:
        # print(d)
        x.append(dic2lst(dic_x[d], x_lst))
        y.append(dic_y[d]["deviation_of_volume"])
    
    # print(x)
    # print(y)
    
    np.save("../data/npy/x_" + target_name + "_dv.npy",x)
    np.save("../data/npy/y_" + target_name + "_dv.npy",y)



def load_npy(path):
    z = np.load(path) 
    return z

def test():
    # dirs_random_l1 =  lr.load_random_dir("L1")
    # seq_csv_path = "../input/input_seq_L1.csv"
    # make_xy(dirs_random_l1, "random_l1_1", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l1_1.npy"))

    # dirs_random_l2 =  lr.load_random_dir("L2")
    # seq_csv_path = "../input/input_seq_L2.csv"
    # make_xy(dirs_random_l2, "random_l2_1", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l2_1.npy"))

    # dirs_random_l3 =  lr.load_random_dir("L3")
    # seq_csv_path = "../input/input_seq_L3.csv"
    # make_xy(dirs_random_l3, "random_l3_1", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l3_1.npy"))
    # dirs_random_l1 =  lr.load_random_dir("L1")

    # dirs_random_l1 =  lr.load_random_dir("L1")
    # seq_csv_path = "../input/input_seq_L1.csv"
    # make_xy(dirs_random_l1, "random_l1_2", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l1_2.npy"))
    # print(load_npy("../data/npy/y_random_l1_2.npy"))

    # dirs_random_l2 =  lr.load_random_dir("L2")
    # seq_csv_path = "../input/input_seq_L2.csv"
    # make_xy(dirs_random_l2, "random_l2_2", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l2_2.npy"))
    # print(load_npy("../data/npy/y_random_l2_2.npy"))

    # dirs_random_l3 =  lr.load_random_dir("L3")
    # seq_csv_path = "../input/input_seq_L3.csv"
    # make_xy(dirs_random_l3, "random_l3_2", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l3_2.npy"))
    # print(load_npy("../data/npy/y_random_l3_2.npy"))


    # dirs_random_l1 =  lr.load_random_dir("oxdna_random_3", "L1")
    # seq_csv_path = "../input/input_seq_L1.csv"
    # make_xy(dirs_random_l1, "random_l1_3", seq_csv_path)
    # print(len(load_npy("../data/npy/x_random_l1_3.npy")))
    # print(len(load_npy("../data/npy/y_random_l1_3.npy")))

    # dirs_random_l1 =  lr.load_random_dir("oxdna_random_6", "L1")
    # seq_csv_path = "../input/input_seq_L1.csv"
    # make_xy(dirs_random_l1, "random_l1_6", seq_csv_path)
    # print(len(load_npy("../data/npy/x_random_l1_6.npy")))
    # print(len(load_npy("../data/npy/y_random_l1_6.npy")))

    # dirs_random_l2 =  lr.load_random_dir("oxdna_random_6", "L2")
    # seq_csv_path = "../input/input_seq_L2.csv"
    # make_xy(dirs_random_l2, "random_l2_6", seq_csv_path)
    # print(len(load_npy("../data/npy/x_random_l2_6.npy")))
    # print(len(load_npy("../data/npy/y_random_l2_6.npy")))

    # dirs_random_l3 =  lr.load_random_dir("oxdna_random_6", "L3")
    # seq_csv_path = "../input/input_seq_L3.csv"
    # make_xy(dirs_random_l3, "random_l3_6", seq_csv_path)
    # print(len(load_npy("../data/npy/x_random_l3_6.npy")))
    # print(len(load_npy("../data/npy/y_random_l3_6.npy")))

    # dirs_random_l1 =  lr.load_random_dir("oxdna_random_6", "L1")
    # seq_csv_path = "../input/input_seq_L1.csv"
    # make_xy_ratioB(dirs_random_l1, "random_l1_6", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l1_6_rb.npy"))
    # print(load_npy("../data/npy/y_random_l1_6_rb.npy"))

    # dirs_random_l1 =  lr.load_random_dir("oxdna_random_6", "L1")
    # seq_csv_path = "../input/input_seq_L1.csv"
    # make_xy_ratioB(dirs_random_l1, "random_l1_6", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l1_6_dv.npy"))
    # print(load_npy("../data/npy/y_random_l1_6_dv.npy"))

    dirs_random_l1 =  lr.load_random_dir("oxdna_random_6", "L2")
    seq_csv_path = "../input/input_seq_L2.csv"
    make_xy_ratioB(dirs_random_l1, "random_l2_6", seq_csv_path)
    print(load_npy("../data/npy/x_random_l2_6_dv.npy"))
    print(load_npy("../data/npy/y_random_l2_6_dv.npy"))

    # dirs_random_l1 =  lr.load_random_dir("oxdna_random_6", "L3")
    # seq_csv_path = "../input/input_seq_L3.csv"
    # make_xy_ratioB(dirs_random_l1, "random_l3_6", seq_csv_path)
    # print(load_npy("../data/npy/x_random_l3_6_dv.npy"))
    # print(load_npy("../data/npy/y_random_l3_6_dv.npy"))


def main():
    dirs_l2 =  lr.load_l2_directory()
    seq_csv_path = "../input/input_seq_L3.csv"
    dirs_l3 =  lr.load_l3_directory()
    # make_xy(dirs_l2, "l2_1", seq_csv_path)
    # print(load_npy("../data/npy/y_l2_1.npy"))
    seq_csv_path = "../input/input_seq_L3.csv"
    make_xy(dirs_l3, "l3_1", seq_csv_path)
    print(load_npy("../data/npy/y_l3_1.npy"))

if __name__ == '__main__':
    # main()
    test()