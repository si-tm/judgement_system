import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
from common import get_target_file as gtf
from common import check_seq as cs
from classify_seq import make_input_seq as mis
from classify_seq import get_structure_from_pil as gsfp
import load_results as lr
import csv
import pickle

def get_target_x(target_dir):
    seq_file_name = gtf.get_seq(target_dir)
    seq = seq_file_name.split("/")[-1][3]
    # print(seq)

    seq_lst = gsfp.seq2structure(seq)
    # print(seq_lst)
    return seq_lst

def get_x1():
    dirs = lr.load_directory()
    x_dic = {}

    x_dic1 = mis.seq_dic1()
    x_dic1_sub = x_dic1

    for d in dirs:
        structurs = get_target_x(d)
        for str in structurs:
            x_dic1_sub[cs.lst2str(str)] = 1
        x_dic[d] = x_dic1_sub
        # print(x_dic1_sub)
        x_dic1_sub = mis.seq_dic1()

    with open("../data/dic/x_1_1.pkl", "wb") as tf:
        pickle.dump(x_dic,tf)
    
    return x_dic

def get_x2():
    dirs = lr.load_directory()
    x_dic = {}

    x_dic2 = mis.seq_dic2()
    x_dic2_sub = x_dic2


    for d in dirs:
        structurs = get_target_x(d)
        for str in structurs:
            print(cs.lst2str(str))
            x_dic2_sub[cs.lst2str(str)] = 1
        x_dic[d] = x_dic2_sub
        print(len(x_dic2_sub))
        x_dic2_sub = mis.seq_dic2()
    
    with open("../data/dic/x_1_2.pkl", "wb") as tf:
        pickle.dump(x_dic,tf)
    
    return x_dic

def get_x(csv_path, result_path):
    dirs = lr.load_l1_directory()
    x_dic = {}

    x_dic_l1 = mis.seq_dic(csv_path)
    x_dic_l1_sub = x_dic_l1


    for d in dirs:
        structurs = get_target_x(d)
        for str in structurs:
            print(cs.lst2str(str))
            x_dic_l1_sub[cs.lst2str(str)] = 1
        x_dic[d] = x_dic_l1_sub
        x_dic_l1_sub =mis.seq_dic(csv_path)
    
    with open(result_path, "wb") as tf:
        pickle.dump(x_dic,tf)
    
    return x_dic

def load_x1():
    with open("../data/dic/x_1_1.pkl", "rb") as tf:
        new_dict = pickle.load(tf)
        return new_dict

    print(new_dict)

def load_x2():
    with open("../data/dic/x_1_2.pkl", "rb") as tf:
        new_dict = pickle.load(tf)
        return new_dict

    print(new_dict)

def load_x(path):
    with open(path, "rb") as tf:
        new_dict = pickle.load(tf)
        print(len(new_dict))
        return new_dict


# def load_x():

def main():
    # get_target_x("../input/results/oxdna_ked/seqA/A4/test_a4_200000_1")
    # print(get_x()["../input/results/oxdna_ked/seqA/A4/test_a4_200000_1"])
    # print(get_x1())
    # get_x1()
    # get_x2()
    l1_csv_path = "../input/input_seq_L1.csv"
    l1_result_path = "../data/dic/x_l1_1.pkl"
    l2_csv_path = "../input/input_seq_L2.csv"
    l2_result_path = "../data/dic/x_l2_1.pkl"
    l3_csv_path = "../input/input_seq_L3.csv"
    l3_result_path = "../data/dic/x_l3_1.pkl"
    # get_x(l1_csv_path, l1_result_path)
    # get_x(l2_csv_path, l2_result_path)
    # get_x(l3_csv_path, l3_result_path)

    # load_x1()
    # load_x2()
    load_x(l1_result_path)
    load_x(l2_result_path)
    load_x(l3_result_path)

if __name__ == '__main__':
    main()
