import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
from common import get_target_file as gtf
from common import check_seq as cs
from common import check_dir as cd
from classify_seq import make_input_seq as mis
from classify_seq import get_structure_from_pil as gsfp
from classify_seq import get_structure_from_req as gsfq
import load_results as lr
import csv
import pickle

def get_target_x(target_dir):
    seq_file_name = gtf.get_seq(target_dir)
    seq = seq_file_name.split("/")[-1][3]
    # print(seq)

    seq_lst = []

    if cd.is_random(target_dir):
        seq_lst = gsfq.seq2structure(target_dir)
    else:
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

def get_x(csv_path, result_path, dirs):
    
    x_dic = {}

    x_dic_l1 = mis.seq_dic(csv_path)
    x_dic_l1_sub = x_dic_l1

    # print(x_dic_l1)


    for d in dirs:
        print(d)
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
        print(new_dict)
        return new_dict


# def load_x():

def test():
    target_dir = "../input/results/oxdna_ked_2/seqA/A1/test_a1_10_4"
    get_target_x(target_dir)

def make_l():
    l1_csv_path = "../input/input_seq_L3.csv"
    l1_result_path = "../data/dic/x_random_l3_5.pkl"
    dirs_l1 = lr.load_random_dir("oxdna_random_5", "L3")
    # print(dirs_l1)
    get_x(l1_csv_path, l1_result_path, dirs_l1)
    load_x(l1_result_path)

def main():
    # get_target_x("../input/results/oxdna_ked/seqA/A4/test_a4_200000_1")
    # print(get_x()["../input/results/oxdna_ked/seqA/A4/test_a4_200000_1"])
    # print(get_x1())
    # get_x1()
    # get_x2()
    l1_csv_path = "../input/input_seq_L1.csv"
    # l1_result_path = "../data/dic/x_l1_1.pkl"
    l1_result_path = "../data/dic/x_l1_random_1.pkl"
    l2_csv_path = "../input/input_seq_L2.csv"
    # l2_result_path = "../data/dic/x_l2_1.pkl"
    l2_result_path = "../data/dic/x_l2_random_1.pkl"
    l3_csv_path = "../input/input_seq_L3.csv"
    # l3_result_path = "../data/dic/x_l3_1.pkl"
    l3_result_path = "../data/dic/x_l3_random_1.pkl"

    # dirs_l1 = lr.load_l1_directory()
    # dirs_l2 = lr.load_l2_directory()
    # dirs_l3 = lr.load_l3_directory()

    dirs_l1 = lr.load_random_dir("L1")
    # print(dirs_l1)
    dirs_l2 = lr.load_random_dir("L2")
    dirs_l3 = lr.load_random_dir("L3")

    get_x(l1_csv_path, l1_result_path, dirs_l1)
    get_x(l2_csv_path, l2_result_path, dirs_l2)
    get_x(l3_csv_path, l3_result_path, dirs_l3)

    # load_x1()
    # load_x2()
    load_x(l1_result_path)
    load_x(l2_result_path)
    load_x(l3_result_path)

if __name__ == '__main__':
    # main()
    # test()
    make_l()
