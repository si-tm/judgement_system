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

# def load_x():

def main():
    # get_target_x("../input/results/oxdna_ked/seqA/A4/test_a4_200000_1")
    # print(get_x()["../input/results/oxdna_ked/seqA/A4/test_a4_200000_1"])
    # print(get_x1())
    # get_x1()
    # get_x2()
    # load_x1()
    load_x2()

if __name__ == '__main__':
    main()
