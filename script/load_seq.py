import numpy as np
import matplotlib.pylab as plt
import csv
# import pickle # pickleshare  0.7.5
import glob
from common import get_target_file as gf

class input_seq_data:
    def __init__(self, seq_file_name):
        self.seq_file_name = seq_file_name
        self.seq = {}
        f = open("../data/sequences4.csv", "r")
        csvreader = csv.reader(f)
        for s in csvreader:
            # print(s)
            self.seq[s[0]] = 0

    def get_seq(self):
        return self.seq
    
    def input_seq_num(self):
        f = open(self.seq_file_name)
        for s2 in f:
            # print(s2[:-1]) #eliminate \n
            for s1 in self.seq:
                if self.is_include(s2[:-1], s1):
                    self.seq[s1] += 1
        f.close()

    def verification(self):
        for s in self.seq:
            if self.seq[s] != 0:
                print(s, self.seq[s])
    
    def is_include(self, seq1, seq2):
        # if seq1 includes seq2, it return True
        if seq2 in seq1:
            return True
        return False

def get_seq_data(target_dir):
    seq_file = gf.get_seq(target_dir)
    seq_dic_data = input_seq_data(seq_file)
    return seq_dic_data.seq


def main():
    dic = get_seq_data("../input/results/oxdna_ked/seqA/A2/test_a2_200000_2")
    for d in dic:
        print(d, dic[d])

if __name__ == '__main__':
  main()
