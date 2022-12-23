import numpy as np
import matplotlib.pylab as plt
import csv
# import pickle # pickleshare  0.7.5
import glob

# def load_file():
#     filename = '../data/sequences3.csv'
#     f2 = open("../data/sequence4.csv", "w")
#     with open(filename, encoding='utf8', newline='') as f:
#         csvreader = csv.reader(f)
#         for row in csvreader:
#             print(row)
#             for seq in row:
#                 if seq != "":
#                     print(seq)
#                     f2.write(seq + "\n")
    
#     f2.close()

# def make_dic():
#     seq = {}
#     f = open("../data/sequences4.csv", "r")
#     csvreader = csv.reader(f)
#     for s in csvreader:
#         print(s)
#         seq[s[0]] = 0

# def load_seq():
#     f = open("../input/sequences/seqA/seqA-GA100000-0.80_final_20200904131038.dat")
#     for seq in f:
#         print(seq[:-1]) #eliminate /n
#     f.close()

class input_seq_data:
    def __init__(self, seq_file_name):
        self.seq_file_name = seq_file_name
        self.seq = {}
        f = open("../data/sequences4.csv", "r")
        csvreader = csv.reader(f)
        for s in csvreader:
            # print(s)
            self.seq[s[0]] = 0
    
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

def main():
    # seq_a_1
    seq_a_1 = input_seq_data("../input/sequences/seqA/seqA-GA100000-0.80_final_20200904131038.dat")
    seq_a_1.input_seq_num()
    seq_a_1.verification()
    # get all seq file
    # seq_files = glob.glob("../input/sequences/*/*")

if __name__ == '__main__':
  main()
