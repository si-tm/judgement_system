import csv
import sys
sys.path.append('../')
sys.path.append('.')
from common import check_seq as cs

def seq_lst1():
    input_seq_f = open("../input/input_seq1.csv", "r")
    input_seq_csv = csv.reader(input_seq_f)

    seq_lst = []

    for lst in input_seq_csv:
        seq = cs.lst2str(lst)
        seq_lst.append(seq)

    return seq_lst

def seq_dic1():
    input_seq_f = open("../input/input_seq1.csv", "r")
    input_seq_csv = csv.reader(input_seq_f)

    # seq_dic[string of seq] = 0 or 1
    seq_dic = {}

    for lst in input_seq_csv:
        seq = cs.lst2str(lst)
        seq_dic[seq] = 0

    return seq_dic

def seq_lst2():
    input_seq_f = open("../input/input_seq2.csv", "r")
    input_seq_csv = csv.reader(input_seq_f)

    seq_lst = []

    for lst in input_seq_csv:
        seq = cs.lst2str(lst)
        seq_lst.append(seq)

    return seq_lst

def seq_dic2():
    input_seq_f = open("../input/input_seq2.csv", "r")
    input_seq_csv = csv.reader(input_seq_f)

    seq_dic = {}

    for lst in input_seq_csv:
        seq = cs.lst2str(lst)
        seq_dic[seq] = 0

    return seq_dic

def seq_lst(path):
    input_seq_f = open(path, "r")
    input_seq_csv = csv.reader(input_seq_f)

    seq_lst = []

    for lst in input_seq_csv:
        seq = cs.lst2str(lst)
        seq_lst.append(seq)

    return seq_lst

def seq_dic(path):
    input_seq_f = open(path, "r")
    input_seq_csv = csv.reader(input_seq_f)

    seq_dic = {}

    for lst in input_seq_csv:
        seq = cs.lst2str(lst)
        seq_dic[seq] = 0

    return seq_dic



def main():
    print(seq_lst("../input/input_seq_L1.csv"))
    print(seq_dic("../input/input_seq_L1.csv"))

if __name__ == '__main__':
  main()