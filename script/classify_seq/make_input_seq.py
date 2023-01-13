import csv
import sys
sys.path.append('../')
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

def main():
    seq_lst1()
    seq_lst2()

if __name__ == '__main__':
  main()