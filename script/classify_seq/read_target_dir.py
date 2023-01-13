import numpy as np
import matplotlib.pylab as plt
import make_input_seq as mis
import sys
sys.path.append('../')
from common import get_target_file as gtf

def make_seq_dic(target_dir):
    seq_f = open(gtf.get_seq(target_dir))
    seq_dic = mis.seq_dic()

    input_seq_set = {}
    for s in seq_f:
        seq = s.split("\n")[0]
        input_seq_set.add(seq)
        

def main():
    target_dir = "../../input/results/oxdna_ked/seqL/L1/test_l1_200000_1"
    make_seq_dic(target_dir)

if __name__ == '__main__':
    main()
