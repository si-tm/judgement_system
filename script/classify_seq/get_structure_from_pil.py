import numpy as np
import matplotlib.pylab as plt

class pil:
    def __init__(self):
        self.class_Alph = ["A", "B", "C", "D", "E", "F", \
            	"G", "H", "I", "J",	"K", "L"]

        self.pil_files = [ \
        "L1-entropyRT-1-22-39.pil",
 	    "L1-entropyRT-1-24-20.pil",
     	"L1-meanStruct-0-23-5.pil",
        "L1-meanStruct-2-13-3.pil",
        "L2-entropyRT-0-21-5.pil",
        "L2-entropyRT-1-36-11.pil",
        "L2-meanStruct-1-36-12.pil",
        "L2-meanStruct-2-31-18.pil",	
        "L3-entropyRT-0-25-2.pil",	
        "L3-entropyRT-1-28-33.pil",	
        "L3-entropyRT-1-29-46.pil",
        "L3-meanStruct-1-30-35.pil"
        ]

        self.seq2pil_dic = {}
        self.pil2seq_dic = {}

        for i in range(len(self.class_Alph)):
            self.pil2seq_dic[self.pil_files[i]] = self.class_Alph[i]
            self.seq2pil_dic[self.class_Alph[i]] = self.pil_files[i]

    def get_structure_pil(self, pil_file_name):
        pf = open("../../input/pils/" + pil_file_name)
        structures = []
        for l in pf:
            if l[0] == "s":
                # "="と"@initial"の間を出力
                lst = l.split(" ")
                structure = lst[lst.index("=") + 1: lst.index("@initial")]
                structures.append(structure)
        return structures

    def get_structure_seq(self, seq):
        pil_file_name = self.seq2pil_dic[seq]
        pf = open("../input/pils/" + pil_file_name)
        structures = []
        for l in pf:
            if l[0] == "s":
                # "="と"@initial"の間を出力
                lst = l.split(" ")
                structure = lst[lst.index("=") + 1: lst.index("@initial")]
                structures.append(structure)
        return structures


    def seq2pil(self, seq_name):
        return self.seq2pil_dic[seq_name]

    def seq2pil(self, pil_name):
        return self.pil2seq_dic[pil_name]

def seq2structure(target_seq):
    pil_data = pil()
    return pil_data.get_structure_seq(target_seq)


def test():
    pil_data = pil()
    # target_seq = "A"
    # target_pil = pil_data.seq2pil(target_seq)
    for target_pil in pil_data.pil_files:
        target_seq = pil_data.pil2seq_dic[target_pil]
        print(target_pil, target_seq)
        print(pil_data.get_structure_seq(target_seq))
        # print(pil_data.get_structure_pil(target_pil))


# def main():

if __name__ == '__main__':
#   main()
    test()