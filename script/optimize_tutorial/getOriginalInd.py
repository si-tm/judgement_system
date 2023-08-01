import pickle
import numpy as np
import csv
import os
import glob
import sys
# from optimize_l1_original_ind import L1Individual
# from optimize_l2_original_ind import L2Individual
from optimize_l3_original_ind import L3Individual

# final.pからスコアの高いindを取り出す

def make_req(type_of_l, filename, lst, target):
    temp_f = open("../../input/template_seq/requirement_" + type_of_l + ".txt", "r")
    os.mkdir("r" + target + "/r" + filename ) # ここ変える
    new_f = open("r" + target + "/r" + filename + "/req_r" + filename + ".txt", "w") # ここ変える
    tmp_theme = ""
    for temp_l in temp_f:
        if temp_l[0] == '#':
            tmp_theme = temp_l[:-1]
            if temp_l[:-1] == "# structure":
                new_f.write("# structure\n")
                new_f.write("number_of_types = " + str(len(lst)) + "\n")
                for l in lst:
                    new_f.write(l + "\n")

        if "# structure" == tmp_theme:
            pass
        else:
            new_f.write(temp_l[:-1] + "\n")
            

def req(comp):
    new_comp = []
    for i, c in enumerate(comp):
        new_comp.append("s" + str(i) + " = " + c + " @initial 1.0 M")
    return new_comp

def Ind2complexes(lst, type_of_l):
    f = open("../../input/structure_seq/input_seq_" + type_of_l + ".csv", "r")
    r = csv.reader(f)

    comp = []

    seq_lst = []
    for l in r:
        for e in l:
            seq_lst.append(e)
    
    for i, e in enumerate(lst):
        if e == 1:
            # print(seq_lst[i])
            comp.append(seq_lst[i])
    
    return comp


# https://gitlab.com/leo.cazenille/qdpy/-/blob/master/qdpy/containers.py
def readFinal(path, type_of_l, target):
    with open(path, "rb") as f:
        data = pickle.load(f)
    # ``data`` is now a dictionary containing all results, including the final container, all solutions, the algorithm parameters, etc.
    grid = data['container']
    for ind in grid:
        if 1:
        # if ind.features[1] >= 2 and ind.features[1] <= 6: # ここ変える
            lst = ind.indexes
            print(np.sum(lst))
            # comp = Ind2complexes(lst, type_of_l)
            # comp = req(comp)
            # print(len(comp))
            # make_req(type_of_l="L3", filename=ind.name, lst=comp, target=target) # ここ変える

            # print(ind.fitness)
            # print(ind.features)

def readFinals(path, type_of_l):
    # script/optimize_tutorial/results/optimizationresults_20230724060622/final
    finals = glob.glob(path + "final/*")
    for final in finals:
        readFinal(final, type_of_l)

def test():
    path="results/optimizationresults_20230731054418/final.p" # ここ変える
    readFinal(path, "L2") # ここ変える

def test_l3():
    # path="results/optimizationresults_20230724060622/final.p" # ここ変える
    # readFinal(path, "L3") # ここ変える
    path="results/optimizationresults_20230725032045/" # ここ変える
    readFinals(path, "L3") # ここ変える

def main(target, type_of_l):
    path="results/optimizationresults_" + target + "/final.p" # ここ変える
    readFinal(path, type_of_l, target) # ここ変える

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("usage : python3 getOriginalInd.py <target name> <type of l>")
    target = sys.argv[1]
    type_of_l = sys.argv[2]
    main(target, type_of_l)
    # test()
    # test_l3()
