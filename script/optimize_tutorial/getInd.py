import pickle
import numpy as np
import csv
import os
import glob

# final.pからスコアの高いindを取り出す

def make_req(type_of_l, filename, lst):
    temp_f = open("../../input/template_seq/requirement_" + type_of_l + ".txt", "r")
    os.mkdir("r20230730163921/r" + filename ) # ここ変える
    new_f = open("r20230730163921/r" + filename + "/req_r" + filename + ".txt", "w") # ここ変える
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

def Ind2lst(ind):
    indexes = np.array(ind[:-1]) > ind[-1]
    strands = [1 if a else 0 for a in indexes]
    return strands


# https://gitlab.com/leo.cazenille/qdpy/-/blob/master/qdpy/containers.py
def readFinal(path, type_of_l):
    with open(path, "rb") as f:
        data = pickle.load(f)
    # ``data`` is now a dictionary containing all results, including the final container, all solutions, the algorithm parameters, etc.
    grid = data['container']
    for ind in grid:
        # if ind.features[1] >= 2 and ind.features[1] <= 6: # ここ変える
        if ind.features[1] >= 2 and ind.features[1] <= 6: # ここ変える
            lst = Ind2lst(ind)
            comp = Ind2complexes(lst, type_of_l)
            comp = req(comp)
            print(comp)
            # make_req(type_of_l="L3", filename=ind.name, lst=comp) # ここ変える

            # print(ind.fitness)
            # print(ind.features)

def readFinals(path, type_of_l):
    # script/optimize_tutorial/results/optimizationresults_20230724060622/final
    finals = glob.glob(path + "final/*")
    for final in finals:
        readFinal(final, "L3")

def test():
    path="results/optimizationresults_20230730163921/final.p" # ここ変える
    readFinal(path, "L1") # ここ変える

def test_l3():
    # path="results/optimizationresults_20230724060622/final.p" # ここ変える
    # readFinal(path, "L3") # ここ変える
    path="results/optimizationresults_20230725032045/" # ここ変える
    readFinals(path, "L3") # ここ変える


if __name__ == '__main__':
    test()
    # test_l3()
    pass
