import numpy as np
import matplotlib.pylab as plt
import sys
import glob
from common import check_dir as cd
from measuring_volume import convexhull_volume as cv

def load_directory():
    folder = glob.glob("../input/results/*/*/*/*/")
    right_dir = []
    for f in folder:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)

    return right_dir

def load_l1_directory():
    folderA = glob.glob("../input/results/oxdna_ked_2/*A/*/*/")
    folderB = glob.glob("../input/results/oxdna_ked_2/*B/*/*/")
    folderC = glob.glob("../input/results/oxdna_ked_2/*C/*/*/")
    folderD = glob.glob("../input/results/oxdna_ked_2/*D/*/*/")
    right_dir = []
    for f in folderA:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderB:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderC:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderD:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)

    return right_dir

def load_l2_directory():
    folderE = glob.glob("../input/results/oxdna_ked_2/*E/*/*/")
    folderF = glob.glob("../input/results/oxdna_ked_2/*F/*/*/")
    folderG = glob.glob("../input/results/oxdna_ked_2/*G/*/*/")
    folderH = glob.glob("../input/results/oxdna_ked_2/*H/*/*/")
    right_dir = []
    for f in folderE:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderF:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderG:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderH:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)

    return right_dir

def load_l3_directory():
    folderI = glob.glob("../input/results/oxdna_ked_2/*I/*/*/")
    folderJ = glob.glob("../input/results/oxdna_ked_2/*J/*/*/")
    folderK = glob.glob("../input/results/oxdna_ked_2/*K/*/*/")
    folderL = glob.glob("../input/results/oxdna_ked_2/*L/*/*/")
    right_dir = []
    for f in folderI:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderJ:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderK:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)
    for f in folderL:
        if cd.included_full_files(f) == True:
            # print(cd.included_full_files(f))
            right_dir.append(f)

    return right_dir

def main():
    # dirs = load_directory()
    dirs_l3 = load_l3_directory()
    print(dirs_l3)

if __name__ == '__main__':
    main()