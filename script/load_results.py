import numpy as np
import matplotlib.pylab as plt
import sys
import glob
import sys
sys.path.append('common')
sys.path.append('measuring_volume')
import check_dir as cd
import measuring_volume.convexhull_volume2 as cv
import get_target_file as gtf


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

def load_random_dir(target_dir, type_of_l):
    folder = glob.glob("../input/results/oxdna_random_2/" + type_of_l + "/*/*/*/")
    folder = glob.glob("../input/results/" + target_dir + "/" + type_of_l + "/*/*/*/")

    right_dir = []
    for f in folder:
        
        if cd.included_full_files(f) == True:
            right_dir.append(f)

    return right_dir

def load_diffseq_dir(path="../input/results/oxdna_random_6_diffseq/L1"):
    # input/results/oxdna_random_6_diffseq/L1/d-0-1/L1_d-0-1_2023-01-30-152202/L1_d-0-1_2023-01-30-152202/energy_L1_d-0-1_2023-01-30-152202.dat
    folder = glob.glob(path + "/*/*/*/")

    right_dir = []
    for f in folder:
        fd = gtf.file_dic(f).keys()
        if "last_conf" in fd and "input" in fd:
            if cd.included_full_files(f) == True:
                right_dir.append(f)

    return right_dir


def main():
    # dirs = load_directory()
    # dirs_l3 = load_l3_directory()
    # print(dirs_l3)
    # dirs_random_l1 = load_random_dir("L1")
    print(len(load_diffseq_dir(path="../input/results/oxdna_random_6_diffseq_2/L1")))

if __name__ == '__main__':
    main()
    