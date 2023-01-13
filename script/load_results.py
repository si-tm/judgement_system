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

def main():
    dirs = load_directory()

if __name__ == '__main__':
    main()