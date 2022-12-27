import numpy as np
import matplotlib.pylab as plt
import csv

def count_csv(path):
    f = open(path, "r")
    f_csv= csv.reader(f)

    count = 0

    for s in f_csv:
        count += len(s)

    print(count)

    return count

def seq_count():
    return count_csv("../../data/sequences4.csv")

def main():
    seq_count()

if __name__ == '__main__':
  main()