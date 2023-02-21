import get_x as gx
import get_y as gy
import load_results as lr
import pickle
from classify_seq import make_input_seq as mis
import numpy as np
import make_xy as mxy

# x_pkl : {dir_name : [0, 1, 0 ... 0]}
# y_pkl : {dir_name : volume}

def dir2target(dir):
    # input/results/oxdna_random_6_diffseq/L1/d-0-1/L1_d-0-1_2023-01-30-152202/L1_d-0-1_2023-01-30-152202/energy_L1_d-0-1_2023-01-30-152202.dat
    return dir.split("/")[-4]

def mean_xy_dic(x_dic_path, y_dic_path, dirs):
    new_x_dic = {}
    new_y_dic = {}
    new_dirs_lst = []

    x_dic = gx.load_x(x_dic_path)
    y_dic= gy.load_y(y_dic_path)

    target2dir_dic = {}
    dir2target_dic = {}


    # { d_name : dirs }を作る
    # x → {d_name : [0, 1, ... 0]}
    # y → {d_name : mean volume}

    for dir in dirs:
        target_name = dir2target(dir)
        new_dirs_lst.append(target_name)
        if not target_name in target2dir_dic:
            target2dir_dic[target_name] = []
        target2dir_dic[target_name].append(dir)
        dir2target_dic[dir] = target_name

    # print(len(target2dir_dic))
    # print(len(dir2target_dic))
    
    # # x → {d_name : [0, 1, ... 0]}
    for d in dirs:
        new_x_dic[dir2target_dic[d]] = x_dic[d]

    # y → {d_name : mean volume}
    for new_d in new_dirs_lst:
        count = 0.0
        sum_volume = 0.0
        for d in target2dir_dic[new_d]:
            count += 1.0
            sum_volume += y_dic[d]
        new_y_dic[new_d] = count/sum_volume

    print(len(new_x_dic))
    print(len(new_x_dic))
    
    x_result_path = "../data/dic/x_random_6_diffseq_mean_l1_3.pkl"
    y_result_path = "../data/dic/y_random_6_diffseq_mean_l1_3.pkl"
    
    with open(x_result_path, "wb") as tf:
        pickle.dump(new_x_dic,tf)
    
    with open(y_result_path, "wb") as tf:
        pickle.dump(new_y_dic,tf)
    
    new_dirs_lst = list(target2dir_dic.keys())
    print(type(new_dirs_lst))
    
    return new_dirs_lst


def mean_make_xy(dirs):


    seq_csv = "../input/input_seq_L1.csv"
    dic_x = gx.load_x("../data/dic/x_random_6_diffseq_mean_l1_3.pkl")
    dic_y = gy.load_y("../data/dic/y_random_6_diffseq_mean_l1_3.pkl")
    x_lst = mis.seq_lst(seq_csv)

    x = []
    y = []

    print("dir_num", len(dirs))
    print(len(dic_x))
    print(len(dic_y))

    for d in dirs:
        # print(d)
        x.append(mxy.dic2lst(dic_x[d], x_lst))
        y.append(dic_y[d])
    
    print(len(x))
    print(len(y))
    
    np.save("../data/npy/x_random_6_diffseq_mean_l1_3.npy",x)
    np.save("../data/npy/y_random_6_diffseq_mean_l1_3.npy",y)


def main():
    x_dic_path="../data/dic/x_random_6_diffseq_l1_3.pkl"
    y_dic_path="../data/dic/y_random_6_diffseq_l1_3.pkl"
    # dirs = lr.load_diffseq_dir(path="../input/results/oxdna_random_6_diffseq_2/L1")
    dirs = lr.load_diffseq_dir(path="../input/results/oxdna_random_6_diffseq_3/L1")
    print(len(dirs))
    new_dirs = mean_xy_dic(x_dic_path, y_dic_path, dirs=dirs)
    # print("d-3-6-11-12-13" in new_dirs)
    mean_make_xy(new_dirs)


if __name__ == '__main__':
    main()
