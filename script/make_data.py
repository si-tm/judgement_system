import sys
sys.path.append("analyze")
sys.path.append("measuring_volume")
sys.path.append("common")
import pickle
import glob
import analyze.count_strands as cs
# import convexhull_volume as cv
import measuring_volume.convexhull_volume2 as cv
import measuring_volume.run_output_bonds as rob
import check_dir as cd
import get_target_file as gtf
import datetime

def make_data(path, type_of_l, version=4):
    # input/results/oxdna_random_6/L1/d-0-1/L1_d-0-1_2023-01-27-083608/L1_d-0-1_2023-01-27-083608/bonds
    dirs = glob.glob(path + "/" + type_of_l + "/*/*/*/")
    dic = {}
    for index, d in enumerate(dirs):
        print(d)
        print(datetime.datetime.now(), " : ", index, "/", len(dirs))
        if cd.included_full_files(d) == False:
            continue
        # make bonds
        rob.make_bonds(d)
        before, after = cs.count_strands(d)
        meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
        dic[d] = {}
        dic[d]["ratio_of_volume"] = after/before
        dic[d]["mean_volume"] = meanv
        dic[d]["deviation_of_volume"] = devv
        
        # print(before, after)
        print(dic[d])
    
    result_path =  "../data/dic/" + type_of_l + "_data_" + str(version) + ".pkl"
    
    with open(result_path, "wb") as tf:
        pickle.dump(dic,tf)

    return dic


def make_data_fromQD(path, type_of_l, version=4):
    # input/results/fromQD/r20230613134109/r1686027494794-20/*
    dirs = glob.glob(path + "*")
    dic = {}
    for index, d in enumerate(dirs):
        print(datetime.datetime.now(), " : ", index + 1, "/", len(dirs))
        if cd.included_full_files(d) == False:
            continue
        print(d)
        # make bonds
        rob.make_bonds(d)
        before, after = cs.count_strands(d)
        meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
        dic[d] = {}
        dic[d]["ratio_of_volume"] = after/before
        dic[d]["mean_volume"] = meanv
        dic[d]["deviation_of_volume"] = devv
        
        # print(before, after)
        print(dic[d])
    
    result_path =  "../data/dic/" + type_of_l + "_data_fromQD_" + str(version) + ".pkl"
    
    with open(result_path, "wb") as tf:
        pickle.dump(dic,tf)

    return dic


def test():
    path="../input/results/oxdna_random_6"
    type_of_l="L1"
    make_data(path, type_of_l)

def make_fromQD():
    path="../input/results/fromQD/r20230613134109/"
    make_data_fromQD(path, type_of_l="L1", version=1)

def main():
    if len(sys.argv) != 4:
        print("usage : python make_data.py [path] [type_of_path] [version]")
    
    path=sys.argv[1]
    type_of_l=sys.argv[2]
    version=int(sys.argv[3])
    # print(path, type_of_l, version)
    make_data(path, type_of_l, version)
    

if __name__ == '__main__':
    # main()
    # test()
    make_fromQD()