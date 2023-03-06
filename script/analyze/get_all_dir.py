import sys
sys.path.append("analyze")
sys.path.append("measuring_volume")
sys.path.append("common")
import glob
import count_strands as cs
import convexhull_volume as cv
import pickle
import run_output_bonds as rob
import check_dir as cd


def apply_for_all(path, type_of_l):
    # input/results/oxdna_random_6/L1/d-0-1/L1_d-0-1_2023-01-27-083608/L1_d-0-1_2023-01-27-083608/bonds
    dirs = glob.glob(path + "/" + type_of_l + "/*/*/*/")
    dic = {}
    for d in dirs:
        if cd.included_full_files(d) == False:
            continue
        rob.make_bonds(d)
        before, after = cs.count_strands(d)
        meanv, devv = cv.convexhull_volume_all_strands_meandev(d)
        dic[d] = {}
        dic[d]["ratio_of_volume"] = after/before
        dic[d]["mean_volume"] = meanv
        dic[d]["deviation_of_volume"] = devv
        
        # print(before, after)
        print(dic[d])
    
    result_path =  "../data/dic/" + type_of_l + "_data_2.pkl"
    
    with open(result_path, "wb") as tf:
        pickle.dump(dic,tf)

    return dic


def test():
    path="../input/results/oxdna_random_6"
    type_of_l="L3"
    apply_for_all(path, type_of_l)

# def main():

if __name__ == '__main__':
    # main()
    test()