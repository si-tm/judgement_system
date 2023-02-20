# import numpy as np
# import matplotlib.pylab as plt
import glob
import sys

# judge if this directory has an error.


# extracting files from directory
def file_dic(dir_path):

    if "random" in dir_path:
        return random_file_dic(dir_path)

    if dir_path[-1] == "/":
        files = glob.glob(dir_path + "*")
    else:
        files = glob.glob(dir_path + "/*")

    # target seq file
    target_c = ""
    for d in dir_path.split("/"):
        if "seq" in d:
            target_c = d[-1]

    fd = {}

    for f in files:
        key = f.split("/")[-1][:f.split("/")[-1].find("_seq" + target_c + "-")]
        fd[key] = f
        if f[-4:] == ".top":
            fd["topology"] = f
        if f.split("/")[-1][:4] == "seq" + target_c:
            fd["seq"] = f
    
    # for f in fd:
    #     print(f)

    return fd

def random_file_dic(dir_path):
    # print("random")
    # input/results/oxdna_random_1/L1/d-0-6-7-4/L1_d-0-6-7-4_0/L1_d-0-6-7-4_0/
    if dir_path[-1] == "/":
        files = glob.glob(dir_path + "*")
        dir_path = dir_path[:-1]
    else:
        files = glob.glob(dir_path + "/*")
    
    # target file
    target = dir_path.split("/")[-1]
    print(dir_path)
    print(target)

    fd = {}

    for f in files:
        key = f.split("/")[-1]
        key = key[:key.find("_" + target)]
        fd[key] = f
        if f[-4:] == ".top":
            fd["topology"] = f
        if key == "hb":
            fd["hb_energy"] = f
        if "seq_req" in key:
            fd["seq"] = f
        if "req_L" in key:
            fd["req"] = f
        # if "seq" in key:
        #     fd["seq"] = f

        # if f.split("/")[-1][:4] == "seq" + target_c:
        #     fd["seq"] = f
    
    # for f in fd:
    #     print(f, fd[f])
    
    return fd

def get_conf(dir_path):
    d = file_dic(dir_path)
    # print(d["last_conf"])
    return d["last_conf"]

def get_input(dir_path):
    d = file_dic(dir_path)
    # print(d["input"])
    return d["input"]

def get_new_input(dir_path):
    d = file_dic(dir_path)
    # print(d["new_input"])
    return d["new_input"]

def get_top(dir_path):
    d = file_dic(dir_path)
    return d["topology"]

def get_seq(dir_path):
    d = file_dic(dir_path)
    # print(d["seq"])
    return d["seq"]

def get_bonds(dir_path):
    if dir_path[-1] == "/":
        return dir_path + "bonds"
    else:
        return dir_path + "/bonds"

def get_req(dir_path):
    d = file_dic(dir_path)
    # print(d)
    return d["req"]

def test():
    dir_path = "../input/results/oxdna_random_1/L1/d-0-6-7-4/L1_d-0-6-7-4_0/L1_d-0-6-7-4_0/"
    dir_path = "../input/results/oxdna_random_6_diffseq_2/L1/d-0-3-4-6-8-14/L1_d-0-3-4-6-8-14_2023-01-31-044547/L1_d-0-3-4-6-8-14_2023-01-31-044547/"
    # last_conf = get_conf(dir_path)
    # print(last_conf)
    dic = file_dic(dir_path)
    # for k in dic:
    #     # if "../input/results/oxdna_random_1/L1/d-0-6-7-4/L1_d-0-6-7-4_0/L1_d-0-6-7-4_0/trajectory_L1_d-0-6-7-4_0.dat" == dic[k]:
    #     print(k)
    #     print(dic[k])
    #     print()
    
    print(len(dic))
    # last_conf = get_conf("../../input/results/oxdna_ked/seqA/A3/test_a3_200000_1")
    # input = get_input("../../input/results/oxdna_ked/seqA/A3/test_a3_200000_1")
    # top = get_top("../../input/results/oxdna_ked/seqA/A3/test_a3_200000_1")


def main():
    # print(sys.argv[1])
    if len(sys.argv) != 3:
        print("usage : python get_target_file.py [target directory] [option1 which_file_selected]")
        print("option1 : last_conf, new_input, input, top")

    if sys.argv[2] == "last_conf":
        last_conf = get_conf(sys.argv[1])
        print(last_conf)
    if sys.argv[2] == "new_input":
        new_input = get_new_input(sys.argv[1])
        print(new_input)
    if sys.argv[2] == "input":
        input = get_input(sys.argv[1])
        print(input)
    if sys.argv[2] == "top":
        top = get_top(sys.argv[1])
        print(top)

if __name__ == '__main__':
#   main()
  test()
