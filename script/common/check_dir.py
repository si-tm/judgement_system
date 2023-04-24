# import numpy as np
# import matplotlib.pylab as plt
import sys
sys.path.append('common')
import common.get_target_file as gtf

def is_random(target_dir):
    if "random" in target_dir:
        return True
    else:
        return False

def included_full_files(target_dir):
    fd = gtf.file_dic(target_dir)
    return "last_conf" in fd and "input" in fd and "trajectory" in fd

def check_file(file_name):
    file = open(file_name, "r")
    for l in file:
        print(l)
    file.close()

def new_input(target_dir):
    fd = gtf.file_dic(target_dir)

    new_input_file_name = "/".join(fd["input"].split("/")[:-1]) + "/new_" + fd["input"].split("/")[-1]
    print(new_input_file_name + " is created")
    new_input_file = open(new_input_file_name, "w")
    input_file = open(fd["input"], "r")

    for l in input_file:
        new_l = ""
        lst = l.split(" ")
        if lst[0] == "topology":
            new_l = "topology" + " = " + fd["topology"] + "\n"
        elif lst[0] == "conf_file":
            new_l = "conf_file" + " = " + fd["generated"] + "\n"
        elif lst[0] == "lastconf_file":
            new_l = "lastconf_file" + " = " + fd["last_conf"] + "\n"
        elif lst[0] == "trajectory_file":
            new_l = "trajectory_file" + " = " + fd["trajectory"] + "\n"
        elif lst[0] == "log_file":
            new_l = "log_file" + " = " + fd["log"] + "\n"
        elif lst[0] == "energy_file":
            new_l = "energy_file" + " = " + fd["energy"] + "\n"
        elif lst[0] == "\tname":
            new_l = "\tname" + " = " + fd["hb_energy"] + "\n"
        else:
            new_l = l
        # print(new_l, end="")
        new_input_file.write(new_l)
    
    new_input_file.close()
    input_file.close()

    return new_input_file_name
    
# def test():
    # included_full_files("../results_KakenhiEvolveDNA/seqA/A4/test_a4_200000_1/")
    # included_full_files("../results_KakenhiEvolveDNA/seqJ/J1/test_j1_200000_1/")
    # new_input("../results_KakenhiEvolveDNA/seqA/A4/test_a4_200000_1/")

def main():
    if len(sys.argv) != 2:
        print("usage : python check_dir.py [target directory]")
    else:
        target_dir = sys.argv[1]

        if included_full_files(target_dir):
            print("this directory includes full files")
        else:
            print("this directory does not includes last_conf files")
        
        new_input_file = new_input(target_dir)
        # check_file(new_input_file)

if __name__ == '__main__':
    # test()
    main()
    