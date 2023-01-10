import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
from common import get_target_file as gtf

def get_conf_data(target_dir):

    conf_f = open(gtf.get_conf(target_dir))
    particle_id = -3

    conf_header_elem = ["rx", "ry", "rz", "bx", "by", "bz", "nx", "ny", "nz", "vx", "vy", "vz", "Lx", "Ly", "Lz"]
    conf_header = ["center-of-mass_position_r", "base_vector_a1", "base_normal_vector_a3", "Velocity", "Angular_velocity"]

    # conf_dic[particle_id]["rz"] = val
    conf_dic = {}

    for l in conf_f:
        if particle_id >= 0:
            lst = l.split(" ")
            lst[-1] = lst[-1][:-1]
            conf_dic[particle_id] = {}

            for i in range(len(lst)):
                conf_dic[particle_id][conf_header_elem[i]] = float(lst[i])
            
        particle_id += 1
    
    return conf_dic

def main():
    target_dir = "../../input/results/oxdna_ked/seqA/A4/test_a4_200000_1"
    print(get_conf_data(target_dir))

if __name__ == '__main__':
  main()
