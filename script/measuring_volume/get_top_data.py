import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
sys.path.append('common')
# from common import get_target_file as gtf
import common.get_target_file as gtf
import csv 
import glob

def make_initial_strands_data(target_dir):
    # particle2strands[particle_id] = strand_id
    particle2strand = {} 
    # strands2particle[strand_id] = {particle_ids}
    strands2particle = {} 
    top_f = open(gtf.get_top(target_dir))

    col = 0
    particle_id = 0
    initial_strand_num = -1

    for l in top_f:
        if col == 0:
            initial_strand_num = int(l.split(" ")[-1].rstrip('\n'))
        else:
            strand_id = int(l.split(" ")[0])
            if strand_id in strands2particle:
                strands2particle[strand_id].add(particle_id)
            else:
                strands2particle[strand_id] = {particle_id}
            particle2strand[particle_id] = strand_id
            particle_id += 1
        col += 1
    
    # print(initial_strand_num)
    
    top_f.close()
    return strands2particle, particle2strand

def get_bonds_data(bonds_name):
    bonds_f = open(bonds_name)
    line = 0

    header = ["id1", "id2", "FENE", "BEXC", "STCK", "NEXC", "HB", "CRSTCK", "CXSTCK", "total"]
    # ex. bonds_dic[((id1, id2), "FENE")] = value
    bonds_dic = {}

    for l in bonds_f:
        if line > 0 and l[0] != '#':
            lst = l.split(" ")
            lst[-1] = lst[-1][:-2]
            id1 = int(lst[0])
            id2 = int(lst[1])
            for index, h in enumerate(header[2:]):
                # print(id1, id2, h, float(lst[index + 2]))
                if lst[index + 2] == '-':
                    continue
                bonds_dic[((id1, id2), h)] = float(lst[index + 2])
        line += 1
    return bonds_dic

def get_connection_strands(bonds_name, strands2particle, particle2strand):
    
    bonds_dic = get_bonds_data(bonds_name)
    # HB < 0.0同士のparticle id1, id2である時、strand idが小さい方に合体させる
    
    for b in bonds_dic:
        
        if b[1] == "HB" and bonds_dic[b] < 0.0:
            particle_id1 = b[0][0]
            particle_id2 = b[0][1]
            # print("bonds : ", particle_id1, particle_id2)
            strand_id1 = particle2strand[particle_id1]
            strand_id2 = particle2strand[particle_id2]
            # strand idが異なる場合、小さい方に合わせる
            # if strand_id2 == strand_id1:
            if strand_id2 != strand_id1:
                print(strand_id2, strand_id1)
                s_id1 = min(strand_id1, strand_id2)
                s_id2 = max(strand_id1, strand_id2)
                if s_id2 in strands2particle:
                    # s_id2 → s_id1へ
                    for particle in list(strands2particle[s_id2]):
                        particle2strand[particle] = s_id1

                    strands2particle[s_id1] |= strands2particle[s_id2]
                    strands2particle.pop(s_id2)

    return strands2particle, particle2strand


def get_connection_strands2(bonds_name, strands2particle, particle2strand):
    
    bonds_dic = get_bonds_data(bonds_name)
    # HB < 0.0同士のparticle id1, id2である時、接続している
    # {1:[2], 2: [1], 10:[14, 20], ...} strand1とstrand2が接続している場合
    connected_strands = {}
    
    for b in bonds_dic:
        
        if b[1] == "HB" and bonds_dic[b] < 0.0:
            particle_id1 = b[0][0]
            particle_id2 = b[0][1]
            strand_id1 = particle2strand[particle_id1]
            strand_id2 = particle2strand[particle_id2]
            # strand idが異なる場合
            if strand_id2 != strand_id1:

                # print("Wow")
               
                if strand_id1 in connected_strands and not (strand_id2 in connected_strands[strand_id1]):
                    connected_strands[strand_id1].append(strand_id2)
                else:
                    connected_strands[strand_id1] = [strand_id2]

                if strand_id2 in connected_strands and not (strand_id1 in connected_strands[strand_id2]):
                    connected_strands[strand_id2].append(strand_id1)
                else:
                    connected_strands[strand_id2] = [strand_id1]
               
    return  connected_strands

def get_particle_strands_data(target_dir):
    strands2particle, particle2strand = make_initial_strands_data(target_dir)
    # bondsファイルを取得する
    # print(target_dir)
    # print(gtf.get_bonds(target_dir))
    strands2particle, particle2strand = get_connection_strands(gtf.get_bonds(target_dir), strands2particle, particle2strand)
    return strands2particle, particle2strand 

def get_connected_strands_data(target_dir):
    strands2particle, particle2strand = make_initial_strands_data(target_dir)
    # bondsファイルを取得する
    connected_strands = get_connection_strands2(gtf.get_bonds(target_dir), strands2particle, particle2strand)
    return connected_strands

def main():
    # target_dir = "../results_KakenhiEvolveDNA/seqA/A4/test_a4_200000_1"
    # target_dir = "../results_KakenhiEvolveDNA/seqL/L2/test_l2_200000_1"
    # target_dir="../results_KakenhiEvolveDNA/seqL/L14/test_l14_200000_1"
    # target_dir="../input/results/oxdna_ked/seqA/A1/test_a1_200000_1"

    target_dirs = glob.glob("../input/results/oxdna_random_6_diffseq_2/L*/*/*/*/")

    for target_dir in target_dirs:
        print(target_dir)
        get_particle_strands_data(target_dir)
    
    # get_particle_strands_data(target_dir)
    # get_connected_strands_data(target_dir)
    
if __name__ == '__main__':
  main()