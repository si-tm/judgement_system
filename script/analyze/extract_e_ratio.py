import sys
sys.path.append("common")
import get_target_file as gtf

def seq2domain(target_dir):
    seq_path = gtf.get_seq(target_dir)
    seq_set = {}
    seq_f = open(seq_path, "r")
    for s in seq_f:
        seq_set.add(s.replace("\n", ""))
    
    print(seq_set)

def pil2req(pil_path, target_dir):
    req_path = target_dir + "/req"
    req = open(req_path, "w")
    all_e_path = target_dir + "/all_e"
    all_e = open(all_e_path , "w")
    pil_f = open(pil_path, "r")

    for l in pil_f:
        if l[0] != 'e':
            req.write(l)
        else:
            all_e.write(l)
    
    seq2domain(target_dir)
    
    print(req_path)
    print(all_e_path)
    
    req.close()
    all_e.close()
    pil_f.close()


def extract_e_ratio():
    pil_path = "../input/pils/" + "L1-entropyRT-1-22-39.pil"
    target_dir = "../input/results/oxdna_ked_2/seqA/" + "A3/test_a3_200000_2"

    # pil to req
    pil2req(pil_path, target_dir)

    # what kind of strand
    # which strands are bonded
    # which e* is corresponded?

def test():
    extract_e_ratio()

if __name__ == '__main__':
    test()
    pass