import sys
sys.path.append("common")
import get_target_file as gtf
import csv

def seq2domain(target_dir, domain_path, file_no):
    seq_path = gtf.get_seq(target_dir)
    seq_set = set()
    seq_f = open(seq_path, "r")
    for s in seq_f:
        seq_set.add(s.replace("\n", ""))
    
    seq_lst = list(seq_set)
    print(seq_lst)

    domain_f = open(domain_path, "r")
    domain_f = csv.reader(domain_f)

    domain_dic = {}
    header = []

    for (index, domain) in enumerate(domain_f):
        if index == 0:
            header = domain
        if domain[0] == file_no:
            for (index, elem) in enumerate(domain):
                domain_dic[header[index]] = elem

    return  seq_lst, domain_dic


def pil2req(pil_path, target_dir):
    req_path = target_dir + "/req"
    req = open(req_path, "w")
    all_e_path = target_dir + "/all_e"
    all_e = open(all_e_path , "w")
    pil_f = open(pil_path, "r")
    domain_path = "../input/domain" + "/a_domain.csv"

    type_of_a2l = target_dir.split("/")[-2][0].lower()
    file_no = target_dir.split("/")[-2][-1]

    for l in pil_f:
        if l[0] != 'e':
            req.write(l)
        else:
            all_e.write(l)
    
    seq_lst, domain_dic = seq2domain(target_dir, domain_path, file_no)

    req.write("\n")
    req.write("# sequence of domains\n")
    for key in domain_dic:
        if len(key) == 1:
            req.write("domain " + key + " = " + domain_dic[key] + "\n")
    
    print(req_path)
    print(all_e_path)
    
    req.close()
    all_e.close()
    pil_f.close()

def reqseq2eratio(target_dir):
    req_path = target_dir + "/req"
    seq_path = gtf.get_seq(target_dir)

    req_f = open(req_path, "r")
    for l in req_f:
        if "domain" == l[:6]:
            domain_name = l.split(" ")[1]
            domain_sequence = l.split(" ")[-1].replace("\n", "")
            print(domain_name, domain_sequence)



def extract_e_ratio():
    pil_path = "../input/pils/" + "L1-entropyRT-1-22-39.pil"
    target_dir = "../input/results/oxdna_ked_2/seqA/" + "A3/test_a3_200000_2"

    # pil to req
    pil2req(pil_path, target_dir)

    # what kind of strand
    reqseq2eratio(target_dir)

    # which strands are bonded
    # which e* is corresponded?

def test():
    extract_e_ratio()

if __name__ == '__main__':
    test()
    pass