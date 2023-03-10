import sys
sys.path.append("common/")
sys.path.append("measuring_volume/")
import run_output_bonds as rob
import get_target_file as gtf
import convexhull_volume as cv
import get_top_data as gtd

# strand数の前後を書く．
def count_strands(target_dir):
    strands2particle, particle2strand = gtd.make_initial_strands_data(target_dir)
    # print("before : ", len(strands2particle))
    before = len(strands2particle)
    strands2particle, particle2strand = gtd.get_connection_strands(gtf.get_bonds(target_dir), strands2particle, particle2strand)
    # print("after : ", len(strands2particle))
    after = len(strands2particle)

    return (before, after)

# strand数の前後を書く．
def analyze_strands(target_dir):
    strands2particle, particle2strand = gtd.make_initial_strands_data(target_dir)
    strands2particle, particle2strand = gtd.get_connection_strands(gtf.get_bonds(target_dir), strands2particle, particle2strand)

    # which paricle are bonded
    # for strand in strands2particle:
        # print(strand)
        # print(strand, strands2particle[strand])
    # print(particle2strand)


def test():
    # target_dir = "../input/results/oxdna_random_6_diffseq_3/L1/d-0-1/L1_d-0-1_2023-01-30-152202/L1_d-0-1_2023-01-30-152202/"
    target_dir = "../input/results/oxdna_ked_2/seqD/D1/test_d1_200000_1"
    before, after = count_strands(target_dir)
    print(before, after)
    analyze_strands(target_dir)

# def main():

if __name__ == '__main__':
    # main()
    test()