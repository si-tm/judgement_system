#!/bin/zsh

get_conf_input () {
    python3 ../common/check_dir.py $1
    traj=`python3 ../common/get_target_file.py $1 last_conf`
    input=`python3 ../common/get_target_file.py $1 new_input`
    top=`python3 ../common/get_target_file.py $1 top`
}

# target_dir="../results_KakenhiEvolveDNA/seqA/A4/test_a4_200000_1/"
# target_dir="../results_KakenhiEvolveDNA/seqL/L2/test_l2_200000_1"
target_dir="../../input/results/oxdna_ked/seqA/A4/test_a4_200000_1"

get_conf_input $target_dir
echo "\n"
echo "trajectory file path is "$traj"\n"
echo "input file path is "$input"\n"
echo "topology file path is "$top"\n"

python3 output_bonds.py $input $traj $top
# python3 --version
