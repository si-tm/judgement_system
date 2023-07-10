#!/bin/bash

lst=(
1686027494794-20
1686027495119-34
1686027495236-39
1686027495452-49
1686027495996-74
1686027496296-88
1686027498133-172
1686027501545-322
1686027502133-349
1686027502305-357
1686027502889-384
1686027503376-406
1686027503530-413
1686027504046-437
1686027504372-452
1686027505695-512
1686027506308-540
1686027507745-606
1686027507766-607
1686027508463-639
1686027510914-752
1686027512448-823
1686027514637-923
1686027514723-927
1686027514919-936
1686027515552-965
1686027515913-977
)

files=(
energy_r20230613134109.dat
energy_seq_dep_r20230613134109.dat
energy_trap_r20230613134109.dat
forces_r20230613134109.dat
generated_r20230613134109.dat
generated_r20230613134109.top
hb_energy_seq_dep_r20230613134109.dat
hb_energy_trap_r20230613134109.dat
hb_r20230613134109.dat
input_r20230613134109
input_seq_dep_r20230613134109
input_trap_r20230613134109
last_conf_r20230613134109.dat.dat
last_conf_seq_dep_r20230613134109.dat
last_conf_trap_r20230613134109.dat
log_r20230613134109.dat
log_seq_dep_r20230613134109.dat
log_trap_r20230613134109.dat
trajectory_r20230613134109.dat
trajectory_seq_dep_r20230613134109.dat
trajectory_trap_r20230613134109.dat
)

# traj=(
#     input/results/fromQD/r20230613134109/r1686027494794-20/reqs/fromQD/L1/r20230613134109/r20230613134109/trajectory_r20230613134109.dat
#     input/results/fromQD/r20230613134109/r1686027494794-20/reqs/fromQD/L1/r20230613134109/r20230613134109/trajectory_seq_dep_r20230613134109.dat
#     input/results/fromQD/r20230613134109/r1686027494794-20/reqs/fromQD/L1/r20230613134109/r20230613134109/trajectory_trap_r20230613134109.dat
# )


# for l in ${lst[@]}
# do
#     traj="r${l}/reqs/fromQD/L1/r20230613134109/r20230613134109/trajectory_r20230613134109.dat"
#     mved="r${l}/trajectory_r20230613134109.dat"

#     mv $traj $mved
    
#     traj="r${l}/reqs/fromQD/L1/r20230613134109/r20230613134109/trajectory_seq_dep_r20230613134109.dat"
#     mved="r${l}/trajectory_seq_dep_r20230613134109.dat"

#     mv $traj $mved
    
#     traj="r${l}/reqs/fromQD/L1/r20230613134109/r20230613134109/trajectory_trap_r20230613134109.dat"
#     mved="r${l}/trajectory_trap_r20230613134109.dat"

#     mv $traj $mved
    
# done

# for l in ${lst[@]}
# do
#     for f in ${files[@]}
#     do
#         before_filename=r${l}/${f}
#         echo $before_filename
#         after_filename=`echo ${before_filename//r20230613134109/r${l}}`
#         echo $after_filename
#         mv $before_filename $after_filename
#     done
# done

# for l in ${lst[@]}
# do
#     # input/results/fromQD/r20230613134109/seq_r1686027508463-639
#     cp seq_r${l} r${l}/
# done

# cp script/optimize_tutorial/r20230613134109/req_1686027494794-20.txt r1686027494794-20/req_1686027494794-20.txt
for l in ${lst[@]}
do
    cp ../../../../script/optimize_tutorial/r20230613134109/req_${l}.txt r${l}/req_r${l}.txt
done
