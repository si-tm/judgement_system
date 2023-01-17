# make execusion directory and results directory
if [ -e {../execusion/} ]; then
    echo "execution directory exit"
else
    mkdir ../execusion/
fi
if [ -e {../results/} ]; then
    echo "results directory exit"
else
    mkdir ../results/
fi
# make top, conf file from sequence
python ../utils/generate_sa.py 
# move input and sequence files to execusion directory
# modify input files
# run files
# move results to result directory
