#!/usr/bin/bash
# # !/bin/zsh
path="script/optimize_tutorial/"
file_LIST=(
    "activityGrid.pdf"     
    "evals_contsize.pdf"     
    "evals_fitnessmax0.pdf"    
    "final.p"   
    "iterations_nbupdated.pdf"    
    "performancesGrid.pdf"
)
DATE=`/bin/date '+%Y%m%d%H%M%S'`
dir="optimizationresults_$DATE"

# echo $file_LIST
# echo $DATE
# echo $dir

/bin/mkdir $dir

# echo $file_LIST

for f in $file_LIST
do 
    echo $dir"/"$f
    /bin/mv $f $dir"/"$f
done