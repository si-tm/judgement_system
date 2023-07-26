#!/usr/bin/bash
echo `ls`
# cd judgement_system/script/optimize_tutorial/
cd script/optimize_tutorial/
python3 optimize_l2.py
chmod 777 collect_optimizationresult.sh
results=`./collect_optimizationresult.sh`

echo $results

