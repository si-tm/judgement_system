#!/usr/bin/bash

cd script/optimize_tutorial/
# cd judgement_system/script/optimize_tutorial/
python3 optimize_l1.py
chmod 777 collect_optimizationresult.sh
results=`./collect_optimizationresult.sh`

echo $results

