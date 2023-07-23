#!/usr/bin/bash

# cd judgement_system/script/optimize_tutorial/
cd script/optimize_tutorial/
python3 optimize_l3.py
chmod 777 collect_optimizationresult.sh
results=`./collect_optimizationresult.sh`

echo $results

