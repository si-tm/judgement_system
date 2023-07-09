#!/usr/bin/bash

cd script/optimize_tutorial/
python3 optimize_l2.py
chmod 777 collect_optimizationresult.sh
results=`./collect_optimizationresult.sh`

echo $results

