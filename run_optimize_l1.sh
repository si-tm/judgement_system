#!/usr/bin/bash

/usr/bin/cd script/optimize_tutorial/
python3 optimize_l1.py
chmod 777 collect_optimizationresult.sh
results=`./collect_optimizationresult.sh`

echo $results

