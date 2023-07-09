#!/bin/bash

lst=`ls`
req_lst=(${lst//,/ })

for req in ${req_lst[@]}
do
    elim_text=${req: 0:${#req}-4}
    dirname=r${elim_text: 4:${#elim_text}}
    echo $elim_text/$req
    echo $dirname/$req
    mkdir $dirname
    cp $req $dirname/$req
done

