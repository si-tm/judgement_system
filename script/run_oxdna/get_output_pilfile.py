#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import sys
import shutil
import glob
import datetime
import subprocess as sp
import pathlib

import config as cfg
import make_mean as mm
import importlib
importlib.reload(cfg)
importlib.reload(mm)


# In[2]:


#条件ごとのpeppercorn実行結果のディレクトリ一覧を"パス名で"得る
def get_result_dir_path(all_results_dir):
    peppercorn_results = glob.glob("{}/peppercorn*".format(all_results_dir))
    return peppercorn_results
    #output : list of "resluts/(peppercorn output libraries)"


# In[3]:


def cut_parent_dir(filepath, parent_dir):
    newfilepath = filepath.replace(parent_dir+"/", "")
    return newfilepath
#cut_parent_dir("../results/output.pil", "../results")


# In[4]:


def get_lastconf(output_folder):
    p_temp = pathlib.Path(output_folder).glob('*_lastconf.dat')
    lastconf = [os.path.basename(str(p)) for p in list(p_temp)]
    return lastconf


# In[5]:


def run_main_for_eachfiles(results_dir_path, empty_dir_list):
    
    d = {' ' :  '_', '.' :  '', ':' : '_'}
    tbl = str.maketrans(d)

    #empty_dir_list = []#結果ファイルがなかったディレクトリを記録する。

    result_files = glob.glob("{}/output*.pil".format(results_dir_path))
    
    
    if len(result_files) == 0:#result_filesがemptyなら、それを記録する
        print("RESULT FILE NOT FOUND\n")
        empty_dir_list.append(results_dir_path)

    else:#result_filesがemptyじゃない時はシミュレーションする
        for result_file_long in result_files:

            datetime_sim = str(datetime.datetime.now()).translate(tbl)#時刻の文字列。ファイル名の作成に使う
            output_folder = ("./sim_result_peppercorn_"+datetime_sim)# oxDNAの結果が入るディレクトリ
            #if not os.path.exists(output_folder):#毎回時刻は違うのでif notしなくて良い
            os.mkdir(output_folder)#oxDNA結果が入るディレクトリを作る
            #print( "oxdna実行結果出力先：", output_folder, "\n")

            result_file_long_name = os.path.basename(result_file_long)
            print("result_file_long_name : ", result_file_long_name, "\n")
            result_file = shutil.copy(result_file_long, os.path.join(output_folder,result_file_long_name))#pil結果ファイルをoxDNA結果ディレクトリにコピー
            print("結果ファイル：", result_file, "\n")

            exeprogram = "python3"
            #exefile = "./main.py"
            exefile = "./main.py"
            executable = [exeprogram, exefile, result_file, output_folder]
            print("実行コマンド：", executable, "\n")
            #python3 argv[0]:main.py argv[1]:<result pilfile> argv[2]:<output folder>
            logfile = os.path.join(output_folder, datetime_sim+"_log.txt")
            #print("ログ出力先: ",logfile, "\n")
            print("-----------------sim : {} start------------------\n".format(output_folder))
            with open(logfile,"w") as log:
                sp.run(executable, stdout=log, stderr=log, text=True)#sp.STDOUTの場合は標準出力に出る
                print(result_file)
                log.close()
                print("oxdna end : ",output_folder,file = sys.stdout)
            oxdna_results = get_lastconf(output_folder)
            print("oxdna_results : ", oxdna_results)
            print("-------------------sim : {} end------------------\n".format(output_folder))            
        

    #print(empty_dir_list)
    
    
        


# In[6]:


def sim_all_results_dir(all_results_dir = cfg.results_dir):
    
    results_path_list= get_result_dir_path(all_results_dir)
    #print("results_path_list :" , results_path_list)
    empty_dir_list = []
    for results_dir_path in results_path_list:
    
        results_dir_name = cut_parent_dir(results_dir_path, all_results_dir)
        #results/<ある条件でのpeppercorn実行結果ディレクトリ> の１つ
        print("ライブラリパス：", results_dir_path, "\n")
        print("ライブラリ名：", results_dir_name, "\n")
        run_main_for_eachfiles(results_dir_path,empty_dir_list)
        
        empty_dir_log = "empty_dir_log.txt"      
        with open(empty_dir_log,"w") as elog:
            for line in empty_dir_list:
                elog.write(str(line))
                elog.write("\n")
            elog.close()
            
    search_dir_name = "sim_result_peppercorn*"
    mm.make_all_mean_file(search_dir_name)
    print("mean file created\n",file=sys.stdout)
        
    print("🎉🤗All simuration was finished!🤗🎉")


# In[7]:


def main():
    d = {' ' :  '_', '.' :  '', ':' : '_'}
    tbl = str.maketrans(d)
    datetimestr = str(datetime.datetime.now()).translate(tbl)
    os.sys.stdout = open('simlog_{}.txt'.format(datetimestr), 'w')
    sim_all_results_dir()
    #def get_sim_resultlist(output_folder):


# In[8]:


if __name__ == "__main__":
    main()


# 

# In[9]:


# exeprogram = "python3"
# exefile = "./main.py"
# result_file = "/Users/takepy/takeoxdna/kakenhievolvedna2/results/peppercorn30x400x2000-L2-meanStruct1x30-nbActive2x7-log10ReactionCount1.0x5.0-length1x550-Grid3x50x55--SparseBinaryRandom-300000x40-1x2-7/outputPepperCorn20220719191756_59310139722517910325659710441010282641_0.pil"
# output_folder = "test"
# executable = [exeprogram, exefile, result_file, output_folder]
# print(exeprogram,exefile,result_file,output_folder)
# sp.run(executable)


# In[10]:


# import run_output_bonds_func as robf
# import importlib
# importlib.reload(robf)
# output_dir = "test"
# target  = "e0"

# path = robf.run_output_bonds(target, output_oxdna_dir = output_dir, oxdna_utils_dir = os.path.join(cfg.oxDNA_dir, "UTILS"),output_bonds = "output_bonds.py")

# with open(path,"r")as f:
#     print(len(f.readlines()))

# #robf.get_lastconf_data(output_dir, target)


# In[11]:


# import subprocess as sp
# output_filename = "test.txt"
# which = sp.run(["which","python2"],capture_output=True, text=True).stdout.replace("\n","")
# print(which)
# with open(output_filename, 'w') as fp:
#     print("output bonds running 🖥",executable)
#     #sp.run([which,"../../oxDNA/UTILS/output_bonds.py","test/e0_input_relax_1e5","test/e0_lastconf.dat"],stdout = fp,stderr = fp)
#     sp.run(["pyenv","global","2.7.18","3.9.0"])#仮想環境を同時に導入すべき
#     sp.run(["/Users/takepy/.pyenv/shims/python2","../../oxDNA/UTILS/output_bonds.py","test/e0_input_relax_1e5","test/e0_lastconf.dat"],stdout = fp,stderr = fp)


# In[12]:


# sp.run(["pyenv","global","2.7.18","3.9.0"])
# sp.run(["pyenv","versions"])


# In[ ]:




