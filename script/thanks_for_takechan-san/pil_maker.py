#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pprint
import pickle
import glob
import os
import sys
import qdpy

sys.path.append("../scripts")
import submitPepperCorn as sp


# In[2]:


class Custom_unpickler(pickle._Unpickler):

    current_module = {"Individual": "qdpy.phenotype", "Fitness": "qdpy.phenotype"}
    
    def find_class(self, module, name):
        if name in Custom_unpickler.current_module:
            module = Custom_unpickler.current_module[name] #backward compatibility
        sys.audit('pickle.find_class', module, name)
        if self.proto < 3 and self.fix_imports:
            if (module, name) in _compat_pickle.NAME_MAPPING:
                module, name = _compat_pickle.NAME_MAPPING[(module, name)]
            elif module in _compat_pickle.IMPORT_MAPPING:
                module = _compat_pickle.IMPORT_MAPPING[module]
        __import__(module, level=0)
        if self.proto >= 4:
            return pickle._getattribute(sys.modules[module], name)[0]
        else:
            return getattr(sys.modules[module], name)


# In[3]:


def read_pickle(path):
    with open(path,"rb") as f:
        res = Custom_unpickler(f,fix_imports=True, encoding="ASCII", errors="strict").load()
        return res


# In[4]:


def indiv_submit(indiv_list,count_start,result_dir):
    counter = count_start
    for indiv in indiv_list:
        
        number = str(counter).zfill(4)
        
        evalfile = os.path.join(result_dir, "{}_test.pil".format(number))
        outputfile = os.path.join(result_dir,"{}_result.pil".format(number))
        logfile = os.path.join(result_dir,"{}_pepper.log".format(number))
        
        lst_indiv = [(i,x) for i, x in enumerate(indiv) if x > 0.0]
        system = sp.generateDNASystem(lst_indiv)
        
        sp.submitSystem(system, evalName = evalfile, outputFile = outputfile, logFile = logfile)
        
        counter = counter + 1


# In[5]:


def make_pil(res_container,result_dir,max_indiv = 30):
    #max_indivで指定した数だとオーバーする場合、どうしたらよいだろうか？
    
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
        
    indiv_list = list(res_container)
    indiv_list = sorted(indiv_list,key = lambda indiv: indiv.fitness[0], reverse = True)
    
    
    
    indiv_submit(indiv_list[:max_indiv],1,result_dir)#先頭から10個
    indiv_submit(indiv_list[-10:],1+max_indiv,result_dir)#最後から10個


# In[6]:


#use to reset pilfiles and logs!
def trash_pils():
    from send2trash import send2trash

    pils = files = glob.glob("../results/peppercorn*/final*/")
    for path in pils:
        send2trash(path)


# In[7]:


# test
# path = "../results/peppercorn30x400x2000-L1-meanStruct1x30-nbActive2x7-log10ReactionCount1.0x5.0-length1x550-Grid3x50x55--Greycode-65536x64/final_20220727093209.p"
# res = read_pickle(path)
# res_container = res["container"]
# res_container
#make_pil(res_container,result_dir="test")


# In[8]:


def make_pil_for_all():
    lst = glob.glob("../results/peppercorn*/final_*.p")

    for path in lst:
        res = read_pickle(path)
        result_dir = path.replace(".p","")
        print(result_dir)

        res_container = res["container"]
        make_pil(res_container,result_dir)


# In[9]:


def main():
    make_pil_for_all()
    print("completed\n")


# In[10]:


if __name__ == "__main__":
    main()

