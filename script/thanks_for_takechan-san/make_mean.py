#!/usr/bin/env python
# coding: utf-8

# In[36]:


import glob
import os
import pandas as pd


# In[94]:


def make_csv(datafilename,meanfilename,en_df):
    try:
        en_df.to_csv(datafilename,header =True, index = False)
        en_df.drop("id",axis=1).mean().to_csv(meanfilename,header=False)
        print(en_df)
        print(en_df.drop("id",axis=1).mean())
    except:
        with open(datafilename,"w") as log1:
            print("PIL FILE NOT FOUND\n")
        with open(meanfilename,"w") as log2:
            print("PIL FILE NOT FOUND\n")


# In[95]:


def make_each_mean_file(result_dirs):
    for folder in result_dirs:
        sizelogs = glob.glob(os.path.join(folder,"*_sizelog.txt"))
        en_df = pd.DataFrame([])
        for sizelog in sizelogs:
            with open(sizelog,"r") as f:
                df = pd.read_csv(sizelog)
                en_df = pd.concat([en_df,df],axis = 0).reset_index(drop=True)

        print(folder,"\n")
        datafilename = os.path.join(folder,"data.csv")
        meanfilename = os.path.join(folder,"mean.csv")

        make_csv(datafilename,meanfilename,en_df)


# In[96]:


def make_all_mean_file(search_dir_name):
    result_dirs = glob.glob(search_dir_name)
    make_each_mean_file(result_dirs)


# In[97]:


def main():
    search_dir_name = "sim_result_peppercorn*"
    make_all_mean_file(search_dir_name)


# In[98]:


if __name__ == "__main__":
    main()


# In[ ]:




