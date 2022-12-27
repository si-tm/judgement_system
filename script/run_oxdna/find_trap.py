#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd


# In[2]:


#target = "e_test"


# In[3]:


def make_domain_list(linelist, length_tup):
    index_domains = []
    
    count_start = 0
    
    for block in linelist:
        #print("block : ", block)
    
        for b in block:
            count_end = count_start + length_tup[b[0][0]]
            index_domains.append([b[0],b[1],count_start,count_end-1])
            count_start = count_end
        
        #index_domain.append( index_domain1 + index_domain2 + [block])
        print(index_domains)
    return index_domains

# In[4]:


def make_domain_df(lst):
    df2 = pd.DataFrame(lst, columns =["domain", "num", "start_index", "end_index"])
    return df2


# In[5]:


def drop_not_connected(df):
    newdf = df[df['num'] >= 0]
    return newdf


# In[6]:


def make_groups(df):
    pairs = df.groupby(["num"])
    return pairs


# In[7]:


def count_domains(df):
    df_pairs = pd.DataFrame(df["num"].value_counts())
    return len(df_pairs)


# In[8]:

def make_external(df, target, output_dir):
    try : 
        #output_path = "traps/external_{}.conf".format(target)
        output_path = os.path.join(output_dir, "{}_external.conf".format(target))
        external_file = open(output_path, "w")
        domains_num = count_domains(df)
        pairs = make_groups(df)
        #print(pairs.groups)
        
        for x in range (0, domains_num):  
            #print("roop", x)
            group = pairs.get_group(x)
            #print("group:\n", group, "\n")
            start_index_min = group["start_index"].min()
            end_index_max =  group["end_index"].max()
            #print("index : ", start_index_min, " " , end_index_max)
            lines = ["{\n","type = mutual_trap\n","particle = {}\n".format(start_index_min),"ref_particle = {}\n".format(end_index_max),"stiff = 1.\n","r0 = 1.2\n","}\n""{\n", "type = mutual_trap\n","particle = {}\n".format(end_index_max),"ref_particle = {}\n".format(start_index_min),"stiff = 1.\n","r0 = 1.2\n","}\n"]
            external_file.writelines(lines)
        external_file.close()
        print(output_path, " was created\n")
    except KeyError:
        print("group {} does not exist".format(x))


# In[9]:


def make_trap(linelist, length_tup, target, output_dir):
    #print(linelist, "\n")
    domain_list = make_domain_list(linelist, length_tup)
    #print("domain_list : ", domain_list, "\n")
    domain_df = make_domain_df(domain_list)
    #print("domain_df : \n", domain_df)
    connected_df = drop_not_connected(domain_df)
    #print("connected_df : \n", connected_df)
    make_external(connected_df, target, output_dir)
    return domain_list


# In[ ]:





# In[ ]:




