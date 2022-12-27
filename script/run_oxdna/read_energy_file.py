#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import os
import re


# In[1]:


# #test
# #[time (steps * dt)]	[potential energy]	[kinetic energy]	[total energy]
# output_folder = "sim_result_peppercorn_2022-08-27_04_33_17385170"
# target = "e0"
# energy_file = os.path.join(output_folder,"{}_energy.dat".format(target))
# with open(energy_file,"r") as f:
#     d = f.readlines()
#     data = []
#     for line in d:
#         data.append(re.findall(r"([-+]?\d+.?\d*)", line))

# energy_df = pd.DataFrame((data),columns=["time","potential_energy","kinetic_energy","total_energy"])
# energy_df.loc[:,"potential_energy"].astype(float).mean()


# In[8]:


def potential_energy_mean(output_folder,target):
    energy_file = os.path.join(output_folder,"{}_energy.dat".format(target))
    with open(energy_file,"r") as f:
        d = f.readlines()
        data = []
        for line in d:
            data.append(re.findall(r"([-+]?\d+.?\d*)", line))

    energy_df = pd.DataFrame((data),columns=["time","potential_energy","kinetic_energy","total_energy"])
    return energy_df.loc[:,"potential_energy"].astype(float).mean()


# In[2]:


# output_folder = "sim_result_peppercorn_2022-08-27_04_33_17385170"
# target = "e0"
# potential_energy_mean(output_folder,target)

