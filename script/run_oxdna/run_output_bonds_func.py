#!/usr/bin/env python
# coding: utf-8

# In[60]:


import pandas as pd
import os
import subprocess as sp
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
import importlib
import config as cfg
importlib.reload(cfg)


# In[61]:


def run_output_bonds (target, output_oxdna_dir, oxdna_utils_dir = os.path.join(cfg.oxDNA_dir, "UTILS"),output_bonds = "output_bonds.py"):
    input_file = os.path.join(output_oxdna_dir, target + "_" + cfg.oxdna_input)
    trajectory_file = os.path.join(output_oxdna_dir, "{}_lastconf.dat".format(target))
    output_bonds_file = os.path.join(oxdna_utils_dir, output_bonds)
    executable = ["python2", output_bonds_file,  input_file, trajectory_file]
    print("EXEC : ",executable)
    output_filename = os.path.join(output_oxdna_dir, '{}_output_bonds.txt'.format(target))
    with open(output_filename, 'w') as fp:
        print("output bonds running 🖥",executable)
        sp.run(executable, stdout = fp, stderr = sp.DEVNULL)
    fp.close()
    return output_filename


# In[62]:


def get_bonds_data(filename):
    #output_bonds_data = pd.read_csv(filename,  sep = " ", header = None).drop(0).drop(columns = [10, 11, 12])
    #output_bonds_data.columns = ["id1", "id2", "FENE","BEXC", "STCK","NEXC", "HB", "CRSTCK", "CXSTCK",  "TOTAL"]
    #output_bonds_data = output_bonds_data.reset_index(drop = True)
    output_bonds_data = pd.read_csv(
        filename,
        sep = " ",
        names = ["id1","id2","FENE","BEXC","STCK","NEXC","HB","CRSTCK","CXSTCK","total","1","2","3"]
    ).iloc[2:,:10]
    output_bonds_data = output_bonds_data[:-1].astype("float").astype({"id1" : "int", "id2" : "int"})
    
    return output_bonds_data


# In[63]:


def get_lastconf_data(output_dir, target):
    lastconf_filename = os.path.join(output_dir, "{}_lastconf.dat".format(target))
    lastconf_data = pd.read_csv(
        lastconf_filename,
        sep = " ",
        header = None,
        names = ["position_rx", "position_ry", "position_rz",
                 "bb_versor_bx", "bb_versor_by", "bb_versor_bz",
                 "N_versor_nx", "N_versor_ny", "N_versor_nz",
                 "velocity_vx", "velocity_vy", "velocity_vz",
                 "angular_velocity_Lx", "angular_velocity_Ly", "angular_velocity_Lz"])
    lastconf_data = lastconf_data.drop(range(0,3)).reset_index(drop = True).astype(float)
    return lastconf_data


# In[64]:


def get_topology_data(output_dir, target): 
    topology = os.path.join(output_dir, "{}.top".format(target))
    topology_data = pd.read_csv(topology, sep = " ", names = ["strand", "nucleotide", "connection1", "connection2"])
    topology_data = topology_data.drop([0]).reset_index(drop = True)
    topology_data = topology_data.reset_index()
    topology_data = topology_data.rename(columns={"index" : "id1"}).astype({"strand" : int, "connection1" : int, "connection2" : int})
    expected_num_strands = topology_data["strand"].max()
    return topology_data, expected_num_strands


# In[65]:


def get_top_pos_data(data1, bonds_data):
    newdata = pd.merge(bonds_data, data1)
    drop_col = ["FENE", "BEXC", "STCK", "NEXC", "CRSTCK", "CXSTCK", "total",
                "connection1", "connection2", "bb_versor_bx",
                "bb_versor_by", "bb_versor_bz", "N_versor_nx", 
                "N_versor_ny", "N_versor_nz",                 
                "velocity_vx", "velocity_vy", "velocity_vz", 
                "angular_velocity_Lx", "angular_velocity_Ly", "angular_velocity_Lz"]
    
# ['id1', 'id2', 'FENE', 'BEXC', 'STCK', 'NEXC', 'HB', 'CRSTCK', 'CXSTCK',
#        'total', 'strand', 'nucleotide', 'connection1', 'connection2',
#        'position_rx', 'position_ry', 'position_rz', 'bb_versor_bx',
#        'bb_versor_by', 'bb_versor_bz', 'N_versor_nx', 'N_versor_ny',
#        'N_versor_nz', 'velocity_vx', 'velocity_vy', 'velocity_vz',
#        'angular_velocity_Lx', 'angular_velocity_Ly', 'angular_velocity_Lz']
    newdata = newdata.drop(drop_col, axis=1)
    
    
    return newdata


# In[66]:


def add_id2_strand(data, topology_data):
    data["id2_strand"] = " "
    for index, row in data.iterrows():
        id2_value = row["id2"]
        id2_strand = topology_data[ id2_value == topology_data["id1"]]["strand"].values[0]
        data.at[index, "id2_strand"] = id2_strand
    #display(data)
    return data


# In[67]:


def count_strands(top_pos_data, data_groupby):
    
    howmany_strands = pd.DataFrame(top_pos_data['id2_strand'].value_counts()).rename({"id2_strand" : "id2_strand_count"}, axis = 1)
    id2_strands_num = len(howmany_strands)
    return id2_strands_num


# In[68]:


def get_connected_strands_data(data):

    connected_strands = data[data["HB"] < 0.0][["strand", "id2_strand"]].values.tolist() 
    set_connected_strands = set( [ (a,b) for a,b in connected_strands ])#setにすると重複がなくなる
    
    #display(set_connected_strands)
    #「他のストランドのヌクレオチドと水素結合しているヌクレオチド」が存在するストランド番号の組み合わせを特定する。
    pair_list =list ( set_connected_strands)#結合があるストランドの組合せのリスト
    #pair_list = [pair for pair in pair_list if pair[0] != pair[1]]#同じもの同士で繋がっているものは含まない
    #もしpair_listが空であった場合、全てのストランドは繋がっていないことになるため、自動的にstrand = 1であるものが出力に選ばれる。

    connected_strand_table = [[]]
    

    #example: (1,2), (2,3), (3,5) , (4,6)→[1,2,3,5], [4,6]
    #pair_listのソート済みを前提とする
    #結合があるストランドのグループを特定する
    for pair in pair_list:
        for line in connected_strand_table:

            if (pair[0] in line) and (pair[1] not in line):
                line.append(pair[1])

            elif  (pair[1] in line) and (pair[0] not in line):
                line.append(pair[0])

            elif (pair[0] not in line) and (pair[1] not in line):
                connected_strand_table.append(list(pair))

    
    connected_strand_table = [x for x in connected_strand_table if x]#最初の空の列を除外
    connected_strand_table = [line for line in connected_strand_table if line[0] != line[1]]#同じストランド同士のつながりを除外
    connected_strand_table = list(map(list, set(map(tuple, connected_strand_table))))
    
    strand_size_table = [len(line) for line in connected_strand_table]
    #print("strand_size_table: ", strand_size_table)
    #print("connected_strand_table : ", connected_strand_table)

    df = pd.DataFrame(columns=['groups', 'group_size'])
    df['groups'] =pd.Series(connected_strand_table)
    df["group_size"] = np.array(strand_size_table)
    #print(df)
    return_data = pd.DataFrame(columns = data.columns)
    
    if df.size > 0:
        actual_num_strands = df["group_size"].max() 
        biggest_group_array = df.iloc[ df["group_size"].idxmax()]["groups"]
    else:
        #どれも繋がっておらずdfが空になっていた場合の対応
        biggest_group_array = [1]
        actual_num_strands = 1
    return_data = data[[a in biggest_group_array for a in data["strand"]]]    

    return return_data, actual_num_strands


# In[69]:


def create_connection_data(target, output_dir):
    output = run_output_bonds(target, output_dir)#ファイル名
    print("run_output_bonds() finished 🧬🖥 : ", output)
    output_bonds_data = get_bonds_data(output)
    
    lastconf_data = get_lastconf_data(output_dir, target)
    topology_data, expected_num_strands = get_topology_data(output_dir, target)
    newdata1 = pd.concat([topology_data, lastconf_data], axis = 1)
    newdata = get_top_pos_data(newdata1, output_bonds_data)
    data = add_id2_strand(newdata, topology_data)#何か発生？
    connected_data, actual_num_strands = get_connected_strands_data(data)#e77error
    return connected_data, expected_num_strands, actual_num_strands


# In[70]:


def main():
    #mainはデバッグのために直接実行した場合のみ起動する。
    connected_data, expected_num_strands, actual_num_strands = create_connection_data(target = "e0", output_dir = "./sim_result_peppercorn_result1_1")#test
    print(connected_data)
    print("expected number of strands : ", expected_num_strands)
    print("actual number of strands : ", actual_num_strands)


# In[71]:


if __name__ == "__main__":
    main()


# In[ ]:





# In[ ]:




