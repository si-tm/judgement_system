#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import pyvista as pv
import os
import run_output_bonds_func as robf
import get_coordinate as gc
import importlib
import matplotlib.pyplot as plt
importlib.reload(robf)
importlib.reload(gc)


# In[9]:


def make_position_data_ndarray(connection_data):
    position_data = gc.get_coordinate_data(connection_data)
    print("🍌")
    return position_data.values


# In[10]:


def plot_coordinate_vista(connection_data, target, output_dir):
    point_cloud = make_position_data_ndarray(connection_data)
    #print(point_cloud)
    pdata = pv.PolyData(point_cloud)
    #display(pdata)#N pointsの値, int
    n_points = pdata.n_points
    #pdata['orig_sphere'] = np.arange(n_points)
    sphere = pv.Sphere(radius=0.2, phi_resolution=10, theta_resolution=10)
    #display(sphere)
    # create many spheres from the point cloud
    pc = pdata.glyph(scale=False, geom=sphere)
    #display(pc)
    pc.plot(cmap='Reds')
    
    plt.savefig(os.path.join(output_dir, "{}.png".format(target)))


# In[11]:


def main():
    target = "e0"
    connection_data = robf.create_connection_data(target, output_dir = "output_oxDNA")
    plot_coordinate_vista(connection_data, target, output_dir)


# In[12]:


if __name__ == "__main__":
    main()


# In[ ]:




