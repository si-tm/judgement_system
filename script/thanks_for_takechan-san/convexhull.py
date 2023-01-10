#!/usr/bin/env python
# coding: utf-8

# In[32]:


import scipy
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import vista_func as vf
import importlib
importlib.reload(vf)


# In[33]:



#def convexhull_volume(target, output_dir):
def convexhull_volume(connection_data, target, output_dir):
    
    points = vf.make_position_data_ndarray(connection_data)
    #display(points)
    hull = ConvexHull(points)
    x = np.array(points[:, 0])
    y = np.array(points[:, 1])
    z = np.array(points[:, 2])
    fig = plt.figure()
    ax = Axes3D(fig) 
    ax.scatter(x, y, z, color="#aa0000")
    
    ax.set_xlabel("a-value")
    ax.set_ylabel("b-value")
    ax.set_zlabel("L-value")
    
    figfilename = output_dir + "/" + target + ".png"
    for simplex in hull.simplices:

        ax.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2], "o-", color="#00aa00", ms=4, mew=0.5)
    fig.savefig(figfilename)
    return hull.volume


# In[34]:


def test():
    import run_output_bonds_func as robf
    # target = "e0"
    target = "../../input/pils/L1-entropyRT-1-22-39.pil"
    connection_data, num1, num2 = robf.create_connection_data(target, output_dir="output_dir/")
    # hull_volume = convexhull_volume(connection_data)
    # hull_volume = convexhull_volume(connection_data, output_dir="output_dir/")
    # print(hull_volume)


# In[35]:


def main():
    test()
    


# In[36]:


if __name__ == "__main__":
    main()


# e19: chimeraでは３つのストランドが絡み合っているように見えたが、output_bondsのデータを見る限りでは第3ストランドが繋がっていない

# In[ ]:




