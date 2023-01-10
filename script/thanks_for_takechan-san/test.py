import numpy as np
import matplotlib.pylab as plt

import scipy
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import vista_func as vf
import importlib
importlib.reload(vf)

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


def main():
    convexhull_volume(connection_data, target, output_dir)

if __name__ == '__main__':
  main()

