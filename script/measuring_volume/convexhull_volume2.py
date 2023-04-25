import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
sys.path.append('.')
sys.path.append('measuring_volume/')
sys.path.append('common/')
import common.get_target_file as gtf
from mpl_toolkits.mplot3d import Axes3D
import get_top_data as gtd
# import get_conf_data as gcd
from scipy.spatial import ConvexHull
# import k3d
import statistics
import math
import subprocess


def get_r(particle2r, strands):

    target_strands = list(strands)

    x = []
    y = []
    z = []

    for particle_id in target_strands:
        x.append(particle2r[particle_id][0])
        y.append(particle2r[particle_id][1])
        z.append(particle2r[particle_id][2])
    
    return x, y, z

def plot_points(points):
    points = np.array(points)
    hull = ConvexHull(points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot defining corner points
    ax.plot(points.T[0], points.T[1], points.T[2], "ko")
    # ax.plot(x, y, z, "ko")

    # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
    for s in hull.simplices:
        s = np.append(s, s[0])  # Here we cycle back to the first coordinate
        ax.plot(points[s, 0], points[s, 1], points[s, 2], "r-")


def convexhull_volume(x, y, z):
    points = []
    for i in range(len(x)):
        lst = [x[i], y[i], z[i]]
        points.append(lst)
        
    hull = ConvexHull(points)

    plot_points(points)

    return hull.volume

def execute_traj2r(target_dir):
    traj = gtf.get_conf(target_dir)
    top = gtf.get_top(target_dir)
    subprocess.run(["measuring_volume/traj2r.py", "xyz", traj, top])

# 非同期処理？？
def get_particle2r(target_dir):
    # execute traj2r.py in target_dir
    # get result file 
    execute_traj2r(target_dir)
    result_r_filename = gtf.get_conf(target_dir) + ".rxyz"
    # read points
    f = open(result_r_filename, "r")
    new_points = []

    for i, l in enumerate(f):
        if i == 0:
            continue
        tmp_point = []
        for r in l[:-2].split(' '):
            tmp_point.append(float(r))
        new_points.append(tmp_point)

    # make particle2r
    particle2r = {}
    for i, p in enumerate(new_points):
        particle2r[i] = p
    
    return particle2r


def tmp_convexhull_volume_all_strands_meandev2(target_dir):
    strands2particle, particle2strand = gtd.make_initial_strands_data(target_dir)
    particle2r = get_particle2r(target_dir)
    volumes = {}

    strands2particle, particle2strand = gtd.get_particle_strands_data(target_dir)
    print(len(strands2particle))


    for strand in strands2particle:
        x, y, z = get_r(particle2r, strands2particle[strand])
        volumes[strand] = convexhull_volume(x, y, z)
        #         print(x, y, z)
        # plot(x, y, z, target_dir)
        # plot_points(volumes[strand])
    
    print(volumes)
    
    mean_volume = 0.0
    num_of_strands = float(len(strands2particle))

    for strand in volumes:
        mean_volume += volumes[strand]
    
    mean_volume /= num_of_strands
    dev_volume = statistics.pstdev(volumes)
    # print(num_of_strands)

    return mean_volume, dev_volume


def test():
    target_dir="../input/results/oxdna_ked/seqA/A4/test_a4_200000_1"
    mean_volume, dev_volume = tmp_convexhull_volume_all_strands_meandev2(target_dir)
    print("average : ", mean_volume)
    print("deviation : ", dev_volume)
    

def main():
    if len(sys.argv) != 2:
        print("usage : python convexhull_volume.py [target directory name]")
    else:
        tmp_convexhull_volume_all_strands_meandev2(sys.argv[1])

if __name__ == '__main__':
    test()
    # main()
