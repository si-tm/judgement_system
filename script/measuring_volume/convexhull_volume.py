import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('../')
sys.path.append('.')
sys.path.append('measuring_volume/')
sys.path.append('common/')
import common.get_target_file as gtf
from mpl_toolkits.mplot3d import Axes3D
import measuring_volume.get_top_data as gtd
import measuring_volume.get_conf_data as gcd
from scipy.spatial import ConvexHull
# import k3d
import statistics
import math

# 間違ってる

def get_all_r(target_dir):
    conf_name = gtf.get_conf(target_dir)
    conf_f = open(conf_name, "r")
    col = 0


    x = []
    y = []
    z = []
    for l in conf_f:
        col += 1
        if col > 3:
            rx = float(l.split(" ")[0])
            ry = float(l.split(" ")[1])
            rz = float(l.split(" ")[2])
            x.append(rx)
            y.append(ry)
            z.append(rz)

    conf_f.close()
    return x, y, z

def translate_strand(points, box=100):
    # type of points is [[cm_pos_x, cm_pos_y, cm_pos_z]]
    # 2.translate(np.rint((s1.cm_pos - s2.cm_pos - diff1 + diff2) / self._box) * self._box)
    # diff1 = np.rint(s1.cm_pos / self._box ) * self._box

    new_points = []

    for index, p in enumerate(points):
        if index == 0:
            plst = np.array(p)
            first_diff = np.rint(plst/box) * box
            first_p = plst
            new_points.append(p)
        else:
            plst = np.array(p)
            tmp_diff = np.rint(plst/box) * box
            # print(tmp_diff)
            # p += np.rint((first_p - plst - first_diff + tmp_diff) / box) * box
            p += first_diff - tmp_diff
            # print(np.rint((- plst + tmp_diff) / box) * box)
            new_points.append(p)
    
    return new_points



def get_r(target_dir, strands):

    target_strands = list(strands)
    conf_dic = gcd.get_conf_data(target_dir)


    x = []
    y = []
    z = []

    for particle_id in target_strands:
        x.append(conf_dic[particle_id]["rx"])
        y.append(conf_dic[particle_id]["ry"])
        z.append(conf_dic[particle_id]["rz"])
    
    return x, y, z

def convexhull_volume(x, y, z):
    points = []
    for i in range(len(x)):
        lst = [x[i], y[i], z[i]]
        points.append(lst)
    
    print(points)
    # 追加
    points = translate_strand(points, box=100)
    points = np.array(points)
    
    hull = ConvexHull(points)



    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection="3d")

    # # Plot defining corner points
    # ax.plot(points.T[0], points.T[1], points.T[2], "ko")
    # # ax.plot(x, y, z, "ko")

    # # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
    # for s in hull.simplices:
    #     s = np.append(s, s[0])  # Here we cycle back to the first coordinate
    #     ax.plot(points[s, 0], points[s, 1], points[s, 2], "r-")
    #     # ax.plot(x, y, z, "r-")

    # # Make axis label
    # for i in ["x", "y", "z"]:
    #     eval("ax.set_{:s}label('{:s}')".format(i, i))

    # plt.show()

    
    return hull.volume
    
def plot(X, Y, Z, target_dir):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(X, Y, Z)
    plt.show()

def convexhull_volume_all_strands(target_dir):
    strands2particle, particle2strand = gtd.get_particle_strands_data(target_dir)

    volumes = {}

    for strand in strands2particle:
        x, y, z = get_r(target_dir, strands2particle[strand])
        volumes[strand] = convexhull_volume(x, y, z)
        # plot(x, y, z, target_dir)
    
    mean_volume = 0.0
    num_of_strands = float(len(strands2particle))

    for strand in volumes:
        # print(volumes[strand])
        mean_volume += volumes[strand]
    
    mean_volume /= num_of_strands
    print(target_dir + " : mean volume is " + str(mean_volume))
    print(target_dir + " : number of strands is " + str(num_of_strands))
    return mean_volume

def convexhull_volume_all_strands_meandev(target_dir):
    strands2particle, particle2strand = gtd.get_particle_strands_data(target_dir)

    volumes = {}

    for strand in strands2particle:
        x, y, z = get_r(target_dir, strands2particle[strand])
        volumes[strand] = convexhull_volume(x, y, z)
        #         print(x, y, z)
        # plot(x, y, z, target_dir)
    
    print(volumes)
    
    mean_volume = 0.0
    num_of_strands = float(len(strands2particle))

    for strand in volumes:
        mean_volume += volumes[strand]
    
    mean_volume /= num_of_strands
    dev_volume = statistics.pstdev(volumes)
    print(num_of_strands)

    return mean_volume, dev_volume

def test():
    target_dir="../../input/results/oxdna_ked/seqA/A4/test_a4_200000_1"
    convexhull_volume_all_strands(target_dir)

def main():
    if len(sys.argv) != 2:
        print("usage : python convexhull_volume.py [target directory name]")
    else:
        convexhull_volume_all_strands(sys.argv[1])

if __name__ == '__main__':
    # test()
    main()
