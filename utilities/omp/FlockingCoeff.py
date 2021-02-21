import numpy as np
import pandas as pd
import math as m
import statistics as st
import time
import numba as nb
import configparser


# USER SETTINGS IMPORT ==================================================================
fps = 30
skip_factor = 1 # Only change this if you use small DT, skip_factor 1 does not skip any frames

config = configparser.ConfigParser()
config.read("../../boid_data/settings.cfg")

close_alert = int(config["Neighbour Finding"]["too_close_radius"])
duration = float(config["Time"]["duration"])
dt = float(config["Time"]["timestep"])
wall_bounds = [int(config["Bounds"]["wall_bound_lower"]), int(config["Bounds"]["wall_bound_upper"])]
boids = int(config["Program"]["number_boids"])
threads = int(config["Program"]["number_procs"])
# =======================================================================================

def skip_diag_strided(A):
    m = A.shape[0]
    strided = np.lib.stride_tricks.as_strided
    s0,s1 = A.strides
    return strided(A.ravel()[1:], shape=(m-1,m), strides=(s0+s1,s1)).reshape(m,-1)


filename = "B" + str(boids) + "Thr" + str(threads)
frame_count = int(duration//(dt*skip_factor))
data = pd.read_csv("../../boid_data/"+filename+".csv", sep=',')
bound_width = wall_bounds[1] - wall_bounds[0]


frames_m = [i for i in range(int(duration/dt))]
times = [f/max(frames_m) * duration for f in frames_m]
mave = []
pavex = []
pavey = []
pavez = []
for f in frames_m:
    pavex.append((data["X"][f*boids:f*boids+boids]).to_numpy())
    pavey.append((data["Y"][f*boids:f*boids+boids]).to_numpy())
    pavez.append((data["Z"][f*boids:f*boids+boids]).to_numpy())

pavex = np.array(pavex)
pavey = np.array(pavey)
pavez = np.array(pavez)

dists = np.zeros((boids, boids))

@nb.njit(parallel=True)
def FindDists(boids, frame, px, py, pz):
    dists = np.zeros((boids, boids))

    for i in range(boids):
        for j in range(boids):
            dists[i,j] = np.sqrt((px[frame][i]-px[frame][j])**2 + (py[frame][i]-py[frame][j])**2 + (pz[frame][i]-pz[frame][j])**2)
            
    return dists


num = 8
dmax = (2 * bound_width) * ((3 * 0.64) / (4 * m.pi * boids))**(1/3)
start = time.time()
file = open("../../boid_data/FC_B{0}Thr{1}.csv".format(boids, threads), "w")
for frame in frames_m:
    dists = FindDists(boids, frame, pavex, pavey, pavez)
    if (frame % 100) == 0:
        percent = 100* frame / frame_count
        print("{}% Processed.".format(int(percent)))
    dists = skip_diag_strided(dists)
    
    top = np.argpartition(dists, num)[:, :num]
    smallest = dists[np.arange(dists.shape[0])[:, None], top]
    ave = np.mean(smallest, -1).mean()
    
    f = (ave - dmax) / (close_alert - dmax)
    file.write(str(frame) + "," + str(f) + "\n")

file.close()