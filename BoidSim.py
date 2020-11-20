import numpy as np 
from matplotlib import animation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import timeit
import statistics as st

import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import boid_lib as bl

###############################################################################
# CLASSES AND FUCTIONS
###############################################################################

class Boid():
    
    def __init__(self, pos, vel, acc, ident, mass=1.0, b_type="boid"):
        self.pos = np.array(pos)
        self.vel = np.array(vel)
        self.acc = np.array(acc)
        self.id = ident
        self.mass = mass
        
        self.near_list = None
        self.near_dists = None
        self.sep_list = None
        self.sep_dists = None
        
        self.steer = np.zeros((3,1))
        
        self.type = b_type

def genFlock(num):
    l_lim_pos = np.array([0, 0, 0])
    u_lim_pos = np.array([500, 500, 500])
    l_lim_vel = np.array([-10, -10, -10])
    u_lim_vel = np.array([10, 10, 10])
    
    w_pos = u_lim_pos - l_lim_pos
    w_vel = u_lim_vel - l_lim_vel
    
    acc = np.array([0, 0, 0])
    flock = []
    for i in range(num):
        pos = l_lim_pos[:, np.newaxis] + np.random.random((3,1)) * w_pos[:, np.newaxis]
        vel = l_lim_vel[:, np.newaxis] + np.random.random((3,1)) * w_vel[:, np.newaxis]
        flock.append(Boid(pos, vel, acc, i))
    return flock    


###############################################################################
# MAIN
###############################################################################

num_boids = 1000
dt = 1

flock = genFlock(num_boids)

# l_lim_pos = np.array([800, 800, 800])
# u_lim_pos = np.array([900, 900, 900])
# l_lim_vel = np.array([-10, -10, -10])
# u_lim_vel = np.array([10, 10, 10])

# w_pos = u_lim_pos - l_lim_pos
# w_vel = u_lim_vel - l_lim_vel

# acc = np.array([0, 0, 0], dtype=np.double)

# for i in range(num_boids):
#     pos = l_lim_pos[:, np.newaxis] + np.random.random((3,1)) * w_pos[:, np.newaxis]
#     vel = l_lim_vel[:, np.newaxis] + np.random.random((3,1)) * w_vel[:, np.newaxis]
    
#     flock.append(bl.Boid(pos, vel, acc, i))

MAX_VEL = 25
MAX_ACC = 3
MAX_FORCE = 3*9.81
WALL_FORCE = MAX_FORCE

C_STR = 0.5
S_STR = 0.8
A_STR = 0.25

plot = True
timing = False

###############################################################################
# TIMING
###############################################################################

def timeTenSecs(flock, dt):
    time = 0
    
    while time < 10:
        bl.update(flock, dt, MAX_VEL, MAX_ACC, MAX_FORCE, C_STR, S_STR, A_STR)
        time += dt    
     
if timing:
    setup = "from __main__ import timeTenSecs, flock, dt"
    t_tensec = timeit.Timer(stmt = "timeTenSecs(flock,dt)", setup=setup).repeat(1,1)
    print(st.mean(t_tensec))

###############################################################################
# ANIMATION and PLOTTING
###############################################################################

def animate(frame):
    bl.update(flock, dt, MAX_VEL, MAX_ACC, MAX_FORCE, C_STR, S_STR, A_STR)
    positions = bl.concVecs(flock)
    scatter._offsets3d = (positions[0, :], positions[1, :], positions[2, :])


if plot:
    positions = bl.concVecs(flock)
    
    # Setup plot for animation
    plt.ion()
    fig = plt.figure()
    ax = axes3d.Axes3D(fig)
    ax.set_xlim3d(-200, 1200)
    ax.set_ylim3d(-200, 1200)
    ax.set_zlim3d(-200, 1200)
    
    # assign initial positions - ie create frame 1
    scatter = ax.scatter(positions[0, :], positions[1, :], positions[2, :],
                          marker='o', edgecolor='k', lw=0.5)
    
    # This must be performed in main body to work!
    anim = animation.FuncAnimation(fig, animate, frames=480, interval=50)
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=5000)
    anim.save("boids_3d_cyth.mp4", writer=writer)