from matplotlib import animation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import timeit
import statistics as st

import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np

class Boid():
    
    def __init__(self, pos, vel, acc, ident, mass=1.0, b_type="boid"):
        self.pos = np.array(pos, ndmin=2)
        self.vel = np.array(vel, ndmin=2)
        self.acc = np.array(acc, ndmin=2)
        self.id = ident
        self.mass = mass
        
        self.near_list = []
        self.near_dists = []
        self.close_list = []
        self.close_dists = []
        
        self.fc = np.zeros((3,1))
        self.fs = np.zeros((3,1))
        self.fa = np.zeros((3,1))
        
        self.steer = np.zeros((3,1))
        
        self.type = b_type
        
        
def findNear(boid, flock):
    exclusive_flock = flock.copy()
    exclusive_flock.pop(flock.index(boid))
    
    near_alert = 150
    sep_alert = 10
    near_list = []
    near_dists = []
    sep_list = []
    sep_dists = []
    
    for b in exclusive_flock:
        pos_diff = b.pos - boid.pos
        sqr_disp = pos_diff * pos_diff
        dist = np.sqrt(np.sum(sqr_disp, axis=0))
        
        if dist < near_alert:
            near_list.append(b.id)
            near_dists.append(pos_diff)
        
        if dist < sep_alert:
            sep_list.append(b.id)
            sep_dists.append(pos_diff)
            
    boid.near_list = near_list
    boid.near_dists = near_dists
    boid.close_list = sep_list
    boid.close_dists = sep_dists
    
def bulkFindNeighb(pos, flock):
    near_alert = 300
    close_alert = 10
    
    seps = pos[:, np.newaxis, :] - pos[:, :, np.newaxis]
    disps = np.sum(seps * seps, axis=0)
    dists = np.sqrt(disps)

    near_boids = dists < near_alert
    near_boids = np.logical_and(dists < near_alert, dists > close_alert)
    near_dists = np.copy(seps)
    far_boids = np.logical_not(near_boids)
    near_dists[0, :, ::1][far_boids] = 0
    near_dists[1, :, ::1][far_boids] = 0
    near_dists[2, :, ::1][far_boids] = 0
    
    close_boids = dists < close_alert
    close_dists = np.copy(seps)
    nclose_boids = np.logical_not(close_boids)
    close_dists[0, :, ::1][nclose_boids] = 0
    close_dists[1, :, ::1][nclose_boids] = 0
    close_dists[2, :, ::1][nclose_boids] = 0

    bulk_nl = np.argwhere(near_dists[0,:,:] != 0)
    bulk_nd = near_dists[:, near_dists[0]!=0]
    bulk_cl = np.argwhere(close_dists[0,:,:] != 0)
    bulk_cd = close_dists[:, close_dists[0]!=0]
    
    
    num_appd_n = 0
    for [i, j] in bulk_nl:
        flock[i].near_list.append(j)
        flock[i].near_dists.append(np.reshape(bulk_nd[:,num_appd_n], (3,1)))
        num_appd_n += 1

    num_appd_c = 0
    for [k, l] in bulk_cl:
        flock[k].close_list.append(l)
        flock[k].close_dists.append(np.reshape(bulk_cd[:,num_appd_c], (3,1)))
        num_appd_c += 1
        
    

def cohereForce(boid, flock):
    near_indicies = boid.near_list
    pos_sum = np.zeros((3,1))
    m_sum = 0

    if len(near_indicies) == 0:
        return
    
    for indx in near_indicies:
        m = flock[indx].mass
        m_sum += m
        pos_sum = np.add(pos_sum, (1/m)*flock[indx].pos)

    target = np.subtract((1/m_sum)*pos_sum, boid.pos)
    #target = (1/count) * pos_sum 
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    steer = np.subtract(target, boid.vel) * C_STR
    boid.fc += steer
    boid.steer += steer


def sepForce(boid, flock):
    near_indicies = boid.close_list
    pos_sum = np.zeros((3,1))
    m_sum = 0
    
    if len(near_indicies) == 0:
        return
    
    for indx in near_indicies:
        sep = np.subtract(boid.pos, flock[indx].pos)
        m = flock[indx].mass
        m_sum += m
        pos_sum = np.add(pos_sum, m*sep)
    
    target = pos_sum * (1/m_sum)
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    steer = np.subtract(target, boid.vel) * S_STR
    boid.fs += steer
    boid.steer += steer
    
    
def alignForce(boid, flock):
    near_indicies = boid.near_list
    vel_sum = np.zeros((3,1))
    m_sum = 0
    
    if len(near_indicies) == 0:
        return
    
    for indx in near_indicies:
        m = flock[indx].mass
        m_sum += m
        vel_sum = np.add(vel_sum, m*flock[indx].vel)
    
    target = (1/m_sum) * vel_sum
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    steer = np.subtract(target, boid.vel) * A_STR
    boid.fa += steer
    boid.steer += steer

def forceCalcs(boid, flock):
    near_indicies = boid.near_list
    sep_indicies = boid.close_list
    pos_sum = np.zeros((3,1))
    vel_sum = np.zeros((3,1))
    m_sum_n = 0
    m_sum_s = 0
    
    if len(near_indicies) == 0:
        return
    
    for indx in near_indicies:
        m = flock[indx].mass
        m_sum_n += m
        pos_sum = np.add(pos_sum, (1/m)*flock[indx].pos)
        vel_sum = np.add(vel_sum, m*flock[indx].vel)
    
    target = np.subtract((1/m_sum_n)*pos_sum, boid.pos)
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    steer = np.subtract(target, boid.vel) * C_STR
    boid.fc += steer
    
    target = (1/m_sum_n) * vel_sum
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    steer = np.subtract(vel_sum, boid.vel) * A_STR
    boid.fa += steer
    
    if len(sep_indicies) == 0:
        boid.steer += steer
        return 
    
    pos_sum = np.zeros((3,1))
    
    for indx in sep_indicies:
        sep = np.subtract(boid.pos, flock[indx].pos)
        m = flock[indx].mass
        m_sum_s += m
        pos_sum = np.add(pos_sum, m*sep)
    
    target = pos_sum * (1/m_sum_s)
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    steer = np.subtract(target, boid.vel) * S_STR
    boid.fs += steer
    

def updatePosEuler(boid, flock, dt):
    f_val = boid.fc + boid.fs + boid.fa
    f_wall = wallForce(boid)
    
    f_sum = limitVec(f_val, MAX_FORCE)
    f_sum = np.add(f_val, f_wall)
    
    m = boid.mass
    
    acc = (1/m) * f_sum
    boid.acc = limitVec(acc, MAX_ACC)
    vel = np.add(boid.vel, boid.acc * dt)
    boid.vel = limitVec(vel, MAX_VEL)
    boid.pos = np.add(boid.pos, boid.vel * dt)
    reset(boid)

    
def reset(boid):
    boid.steer = np.zeros((3,1))
    boid.fc = np.zeros((3,1))
    boid.fs = np.zeros((3,1))
    boid.fa = np.zeros((3,1))
    boid.near_list = []
    boid.near_dists = []
    boid.close_list = []
    boid.close_dists = []

        

#==============================================================================    
    
def gen_flock(num):
    l_lim_pos = np.array([0, 0, 0])
    u_lim_pos = np.array([250, 250, 250])
    l_lim_vel = np.array([-10, -10, -10])
    u_lim_vel = np.array([10, 10, 10])
    
    w_pos = u_lim_pos - l_lim_pos
    w_vel = u_lim_vel - l_lim_vel
    
    acc = np.zeros((3,1))
    flock = []
    for i in range(num):
        pos = l_lim_pos[:, np.newaxis] + np.random.random((3,1)) * w_pos[:, np.newaxis]
        vel = l_lim_vel[:, np.newaxis] + np.random.random((3,1)) * w_vel[:, np.newaxis]
        flock.append(Boid(pos, vel, acc, i))
    return flock        


def wallForce(boid):
    f_wall = np.zeros((3,1))
    
    if np.any(boid.pos > 1000):
        out_bound_pos = -(boid.pos > 800).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_pos * np.abs(boid.pos))
    if np.any(boid.pos < 0):
        out_bound_neg = (boid.pos < 0).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_neg * np.abs(boid.pos))
        
    return f_wall


def limitVec(vec, limit):
    vec_mag = np.linalg.norm(vec)
    if abs(vec_mag) > limit:
        return (limit / vec_mag) * vec
    else:
        return vec


def update(flock, dt):
    positions = concVecs(flock)
    bulkFindNeighb(positions, flock)
    for b in flock:
        #findNear(b, flock)
        cohereForce(b, flock)
        sepForce(b, flock)
        alignForce(b, flock)
        #forceCalcs(b, flock)
        updatePosEuler(b, flock, dt)
        if flock.index(b) == 0:
            print(b.vel, "\n")
        

def update1(flock, dt):
    #positions = concVecs(flock)
    #bulkFindNeighb(positions, flock)
    for b in flock:
        findNear(b, flock)
        cohereForce(b, flock)
        sepForce(b, flock)
        alignForce(b, flock)
        #forceCalcs(b, flock)
        updatePosEuler(b, flock, dt)
        

def concVecs(flock):
    for b in flock:
        if flock.index(b) == 0:
            ini_vec = flock[0].pos
        else:
            ini_vec = np.concatenate((ini_vec, b.pos), axis=1)
    return ini_vec

# def concVecs(flock):
#     x_list = []
#     y_list = []
#     z_list = []
#     for b in flock:
#         x_list.append(b.pos[0,0])
#         y_list.append(b.pos[1,0])
#         z_list.append(b.pos[2,0])
#     mat_list = [np.array(x_list), np.array(y_list), np.array(z_list)]
#     matrix = np.array(mat_list)
#     return matrix
    
#==============================================================================

np.random.seed(11111)

MAX_VEL = 25
MAX_ACC = 2
MAX_FORCE = 10
WALL_FORCE = MAX_FORCE

C_STR = 0.4
S_STR = 0.8
A_STR = 0.5

num_boids = 100
dt = 1

flock = gen_flock(num_boids)
positions = concVecs(flock)


plot = True
timing = False

###############################################################################
# TIMING
###############################################################################

def timeTenSecs(flock, dt):
    time = 0
    while time < 11:
        update(flock, dt)
        #print(time)
        time += dt
          
def timeTenSecs1(flock, dt):
    time = 0
    while time < 11:
        update1(flock, dt)
        time += dt  

     
if timing:
    setup = "from __main__ import timeTenSecs, flock, dt"
    t_tensec = timeit.Timer(stmt = "timeTenSecs(flock,dt)", setup=setup).repeat(1,1)
    print(st.mean(t_tensec))
    
    # setup = "from __main__ import timeTenSecs1, flock, dt"
    # t_tensec = timeit.Timer(stmt = "timeTenSecs1(flock,dt)", setup=setup).repeat(3,1)
    # print(st.mean(t_tensec))
    
    
    # setup = "from __main__ import bulkFindNeighb, flock, positions"
    # t_tensec = timeit.Timer(stmt = "bulkFindNeighb(positions, flock)", setup=setup).repeat(1,1)
    # print(st.mean(t_tensec))
    
    # setup = "from __main__ import bulkFindNeigbNoNumpy, flock, positions"
    # t_tensec = timeit.Timer(stmt = "bulkFindNeigbNoNumpy(positions, flock)", setup=setup).repeat(1,1)
    # print(st.mean(t_tensec))
    

###############################################################################
# ANIMATION and PLOTTING
###############################################################################

def animate(frame):
    update(flock, dt)
    positions = concVecs(flock)
    scatter._offsets3d = (positions[0, :], positions[1, :], positions[2, :])


if plot:
    positions = concVecs(flock)
    
    # Setup plot for animation
    plt.ion()
    fig = plt.figure()
    ax = axes3d.Axes3D(fig)
    ax.set_xlim3d(-200, 1200)
    ax.set_ylim3d(-200, 1200)
    ax.set_zlim3d(-200, 1200)
    
    # assign initial positions - ie create frame 1
    scatter = ax.scatter(positions[0, :], positions[1, :], positions[2, :],
                          marker='o', edgecolor='k', lw=0.5, s=2)
    
    # This must be performed in main body to work!
    anim = animation.FuncAnimation(fig, animate, frames=480, interval=50)
    
    # Writer = animation.writers['ffmpeg']
    # writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=5000)
    # anim.save("boids_3d_cyth.mp4", writer=writer)
