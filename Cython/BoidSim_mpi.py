from itertools import accumulate

import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np 
import pandas as pd
import time as timer

from mpi4py import MPI
#import sys

import boid_lib_mpi as bl

###############################################################################
# CLASSES AND FUCTIONS
###############################################################################

class Boid():
    
    def __init__(self, pos, vel, acc, ident, mass=1.0, b_type="boid"):
        self.pos = pos
        self.vel = vel
        self.acc = acc
        self.id = ident
        self.mass = mass
        
        self.near_list = []
        self.near_pos = []
        self.near_vels = []
        self.near_masses = []
        self.close_list = []
        self.close_pos = []
        self.close_masses = []
        
        self.f_cohere = np.zeros((3,1), np.double)
        self.f_sep = np.zeros((3,1), np.double)
        self.f_align = np.zeros((3,1), np.double)
        self.f_wall = np.zeros((3,1), np.double)
        
        self.type = b_type

def genFlock(num):
    l_lim_pos = np.array([0, 0, 0])
    u_lim_pos = np.array([500, 500, 500])
    l_lim_vel = np.array([-10, -10, -10])
    u_lim_vel = np.array([10, 10, 10])
    
    w_pos = u_lim_pos - l_lim_pos
    w_vel = u_lim_vel - l_lim_vel
    
    acc = list(np.zeros((3,1)))
    flock = []
    for i in range(num):
        pos = list(l_lim_pos[:, np.newaxis] + np.random.random((3,1)) 
                   * w_pos[:, np.newaxis])
        vel = list(l_lim_vel[:, np.newaxis] + np.random.random((3,1)) 
                   * w_vel[:, np.newaxis])
        flock.append(Boid(pos, vel, acc, i))
    return flock    

def extraFlock(num, flock):
    l_lim_pos = np.array([700, 700, 700])
    u_lim_pos = np.array([900, 900, 900])
    l_lim_vel = np.array([-10, -10, -10])
    u_lim_vel = np.array([10, 10, 10])
    w_pos = u_lim_pos - l_lim_pos
    w_vel = u_lim_vel - l_lim_vel
    acc = np.array([0, 0, 0]).tolist()
    
    for i in range(num):
        pos = l_lim_pos[:, np.newaxis] + np.random.random((3,1)) * w_pos[:, np.newaxis]
        vel = l_lim_vel[:, np.newaxis] + np.random.random((3,1)) * w_vel[:, np.newaxis]
        flock.append(Boid(pos, vel, acc, i))

def splitFlock(flock, numworkers):
    num_boids = len(flock)
    # Split flock into equal lists
    floor_div = num_boids // numworkers
    remain_boids = num_boids % numworkers
    split_by = [floor_div for i in range(numworkers)]
    split_by[-1] = floor_div + remain_boids
    
    split_flock = [flock[x-y:x] for x, y in zip(accumulate(split_by), split_by)]

    return split_flock

###############################################################################
# MAIN
###############################################################################

timestart = timer.time()
#np.random.seed(11111)

MAXWORKER = 5
MINWORKER = 1
DIRECTOR = 0

# Simulation settings
DURATION = 121

# Tags
BEGIN_TAGS = [i for i in range(DURATION)]
DONE_TAGS = [i for i in range(DURATION, 2*DURATION)]
DATA_TAGS = [i for i in range(2*DURATION, 3*DURATION)]
POS_BEG_TAGS = [i for i in range(3*DURATION, 4*DURATION)]
POS_DONE_TAGS = [i for i in range(4*DURATION, 5*DURATION)]

#==============================================================================

def main():

    num_boids = 500
    dt = 1
    time = 0
    its = 0
    
    # MPI Setup
    comm = MPI.COMM_WORLD
    taskid = comm.Get_rank()
    numtasks = comm.Get_size()
    numworkers = numtasks-1
    split_flock = None
    
    #--------------------------------------------------------------------------
    # DIRECTOR CODE
    #--------------------------------------------------------------------------

    while time < DURATION:
        
        if taskid == DIRECTOR:
            
            if time == 0:
                if (numworkers > MAXWORKER) or (numworkers < MINWORKER):
                    print("ERROR: the number of tasks must be between %d and %d." % (MINWORKER+1,MAXWORKER+1))
                    print("Quitting...")
                    comm.Abort()
                    
                # Generate flock
                flock = genFlock(num_boids)
                #extraFlock(num_boids, flock)
                
                columns = ["Time", "ID", "X Position", "Y Position", "Z Position"]
                df = pd.DataFrame(columns=columns)
                
            if split_flock is not None:
                flock = []
                for split in split_flock:
                    flock += split
            
            # Find position matrix
            positions = bl.concVecs(flock)
            
            # Find Neighbours
            bl.bulkFindNeighb(positions, flock)
            # Split the flock up into approx equal sub-lists 
            
            split_flock = splitFlock(flock, numworkers)
            
            # Send to workers
            for worker_id in range(1, numworkers+1):
                comm.send(split_flock[worker_id-1], dest=worker_id, tag=BEGIN_TAGS[its])   
     
              # Recieve from workers
            worker_data = [[] for i in range(numworkers)]
            for worker_id in range(1, numworkers+1):
                split_flock[worker_id-1] = comm.recv(source=worker_id, tag=DONE_TAGS[its])
                worker_data[worker_id-1] = comm.recv(source=worker_id, tag=DATA_TAGS[its])
            
            
            # Save position data for later plotting
            _df = pd.concat([pd.DataFrame(i, columns=columns) for i in worker_data],
                            ignore_index=True)
            df = pd.concat([df, _df], ignore_index=True)
            
            print(time)
            
            if time == DURATION-dt:
                df.to_csv("positions.csv", index=False)
                timeend = timer.time()
                print(timeend-timestart)
                print("complete")
                
            # End of Director code
        
        #--------------------------------------------------------------------------
        # WORKER CODE
        #--------------------------------------------------------------------------
        
        elif taskid != DIRECTOR:
            
            boids = comm.recv(source=DIRECTOR, tag=BEGIN_TAGS[its])
            t1 = MPI.Wtime()
            bl.update(boids, dt)
            t2 = MPI.Wtime()
            print(t2-t1)
            
            data = []
            for b in boids:
                data.append([time, b.id, b.pos[0][0], b.pos[1][0], b.pos[2][0]])
            
            comm.send(boids, dest=DIRECTOR, tag=DONE_TAGS[its])
            comm.send(data, dest=DIRECTOR, tag=DATA_TAGS[its])
            
            # End of Worker code
            
        time += dt
        its += 1
    
#==============================================================================
    
main()


