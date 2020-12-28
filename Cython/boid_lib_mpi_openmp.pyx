import cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt, pow
from mpi4py import MPI
from cython.parallel cimport prange
import time
import sys


#--------------------------------------------------------------------------
# OPTIMISATIONS
#--------------------------------------------------------------------------


DTYPE = np.double
ctypedef np.double_t DTYPE_t

cdef:
    float MAX_VEL = 25
    float MAX_ACC = 2
    float MAX_FORCE = 10
    float WALL_FORCE = MAX_FORCE

    float C_STR = 0.4
    float S_STR = 0.9
    float A_STR = 0.8

#--------------------------------------------------------------------------
# NEIGHBOUR FINDING
#--------------------------------------------------------------------------


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef bulkFindNeighb(np.ndarray positions, list flock):
    cdef:
        int near_alert = 275
        int close_alert = 10
        double [:,::1] pos = positions
        Py_ssize_t nrow = pos.shape[0]
        Py_ssize_t ncol = pos.shape[1]
        Py_ssize_t i, j
        double [:,::1] dists = np.zeros((ncol, ncol), np.double)
        long [:,::1] near_boids = np.zeros((ncol, ncol), np.int)
        long [:,::1] close_boids = np.zeros((ncol, ncol), np.int)
        list bulk_nl = []
        list bulk_cl = []
        double dist_ij
        int close_tag, near_tag
        
    # Find matrix of inter-boid distances
    for i in prange(ncol-1, nogil=True, schedule='guided'):
    # for i in range(ncol-1):
        for j in range(i+1, ncol):
            dists[i,j] = sqrt(pow((pos[0,i] - pos[0,j]), 2) 
                              + pow((pos[1,i] - pos[1,j]), 2)
                              + pow((pos[2,i] - pos[2,j]), 2))
            dists[j,i] = dists[i,j]

    # Find indicies of boids that are near or close
    for i in prange(ncol-1, nogil=True, schedule='guided'):
    # for i in range(ncol-1):
        for j in range(i+1, ncol):
            dist_ij = dists[i,j]
            if dist_ij < close_alert:
                with gil:
                    bulk_cl.append([i,j])
                    bulk_cl.append([j,i])
            elif dist_ij < near_alert:
                with gil:
                    bulk_nl.append([i, j])
                    bulk_nl.append([j, i])
    
    t3 = time.time()
    # Append properties of near boids to boid object
    for [i, j] in bulk_nl:
        flock[i].near_list.append(j)
        flock[i].near_masses.append(flock[j].mass)
        flock[i].near_pos.append(flock[j].pos)
        flock[i].near_vels.append(flock[j].vel)

    # Append properties of close boids to boid object
    for [i, j] in bulk_cl:
        flock[i].close_list.append(j)
        flock[i].close_masses.append(flock[j].mass)
        flock[i].close_pos.append(flock[j].pos)


#--------------------------------------------------------------------------
# FORCES
#--------------------------------------------------------------------------

    
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef void cohereForce(object boid):
    cdef Py_ssize_t num_near = len(boid.near_list)
     
    if num_near == 0:
        return
    
    cdef:
        float [::1] near_masses = np.array(boid.near_masses, np.float32)
        list near_pos = boid.near_pos
        list near_indicies
        np.ndarray pos_sum, target, steer
        double target_mag
        float m, m_sum
        int indx
    
    pos_sum = np.zeros((3,1))
    for indx in range(num_near):
        pos_sum = np.add(pos_sum, (1/near_masses[indx])*near_pos[indx])
    
    m_sum = sum(near_masses)
    target = np.subtract((1/m_sum)*pos_sum, boid.pos)
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    boid.f_cohere = np.subtract(target, boid.vel) * C_STR


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef void sepForce(object boid):
    cdef:
        list near_indicies
        np.ndarray pos_sum, target, steer, sep 
        double target_mag
        float m, m_sum
        int indx
        Py_ssize_t num_close = len(boid.close_list)
    
    if num_close == 0:
        return
    
    cdef float [::1] close_masses = np.array(boid.close_masses, np.float32)
    cdef list close_pos = boid.close_pos
    
    pos_sum = np.zeros((3,1))
    for indx in range(num_close):
        sep = np.subtract(boid.pos, close_pos[indx])
        pos_sum = np.add(pos_sum, close_masses[indx]*sep)
    
    m_sum = sum(close_masses)
    target = pos_sum * (1/m_sum)
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    boid.f_sep = np.subtract(target, boid.vel) * S_STR
    

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef void alignForce(object boid):
    cdef Py_ssize_t num_near = len(boid.near_pos)
    
    if num_near == 0:
        return
    
    cdef:
        list near_indicies
        np.ndarray vel_sum, target, steer
        double target_mag
        float m, m_sum
        int indx
        float [::1] near_masses = np.array(boid.near_masses, np.float32)
        list near_vels = boid.near_vels
    
    vel_sum = np.zeros((3,1))
    for indx in range(num_near):
        vel_sum = np.add(vel_sum, near_masses[indx]*near_vels[indx])
    
    m_sum = sum(near_masses)
    target = (1/m_sum) * vel_sum
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    boid.f_align = np.subtract(target, boid.vel) * A_STR


cdef void wallForce(object boid):
    cdef np.ndarray f_wall, out_bound_pos, out_bound_neg
    
    f_wall = np.zeros((3,1))
    
    if np.any(boid.pos > 1500):
        out_bound_pos = -(boid.pos > 1500).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_pos * np.abs(boid.pos))
    if np.any(boid.pos < 0):
        out_bound_neg = (boid.pos < 0).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_neg * np.abs(boid.pos))
        
    boid.f_wall = f_wall


#--------------------------------------------------------------------------
# POS UPDATE AND RESET BOID PROPERTIES
#--------------------------------------------------------------------------


cdef void reset(object boid):
    boid.f_cohere = np.zeros((3,1), np.double)
    boid.f_sep = np.zeros((3,1), np.double)
    boid.f_align = np.zeros((3,1), np.double)
    boid.f_wall = np.zeros((3,1), np.double)
    boid.near_list = []
    boid.near_pos = []
    boid.near_vels = []
    boid.near_masses = []
    boid.close_list = []
    boid.close_pos = []
    boid.close_masses = []


cdef void updatePosEuler(object boid, float dt):
    cdef: 
        np.ndarray f_val, f_wall, f_sum, acc, vel, pos
        float m
    
    f_val = boid.f_cohere + boid.f_sep + boid.f_align
    f_wall = boid.f_wall
    f_sum = limitVec(f_val, MAX_FORCE)
    f_sum = np.add(f_val, f_wall)
    m = boid.mass
    
    acc = (1/m) * f_sum
    boid.acc = limitVec(acc, MAX_ACC)
    vel = np.add(boid.vel, boid.acc * dt)
    boid.vel = limitVec(vel, MAX_VEL)
    boid.pos = np.add(boid.pos, boid.vel * dt)
              

#--------------------------------------------------------------------------
# UPDATE - MAIN
#--------------------------------------------------------------------------


cpdef update(list boids, float dt):
    cdef: 
        object b
        double t1, t2, t3, t4, t5, t6, t7
    
    for b in boids:
        # if boids.index(b) == 0:
        #     t1 = time.time()
        #     cohereForce(b)
        #     t2 = time.time()
        #     sepForce(b)
        #     t3 = time.time()
        #     alignForce(b)
        #     t4 = time.time()
        #     wallForce(b)
        #     t5 = time.time()
        #     updatePosEuler(b, dt)
        #     t6 = time.time()
        #     reset(b)
        #     t7 = time.time()
        #     print(t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6)
        # else:
        cohereForce(b)
        sepForce(b)
        alignForce(b)
        wallForce(b)
        updatePosEuler(b, dt)
        reset(b)


#--------------------------------------------------------------------------
# UTILITIES
#--------------------------------------------------------------------------


cdef np.ndarray limitVec(np.ndarray vec, float limit):
    cdef double vec_mag

    vec_mag = np.linalg.norm(vec)
    if abs(vec_mag) > limit:
        return (limit / vec_mag) * vec
    else:
        return vec


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef concVecs(list flock):
    cdef: 
        list x_list = []
        list y_list = []
        list z_list = []
        list mat_list
        np.ndarray matrix
        object b
        
    for b in flock:
        x_list.append(b.pos[0,0])
        y_list.append(b.pos[1,0])
        z_list.append(b.pos[2,0])
    mat_list = [np.array(x_list), np.array(y_list), np.array(z_list)]
    matrix = np.array(mat_list)
    return matrix






