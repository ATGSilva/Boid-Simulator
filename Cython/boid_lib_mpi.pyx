import cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt, pow
from mpi4py import MPI
from cython.parallel cimport prange
import time
import sys


#--------------------------------------------------------------------------
# OPTIMISATIONS & SETTINGS
#--------------------------------------------------------------------------


ctypedef np.double_t DTYPE_t

cdef:
    float MAX_VEL = 25
    float MAX_ACC = 2
    float MAX_FORCE = 10
    float WALL_FORCE = MAX_FORCE

    float C_STR = 0.9
    float S_STR = 0.6
    float A_STR = 0.7

#--------------------------------------------------------------------------
# NEIGHBOUR FINDING
#--------------------------------------------------------------------------

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef bulkFindNeighb(np.ndarray positions, list flock):
    cdef:
        int near_alert = 300
        int close_alert = 10
        double [:,::1] pos = positions
        Py_ssize_t ncol = len(flock)
        Py_ssize_t i, j
        double [:,::1] dists = np.zeros((ncol, ncol), np.double)
        list bulk_nl = []
        list bulk_cl = []
        double dist_ij
        double t1, t2, t3, t4
    
    #t1 = time.time()
    for i in prange(ncol-1, nogil=True, schedule='guided'):
        for j in range(i+1, ncol):
            dists[i,j] = sqrt(pow((pos[0,i] - pos[0,j]), 2) 
                              + pow((pos[1,i] - pos[1,j]), 2)
                              + pow((pos[2,i] - pos[2,j]), 2))
            dists[j,i] = dists[i,j]
    #t2 = time.time()
    
    for i in prange(ncol-1, nogil=True, schedule='guided'):
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
    #t3 = time.time()
    
    for [i, j] in bulk_nl:
        flock[i].near_list.append(j)
        flock[i].near_masses.append(flock[j].mass)
        flock[i].near_pos.append(flock[j].pos)
        flock[i].near_vels.append(flock[j].vel)
    
    for [i, j] in bulk_cl:
        flock[i].close_list.append(j)
        flock[i].close_masses.append(flock[j].mass)
        flock[i].close_pos.append(flock[j].pos)
    #t4 = time.time()
    #print(t2-t1, t3-t2, t4-t3)

#--------------------------------------------------------------------------
# FORCES
#--------------------------------------------------------------------------


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
cdef void cohereForce(object boid):
    cdef Py_ssize_t num_near = len(boid.near_list)
    
    if num_near == 0:
        return
                                   
    cdef:
        double [:,::1] pos_sum = np.zeros((3,1), np.double)
        double [:,::1] target = np.zeros((3,1), np.double)
        double [:,::1] boid_pos = np.array(boid.pos)
        double [:,::1] boid_vel = np.array(boid.vel)
        float [::1] near_masses = np.array(boid.near_masses, np.float32)
        double [:,::1] f_cohere = np.zeros((3,1), np.double)
        float m_sum
        int indx, j
        double [:,:] near_pos = np.hstack(np.array(boid.near_pos))
    
    m_sum = sum(near_masses)
    for j in prange(3, nogil=True, schedule='dynamic'):
        for indx in range(num_near):
            pos_sum[j,0] = pos_sum[j,0] + (1/near_masses[indx]) * near_pos[j,indx]
        target[j,0] = (1/m_sum)*pos_sum[j,0] - boid_pos[j,0]
    
    target = normTarg(target)
    
    for j in prange(3, nogil=True, schedule='static'):
        f_cohere[j,0] = (target[j,0] - boid_vel[j,0]) * C_STR
    
    f_cohere = limitVec(f_cohere, MAX_ACC)
    
    boid.f_cohere = np.array(f_cohere)


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
cdef void sepForce(object boid):
    cdef Py_ssize_t num_close = len(boid.close_list)
    
    if num_close == 0:
        return
    
    cdef:
        double [:,::1] boid_pos = np.array(boid.pos)
        double [:,::1] boid_vel = np.array(boid.vel)
        double [:,::1] pos_sum = np.zeros((3,1), np.double)
        double [:,::1] sep = np.zeros((3,1), np.double) 
        double [:,::1] target = np.zeros((3,1))
        double [:,::1] f_sep = np.zeros((3,1), np.double)
        float m_sum
        int indx, j
        float [::1] close_masses = np.array(boid.close_masses, np.float32)
        double [:,:] close_pos = np.hstack(np.array(boid.close_pos))
    
    m_sum = sum(close_masses)
    for j in prange(3, nogil=True, schedule='dynamic'):
        for indx in range(num_close):
            sep[j,0] = boid_pos[j,0] - close_pos[indx,j]
            pos_sum[j,0] = pos_sum[j,0] + close_masses[indx]*sep[j,0]
        target[j,0] = (1/m_sum) * pos_sum[j,0]
    
    target = normTarg(target)
    
    for j in prange(3, nogil=True, schedule='static'):
        f_sep[j,0] = (target[j,0] - boid_vel[j,0]) * S_STR
        
    f_sep = limitVec(f_sep, MAX_ACC)
    
    boid.f_sep = np.array(f_sep)
    

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef void alignForce(object boid):
    cdef Py_ssize_t num_near = len(boid.near_pos)
    
    if num_near == 0:
        return
    
    cdef:
        double [:,::1] boid_vel = np.array(boid.vel)
        double [:,::1] vel_sum = np.zeros((3,1), np.double)
        double [:,::1] target = np.zeros((3,1), np.double)
        float [::1] near_masses = np.array(boid.near_masses, np.float32)
        double [:,::1] f_align = np.zeros((3,1), np.double)
        double target_mag
        float m_sum
        int indx, j
        double [:,:] near_vels = np.hstack(np.array(boid.near_vels))
    
    m_sum = sum(near_masses)
    for j in prange(3, nogil=True, schedule='dynamic'):
        for indx in range(num_near):
            vel_sum[j,0] = vel_sum[j,0] + near_masses[indx] * near_vels[j,indx]
        target[j,0] = (1/m_sum) * vel_sum[j,0]
    
    target = normTarg(target)

    for j in prange(3, nogil=True, schedule='static'):
        f_align[j,0] = (target[j,0] - boid_vel[j,0]) * A_STR
    
    f_align = limitVec(f_align, MAX_ACC)
    
    boid.f_align = np.array(f_align)


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef void wallForce(object boid):
    cdef np.ndarray f_wall, out_bound_pos, out_bound_neg
    
    f_wall = np.zeros((3,1))
    boid_pos = np.array(boid.pos)
    
    if np.any(boid_pos > 800):
        out_bound_pos = -(boid_pos > 800).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_pos * np.abs(boid_pos))
    if np.any(boid_pos < 0):
        out_bound_neg = (boid_pos < 0).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_neg * np.abs(boid_pos))
        
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


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef void updatePosEuler(object boid, float dt):
    cdef:
        double [:,::1] f_val = np.zeros((3,1), np.double)
        double [:,::1] f_wall = boid.f_wall
        double [:,::1] f_sum = np.zeros((3,1), np.double)
        double [:,::1] f_cohere = boid.f_cohere
        double [:,::1] f_sep = boid.f_sep
        double [:,::1] f_align = boid.f_align
        double [:,::1] boid_acc = np.array(boid.acc, np.double)
        double [:,::1] boid_vel = np.array(boid.vel, np.double)
        double [:,::1] boid_pos = np.array(boid.pos, np.double)
        Py_ssize_t i
        float m = boid.mass
    
    for i in range(3):
        f_val[i][0] = f_cohere[i][0] + f_sep[i][0] + f_align[i][0]
    
    f_sum = limitVec(f_val, MAX_FORCE)
    
    for i in range(3):
        f_sum[i][0] = f_val[i][0] + f_wall[i][0]
        boid_acc[i][0] = (1/m) * f_sum[i][0]
        
    boid_acc = limitVec(boid_acc, MAX_ACC)
    
    for i in range(3):
        boid_vel[i][0] = boid_vel[i][0] + (dt * boid_acc[i][0])

    boid_vel = limitVec(boid_vel, MAX_VEL)
    
    for i in range(3):
        boid_pos[i][0] = boid_pos[i][0] + (dt * boid_vel[i][0])
    
    boid.acc = list(np.array(boid_acc))
    boid.vel = list(np.array(boid_vel))
    boid.pos = list(np.array(boid_pos))
              

# cdef void updatePosEuler(object boid, float dt):
#     cdef np.ndarray f_val, f_wall, f_sum, acc, vel, pos
#     cdef float m
#     cdef:
#         boid_acc = np.array(boid.acc)
#         boid_vel = np.array(boid.vel)
#         boid_pos = np.array(boid.pos)
    
#     f_val = boid.f_cohere + boid.f_sep + boid.f_align
#     f_wall = boid.f_wall
#     f_sum = limitNumpy(f_val, MAX_FORCE)
#     f_sum = np.add(f_val, f_wall)
    
#     m = boid.mass
    
#     acc = (1/m) * f_sum
#     boid_acc = limitNumpy(acc, MAX_ACC)
#     boid_vel = np.add(boid_vel, boid_acc * dt)
#     boid_vel = limitNumpy(boid_vel, MAX_VEL)
#     boid_pos = np.add(boid_pos, boid_vel * dt)

#     boid.acc = list(boid_acc)
#     boid.vel = list(boid_vel)
#     boid.pos = list(boid_pos)

#--------------------------------------------------------------------------
# UPDATE - MAIN
#--------------------------------------------------------------------------


cpdef update(list boids, float dt):
    cdef object b
    #cdef double t1, t2, t3, t4, t5, t6, t7
    
    for b in boids:
        if boids.index(b) == 0:
            t1 = time.time()
            cohereForce(b)
            t2 = time.time()
            sepForce(b)
            t3 = time.time()
            alignForce(b)
            t4 = time.time()
            wallForce(b)
            t5 = time.time()
            updatePosEuler(b, dt)
            t6 = time.time()
            reset(b)
            t7 = time.time()
            print(t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6)
        else:
            cohereForce(b)
            sepForce(b)
            alignForce(b)
            wallForce(b)
            updatePosEuler(b, dt)
            reset(b)

#--------------------------------------------------------------------------
# UTILITIES
#--------------------------------------------------------------------------


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
cdef double [:,::1] normTarg(double [:,::1] target):
    cdef:
        Py_ssize_t j
        double target_mag, norm_fact
        double [:,::1] norm_dims = np.zeros((3,1), np.double)
    
    for j in prange(3, nogil=True, schedule='static'):
        norm_dims[j,0] = pow(target[j,0], 2)
    target_mag = sqrt(norm_dims[0,0] + norm_dims[1,0] + norm_dims[2,0])
    
    if target_mag != 0:
        norm_fact = MAX_VEL/target_mag
        for j in prange(3, nogil=True):
            target[j,0] = norm_fact * target[j,0]
    
    return target


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
cdef double [:,::1] limitVec(double [:,::1] vec, float limit):
    cdef double vec_mag
    cdef Py_ssize_t i, j
    cdef double [:,::1] norm_dims = np.zeros((3,1), np.double)
    cdef double norm, norm_fact

    for i in range(3):
        norm_dims[i][0] = pow(vec[i][0], 2)
    norm = sqrt(norm_dims[0][0] + norm_dims[1][0] + norm_dims[2][0])
    
    if norm > limit:
        norm_fact = (limit / norm)
        for i in range(3):
            vec[i][0] = norm_fact * vec[i][0]
        return vec
    else:
        return vec


cdef np.ndarray limitNumpy(np.ndarray vec, float limit):
    cdef double vec_mag

    vec_mag = np.linalg.norm(vec)
    if abs(vec_mag) > limit:
        return (limit / vec_mag) * vec
    else:
        return vec


cpdef concVecs(list flock):
    cdef: 
        list x_list = []
        list y_list = []
        list z_list = []
        np.ndarray positions
        object b
        
    for b in flock:
        x_list.append(b.pos[0][0])
        y_list.append(b.pos[1][0])
        z_list.append(b.pos[2][0])
    positions = np.array([x_list, y_list, z_list])
    return positions






