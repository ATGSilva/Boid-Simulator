import cython
cimport openmp

import numpy as np
cimport numpy as np

from libc.math cimport sqrt

from matplotlib import animation
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

DTYPE = np.double
ctypedef np.double_t DTYPE_t


cpdef findNear(object boid, list flock):
    cdef object b
    cdef list exclusive_flock, near_list, near_dists, sep_list, sep_dists
    cdef int near_alert, sep_alert
    cdef np.ndarray pos_diff, sqr_disp
    cdef double dist
    
    exclusive_flock = flock.copy()
    exclusive_flock.pop(flock.index(boid))
    near_alert = 120
    sep_alert = 5
    near_list = []
    near_dists = []
    sep_list = []
    sep_dists = []
    
    for b in exclusive_flock:
        pos_diff = b.pos - boid.pos
        sqr_disp = pos_diff * pos_diff
        dist = sqrt(np.sum(sqr_disp, axis=0))
        
        if dist < near_alert:
            near_list.append(b.id)
            near_dists.append(pos_diff)
        
        if dist < sep_alert:
            sep_list.append(b.id)
            sep_dists.append(pos_diff)
            
    boid.near_list = near_list
    boid.near_dists = near_dists
    boid.sep_list = sep_list
    boid.sep_dists = sep_dists
    
    
cpdef cohereForce(object boid, list flock, float MAX_VEL, float C_STR):
    cdef list near_indicies
    cdef np.ndarray pos_sum, target, steer
    cdef double target_mag
    cdef float m, m_sum
    cdef int indx
    
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
    boid.steer += steer


cpdef sepForce(object boid, list flock, float MAX_VEL, float S_STR):
    cdef list near_indicies
    cdef np.ndarray pos_sum, target, steer
    cdef double target_mag
    cdef float m, m_sum
    cdef int indx
    
    sep_indicies = boid.sep_list
    pos_sum = np.zeros((3,1))
    m_sum = 0
    
    if len(sep_indicies) == 0:
        return
    
    for indx in sep_indicies:
        sep = np.subtract(boid.pos, flock[indx].pos)
        m = flock[indx].mass
        m_sum += m
        pos_sum = np.add(pos_sum, m*sep)
    
    target = pos_sum * (1/m_sum)
    target_mag = np.linalg.norm(target)
    target = target * (MAX_VEL/target_mag) if target_mag != 0 else target
    steer = np.subtract(target, boid.vel) * S_STR
    boid.steer += steer
    
    
cpdef alignForce(object boid, list flock, float MAX_VEL, float A_STR):
    cdef list near_indicies
    cdef np.ndarray vel_sum, target, steer
    cdef double target_mag
    cdef float m, m_sum
    cdef int indx
    
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
    steer = np.subtract(vel_sum, boid.vel) * A_STR
    boid.steer += steer


cpdef reset(object boid):
    boid.steer = np.zeros((3,1))
    boid.near_list = None
    boid.near_dists = None
    boid.sep_list = None
    boid.sep_dists = None


cpdef updatePosEuler(object boid, list flock, float dt, 
                   float MAX_FORCE, float MAX_ACC, float MAX_VEL):
    cdef np.ndarray f_val, f_wall, f_sum, acc, vel, pos
    cdef float m
    
    f_val = boid.steer
    f_wall = wallForce(boid, MAX_FORCE)
    f_sum = limitVec(f_val, MAX_FORCE)
    f_sum = np.add(f_val, f_wall)
    
    m = boid.mass
    
    acc = (1/m) * f_sum
    boid.acc = limitVec(acc, MAX_ACC)
    vel = np.add(boid.vel, boid.acc * dt)
    boid.vel = limitVec(vel, MAX_VEL)
    boid.pos = np.add(boid.pos, boid.vel * dt)
        

#==============================================================================        


cpdef wallForce(object boid, float MAX_FORCE):
    cdef float WALL_FORCE
    cdef np.ndarray f_wall, out_bound_pos, out_bound_neg
    
    WALL_FORCE = MAX_FORCE
    f_wall = np.zeros((3,1))
    
    if np.any(boid.pos > 1000):
        out_bound_pos = -(boid.pos > 800).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_pos * np.abs(boid.pos))
    if np.any(boid.pos < 0):
        out_bound_neg = (boid.pos < 0).astype(int)
        f_wall = np.add(f_wall, WALL_FORCE * out_bound_neg * np.abs(boid.pos))
        
    return f_wall


cpdef limitVec(np.ndarray vec, float limit):
    cdef double vec_mag

    vec_mag = np.linalg.norm(vec)
    if abs(vec_mag) > limit:
        return (limit / vec_mag) * vec
    else:
        return vec


cpdef update(list flock, float dt, float m_vel, float m_acc, float m_for,
           float c_str, float s_str, float a_str):
    cdef float MAX_VEL, MAX_ACC, MAX_FORCE, WALL_FORCE, C_STR, S_STR, A_STR
    cdef object b
    
    MAX_VEL = m_vel
    MAX_ACC = m_acc
    MAX_FORCE = m_for
    WALL_FORCE = MAX_FORCE
    
    C_STR = c_str
    S_STR = s_str
    A_STR = a_str
    
    for b in flock:
        findNear(b, flock)
        cohereForce(b, flock, MAX_VEL, C_STR)
        sepForce(b, flock, MAX_VEL, S_STR)
        alignForce(b, flock, MAX_VEL, A_STR)
        updatePosEuler(b, flock, dt, MAX_FORCE, MAX_ACC, MAX_VEL)

        
cpdef concVecs(list flock):
    cdef object b
    cdef np.ndarray ini_vec
    
    for b in flock:
        if flock.index(b) == 0:
            ini_vec = flock[0].pos
        else:
            ini_vec = np.concatenate((ini_vec, b.pos), axis=1)
    return ini_vec

    

