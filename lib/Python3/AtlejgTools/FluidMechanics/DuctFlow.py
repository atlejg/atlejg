'''
all values are si
here we assume infinite parallel plates.
see White chap. 6.6 p. 359 (4.th ed.)
note1: height is distance between plates (= 2*h in White)
       so, we assume height <= widht
note2: only implemented for laminar flow
'''

import pylab as pl

RE_LAMINAR = 2300

def mean_velocity(q, height, width):
    '''
    all values are si
    we assume height <= widht
    '''
    A = width*height
    return q / A

def reynolds_number(q, height, width, rho, mu):
    '''
    all values are si
    this is based on hydraulic diam
    we assume height <= widht
    '''
    Dh = 2*height
    return rho * mean_velocity(q, height, width) * Dh / mu

def friction_factor(q, height, width, rho, mu, re=None):
    '''
    all values are si
    ref: White p.345 eq.6.55 (turb) and p.342 eq.6.46 (lamin)
    we assume height <= widht
    '''
    if re == None: re = reynolds_number(q, height, width, rho, mu)
    if re < RE_LAMINAR: return 96. / re
    else: raise Exception("cannot handle turbulent flow (yet...)")

def pressure_drop(q, height, width, rho, mu, length):
    '''
    all values are si --> returns Pa
    we assume height <= widht
    '''
    Dh = 2*height
    return 0.5 * rho * mean_velocity(q, height, width)**2 * friction_factor(q, height, width, rho, mu) * length/Dh

def wall_shear_stress(q, height, width, rho, mu):
    '''
    all values are si
    ref: White p.340 eq.6.29
    we assume height <= widht
    '''
    return rho * mean_velocity(q, height, width)**2 * friction_factor(q, height, width, rho, mu) / 8.
