'''
all values are si
'''

import pylab as pl

RE_LAMINAR = 2300

def mean_velocity(q, diam):
    '''
    all values are si
    '''
    A = pl.pi * diam**2 / 4.
    return q / A

def reynolds_number(q, diam, rho, mu):
    '''
    all values are si
    '''
    return rho * mean_velocity(q, diam) * diam / mu

def friction_factor(q=0, diam=0, rho=0, mu=0, re=None):
    '''
    all values are si
    ref: White p.345 eq.6.55 (turb) and p.342 eq.6.46 (lamin)
    '''
    if re == None: re = reynolds_number(q, diam, rho, mu)
    #print "PipeFlow: re =", re
    if   re < RE_LAMINAR: return 64. / re
    elif re < 1e5       : return 0.316 / re**0.25
    else                : return (1.803*pl.log(re/6.9)/pl.log(10.))**-2. # the constant 1.8 adjusted to have better continuity

def pressure_drop(q, diam, rho, mu, length):
    '''
    all values are si
    '''
    return 0.5 * rho * mean_velocity(q, diam)**2 * friction_factor(q, diam, rho, mu) * length/diam

def wall_shear_stress(q, diam, rho, mu):
    '''
    all values are si
    ref: White p.340 eq.6.29
    '''
    return rho * mean_velocity(q, diam)**2 * friction_factor(q, diam, rho, mu) / 8.

if __name__ == '__main__' :
    from pylab import *
    import AtlejgTools.SimulationTools.UnitConversion as u

    d = 4*u.INCH
    v = 3.
    rho = 900.
    mu  = 2*u.CP

    q = v * pi*d**2/4.

    print(friction_factor(q, d, rho, mu))
    print(wall_shear_stress(q, d, rho, mu))
