'''

do some calculations on a open channel flow inside an inclined pipe (like injector A11, Peregrino)

for now, only pipes (TYPE_PIPE) are supported.

the angle theta is like this (the pointed angle between dotted lines):


           %%%    %%%
      %%%              %%%

  %%%                      %%%

 %%%                         %%%
               .
 %%%          . .            %%%
             .   .
 %%%        .     .         %%%
           .       .
    %%%   .         .    %%%
         .           .
          %%%     %%%
'''

from scipy.optimize import fsolve
import pylab as pl

TYPE_PIPE = 'pipe'

def _flowarea(theta, r):
    return r**2/2. * (theta-pl.sin(theta))

def afa(theta, r):
    '''
    this is A*f(A)
    '''
    A = _flowarea(theta, r)
    return A**(5/3.) / (r*theta)**(2/3.)

def _manning_vmean(rh, n, s, k=1.):
    '''
    https://en.wikipedia.org/wiki/Manning_formula
    rh: hydraulic radius
    n: manning-roughness
    s: slope
    k: conversion factor. 1 for Si
    return: mean velocity [m/s]
    '''
    return k*rh**(2/3.)*pl.sqrt(s) / n

def _manning_func(q, n, s, k=1.):
    '''
    https://en.wikipedia.org/wiki/Manning_formula, after some manipulation
    q: flowrate
    n: manning-roughness
    s: slope
    k: conversion factor. 1 for Si
    return: value
    '''
    return q*n / (k*pl.sqrt(s))

def _hydr_rad(theta, r):
    p = r*theta      # wetted perimeter
    return _flowarea(theta, r)/p

def flow_characteristics(gtype, geom, inclin, q, n=0.022, k=1., full_output=False):
    '''
    given channel type, channel goemetry, flow-rate and inclination, it will
    calculate mean flow velocity and some geometric characteristics for the flow.
    # input
      gtype  : geometry type, f.ex. TYPE_PIPE
      geom   : struct with geometry details
      for TYPE_PIPE:
          geom.diam
      inclin : inclination of channel [deg]. 0 means horizontal
      q      : flowrate [m3/s]
      n      : manning friction factor. default is corrugated steel
      k      : unit conversion. 1. for Si
      full_output: use this one to get some info on the numerical solution (converged or not?)
    # output
      v      : mean velocity
      for TYPE_PIPE:
          theta  : angle discribing fluid channel (see above)
          w      : width of fluid channel (upper)
          h      : height of fluid channel
      ans : if full_output, this will be the full output from fsolve
    '''
    #
    s = pl.sin(inclin*pl.pi/180.)            # slope
    #
    if gtype == TYPE_PIPE:
        r = geom.diam/2.
        func = lambda theta: _flowarea(theta, r)*_manning_vmean(_hydr_rad(theta, r), n, s, k) - q
        theta0 = pl.pi/4.              # initial guess (0 got stuck)
        ans = fsolve(func, theta0, full_output=full_output)
        theta = ans[0][0] if full_output else ans[0]
        v = _manning_vmean(_hydr_rad(theta, r), n, s, k)
        h = r - r*pl.cos(theta/2.)
        w = 2*r*pl.sin(theta/2.)
        if full_output:
            return v, theta, w, h, ans   # angle, width, height, output from fsolve
        else:
            return v, theta, w, h        # angle, width, height
