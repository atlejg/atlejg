import AtlejgTools.SimulationTools.UnitConversion as u

'''
assumes all values are si
   r_i   : inner radius [m]
   r_o   : outer radius [m]
   omega : rotational speed [rad/s]
   mu    : dynamic viscosity [kg/m/s]
   rho   : density
assumes outer wall is stationary

References:
   r1: Technical Note 2006-1: Pine Research Instrumentation. http://www.pineinst.com/echem/files/LMECN200601.pdf
   r2: http://astro.berkeley.edu/~jrg/ay202/node140.html (not online. have a printout)
   r3: http://www.pineinst.com/echem/files/LMECN200601.pdf
   r4: D.C. Silverman. Corrosion-vol61, no6, p515ff
'''

def reynolds_number(r_i, omega, rho, mu):
    return rho * omega*r_i * (2*r_i) / mu # r1 eq.1

def taylor_number(r_i, r_o, omega, rho, mu):
    d = r_o - r_i
    d2 = r_o**2 - r_i**2
    Ta = 4 * omega**2 * r_i**2 * d**4 / d2 / (mu/rho)**2
    return Ta

def wall_shear_stress(r_i, r_o, omega, rho, mu):
    ta = taylor_number(r_i, r_o, omega, rho, mu)
    re = reynolds_number(r_i, omega, rho, mu)
    print('ta = %g re = %g' % (ta, re))
    if  ta < critical_taylor_number(omega, 0):
        # now we have couette flow
        d = r_o - r_i
        d2 = r_o**2 - r_i**2
        A = -omega * r_i**2 / d2
        B = omega * r_i**2 * r_o**2 / d2
        #
        dudr_o = A - B/r_o**2 # du/dr @ r=r_o ( == -2*A)
        dudr_i = A - B/r_i*2 # du/dr @ r=r_i
        tau_o = -mu * dudr_o
        tau_i = -mu * dudr_i
        return tau_i, tau_o
    elif re > 200:
        # using r3 eq 5 (turbulent)
        return 0.0791 * re**-0.3 * rho * r_i**2 * omega**2
    else:
        raise Exception("dont know what to do")

def critical_taylor_number(omega_inner, omega_outer):
    return 2 * 1707.762 / (1. + omega_outer/omega_inner)

def bob_radius_ala_silverman(d_pipe, v_pipe, rho, mu, Sc):
    # using r4 eq 5
    d_bob = (8.442 * d_pipe**0.1786 * Sc**0.0857 * (mu/rho)**0.25 * v_pipe**-0.25)**2.333
    return d_bob / 2.

def omega(shear_stress, rho, mu, Sc=None, d_pipe=None, v_pipe=None, r_i=None):
    if r_i == None: r_i = bob_radius_ala_silverman(d_pipe, v_pipe, rho, mu, Sc)
    # new find omega s.t. shear stress on bob is equal to given shear stress
    omega = (shear_stress * (2*r_i**2/(mu/rho))**0.3 / (0.0791*r_i**2*rho))**(1./1.7)
    #stop
    return omega, r_i
