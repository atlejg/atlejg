from pylab import *
import AtlejgTools.SimulationTools.UnitConversion as u

def _exp10(x): return exp(log(10.)*x)

def walther_kinematic_viscosity(T1_ref, T2_ref, nu1_ref, nu2_ref, T):
    '''
    calculates kinematic [cSt] viscosity based on two reference points.
    T1_ref  : first point temp [K]
    T2_ref  : second point temp [K]
    nu1_ref : first point viscosity [cSt]
    nu2_ref : second point viscosity [cSt]
    ref: Honeywell report "Evaluation of Crude Corrsivity of Crude Oils" jan 2010, Equation 1
    see mail from Jorun Albertsen 3/5-2010
    '''
    # gonna solve Mx = y w.r.t. x=[A, B]
    M = matrix([ [1., -log10(T1_ref)],
                 [1., -log10(T2_ref)] ])
    y = matrix([log10(log10(nu1_ref + 0.7)), log10(log10(nu2_ref + 0.7))])
    x = linalg.inv(M) * y.transpose()
    # matrixes more complex than in matlab...
    A = x[0,0]
    B = x[1,0]
    #stop
    # use walthers equation
    return _exp10(_exp10(A - B*log10(T))) - 0.7

def density(T_ref, rho_ref, T):
    '''
    calculates density based on one reference point.
    T1_ref  : reference temp [K or C]
    rho_ref : reference density [any]
    T       : temp [K or C]
    ref: Honeywell report "Evaluation of Crude Corrsivity of Crude Oils" jan 2010, Equation 2
    see mail from Jorun Albertsen 3/5-2010
    the constant 'a' was determined from Table 4. in Honeywell report
    '''
    a = 0.00074
    return rho_ref / (1. + a*(T-T_ref))
