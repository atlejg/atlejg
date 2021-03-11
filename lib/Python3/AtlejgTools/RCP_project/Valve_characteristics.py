import AtlejgTools.Utils as UT
import AtlejgTools.SimulationTools.UnitConversion as U
import pylab as pl

_RHO_CAL = 1000.25
_MU_CAL  = 1.45

USE_STD_AICD_FUNC = True
RHO_EXP           = 1.

def rcp_coeff(rho, mu, rho_cal, mu_cal, cnst, y):
    '''
    for debug purposes, its useful to have this isolated
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    '''
    if USE_STD_AICD_FUNC: return (rho**2/rho_cal) * (mu_cal/mu)**y * cnst
    else                : return (rho**RHO_EXP/rho_cal) * (y-mu/mu_cal) * cnst

def rcp_dp(rho, mu, rho_cal, mu_cal, cnst, x, y, q):
    '''
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    q = q * 24. / 1000. # convert to m3/day
    return rcp_coeff(rho, mu, rho_cal, mu_cal, cnst, y) * q**x

def rcp_dp2(valve_type, rho, mu, q):
    '''
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    q = q * 24. / 1000. # convert to m3/day
    rho_cal, mu_cal, cnst, x, y = rcp_params(valve_type)
    return rcp_coeff(rho, mu, rho_cal, mu_cal, cnst, y) * q**x

def rcp_dp3(rcp_params, rho, mu, q):
    '''
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    rcp_params: (rho_cal, mu_cal, cnst, x, y)
    '''
    q = q * 24. / 1000. # convert to m3/day
    rho_cal, mu_cal, cnst, x, y = rcp_params
    return rcp_coeff(rho, mu, rho_cal, mu_cal, cnst, y) * q**x

def rcp_flowr(rho, mu, rho_cal, mu_cal, cnst, x, y, dp):
    '''
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    return 1000./24. * (dp/rcp_coeff(rho, mu, rho_cal, mu_cal, cnst, y))**(1./x)

def rcp_flowr2(valve_type, rho, mu, dp):
    '''
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    rho_cal, mu_cal, cnst, x, y = rcp_params(valve_type)
    return 1000./24. * (dp/rcp_coeff(rho, mu, rho_cal, mu_cal, cnst, y))**(1./x)

def rcp_flowr3(rcp_params, rho, mu, dp):
    '''
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    rcp_params: (rho_cal, mu_cal, cnst, x, y)
    '''
    rho_cal, mu_cal, cnst, x, y = rcp_params
    return 1000./24. * (dp/rcp_coeff(rho, mu, rho_cal, mu_cal, cnst, y))**(1./x)

def icd_coeff(rho, mu, cnst, rho_cal=_RHO_CAL, mu_cal=_MU_CAL):
    '''
    for debug purposes, its useful to have this isolated
    [rho] = kg/m3
    [mu]  = cp
    '''
    return (rho_cal/rho * mu/mu_cal)**.25 * rho/rho_cal * cnst

def sicd_coeff(rho, mu, bar_strength, rho_cal=_RHO_CAL, mu_cal=_MU_CAL):
    '''
    [rho] = kg/m3
    [mu]  = cp
    '''
    return (rho_cal/rho * mu/mu_cal)**.25 * rho/rho_cal * bar_strength

def icd_dp(rho, mu, cnst, q, rho_cal=_RHO_CAL, mu_cal=_MU_CAL):
    '''
    calibration vals are default eclipse
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    q = q * 24. / 1000. # convert to m3/day
    return icd_coeff(rho, mu, cnst, rho_cal, mu_cal) * q**2

def icd_flowr(rho, mu, cnst, dp, rho_cal=_RHO_CAL, mu_cal=_MU_CAL):
    '''
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    return 1000./24. * (dp/icd_coeff(rho, mu, cnst, rho_cal, mu_cal))**0.5

def icd_dp2(gor, rho_o, rho_g, mu_o, mu_g, cnst, q, rho_cal=_RHO_CAL, mu_cal=_MU_CAL):
    '''
    oil/gas only
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    so = 1. / (1. + gor) # oil saturation
    mu  = so*mu_o  + (1-so)*mu_g
    rho = so*rho_o + (1-so)*rho_g
    return icd_dp(rho, mu, cnst, q, rho_cal, mu_cal)




def icd_strength(bar_strength):
    '''
    values found from troll weldata-file. confirmed by edel's best practice document.
    for flexibility, we interpolate in between'''
    bar_strength = float(bar_strength) # in case it is a string
    x = []; y = []
    x.append(0.0); y.append(0.0)
    x.append(0.2); y.append(0.00033)
    x.append(0.4); y.append(0.00048)
    x.append(0.8); y.append(0.00081)
    x.append(1.6); y.append(0.00159)
    x.append(3.2); y.append(0.00346)
    return UT.interpolate(x, y, bar_strength)

def sicd_dp(bar_strength, rho, mu, q):
    '''
    bar_strength is typically 0.8, 1.6 or 3.2
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    return icd_dp(rho, mu, icd_strength(bar_strength), q)

def sicd_flowr(bar_strength, rho, mu, dp):
    '''
    bar_strength is typically 0.8, 1.6 or 3.2
    [rho] = kg/m3
    [mu]  = cp
    [q]   = l/h
    [dp]  = bar
    '''
    return icd_flowr(rho, mu, icd_strength(bar_strength), dp)

def mu_mix(mu_o, mu_g, gvf):
    '''
    [mu]  = cp
    '''
    return (1-gvf)*mu_o + gvf*mu_g

def rho_mix(rho_o, rho_g, gvf, a=1., c=1.):
    '''
    [rho] = kg/m3
    '''
    return (1-gvf)**a*rho_o + gvf**c*rho_g

def fluent_sicd_coefficients(area, length, rho, mu, bar_strength):
    '''
    calculates the coefficients for the power law of the porous zone in fluent that represents the valve.
    need to do this for each phase in the simulation.
    area        : cross-sectional area of the porous zone [m2]
    length      : streamwise length of  the porous zone [m]
    bar_strength: typically 0.4, 0.8, 3.2 etc
    [rho] = kg/m3
    [mu]  = cp
    '''
    c1 = 2.
    c0 = U.BAR * (_RHO_CAL/rho * mu/_MU_CAL)**.25 * (rho/_RHO_CAL) * icd_strength(bar_strength) * area**c1 * U.DAY**c1 / length
    if c0 > 1e20:
        print("warning: c0 = %e > 1e20 (which is max limit for fluent gui-input)" % c0)
    return (c0, c1)


def fluent_rcp_coefficients(area, length, rho, mu, valve_type=None, rho_cal=0., mu_cal=0., x=0., y=0., a_aicd=0.):
    '''
    calculates the coefficients for the power law of the porous zone in fluent that represents the valve.
    need to do this for each phase in the simulation.
    area      : cross-sectional area of the porous zone [m2]
    length    : streamwise length of  the porous zone [m]
    [rho] = kg/m3
    [mu]  = cp
    valve_type: see rcp_params above
    rho_cal   : calibration density       (needed if valve_type is not given)
    mu_cal    : calibration viscosity     (needed if valve_type is not given)
    x         : flowrate exponent         (needed if valve_type is not given)
    y         : phase separation exponent (needed if valve_type is not given)
    a_aicd    : the linear constant       (needed if valve_type is not given)
    '''
    if not valve_type is None:
        rho_cal, mu_cal, a_aicd, x, y = rcp_params(valve_type)
    c1 = x
    c0 = U.BAR * (rho**2 / rho_cal) * (mu_cal/mu)**y * a_aicd * area**c1 * U.DAY**c1 / length
    if c0 > 1e20:
        print("warning: c0 = %e > 1e20 (which is max limit for fluent gui-input)" % c0)
    return (c0, c1)

def sicd_length(dp, q, segm_length, rho, mu, bar_strength):
    '''
    calculates the distance between sicd's to obtain the given flow & dp on a segment.
    useful for optimizing sicd's.
    [dp] = bar
    [q]  = m3/d
    '''
    c = icd_coeff(rho, mu, icd_strength(bar_strength))
    l = (dp/c)**(1/2.) * segm_length / q
    return l


def rcp_length(dp, q, segm_length, rho, mu, valve_type=None, rho_cal=0., mu_cal=0., x=0., y=0., a_aicd=0.):
    '''
    calculates the distance between sicd's to obtain the given flow & dp on a segment.
    useful for optimizing sicd's.
    [dp]  = bar
    [q]   = m3/d
    [rho] = kg/m3
    [mu]  = cp
    valve_type: implemented: ar3, ar7 (based on oseberg FFR tests)
    rho_cal   : calibration density       (needed if valve_type is not given)
    mu_cal    : calibration viscosity     (needed if valve_type is not given)
    x         : flowrate exponent         (needed if valve_type is not given)
    y         : phase separation exponent (needed if valve_type is not given)
    a_aicd    : the linear constant       (needed if valve_type is not given)
    '''
    if not valve_type is None:
        rho_cal, mu_cal, a_aicd, x, y = rcp_params(valve_type)
    c = rcp_coeff(rho, mu, rho_cal, mu_cal, a_aicd, y)
    l = (dp/c)**(1/x) * segm_length / q
    return l

def dpoil_aicv_1(q):
    '''
    input from G:\TPD\RD_PLAB Testdata\AICD High Pressure\AICV\Testmatrix AICV 6mm Bressay oil 115 cP 2014.xlsx
    units:
    q  : m3/d
    dp : bar
    '''
    q /= 24. # m3/d --> m3/h
    return 3.4279*q**2 + 0.2892*q

def dpwater_aicv_1(q):
    '''
    input from G:\TPD\RD_PLAB Testdata\AICD High Pressure\AICV\Testmatrix AICV 6mm Bressay oil 115 cP 2014.xlsx
    q  : m3/d
    dp : bar
    '''
    q /= 24. # m3/d --> m3/h
    return 440.*q**2.5

def nozzle_area(vps, diam):
    '''
    vps  : valves pr. segm.
    diam : diameter of nozzle hole [mm]
    '''
    r = diam/2. * 1e-3  # radius in m
    return vps * pl.pi*r*r

def nozzle_dp(q, rho, diam, vps, Cv=0.686):
    '''
    q    : flow rate [m3/d] (pr. segment)
    rho  : fluid density [kg/m3]
    diam : nozzle diameter [mm]
    vps  : valves pr. segm.
    '''
    area = nozzle_area(vps, diam)
    v = q/U.DAY / area
    dp = 0.5*Cv*rho*v*v / U.BAR
    return dp

def rcp_dp4(rcp_params, bl, vl, bg, vg, mu_o, mu_g, mu_w, rho_o, rho_g, rho_w, wct, gvf, q):
    '''
    this is juan's improved dp-function
    rcp_params: (rho_cal, mu_cal, cnst, x, y)
    bl - sigmoide growth factor liquid
    vl - skewness to asymptote liquid
    bg - sigmoide growth factor gas
    vg - skewness to asymptote gas
    '''
    rho_cal, mu_cal, cnst, x, y = rcp_params
    #
    # liquid transition function
    fw_num1 = (1-pl.exp(-bl*(1-wct)))**vl*(1+pl.exp(-bl))**vl
    fw_den1 = (1+pl.exp(-bl*(1-wct)))**vl*(1-pl.exp(-bl))**vl
    fla  = 1-fw_num1/fw_den1
    #
    # gas transition function
    fg_num1 = (1-pl.exp(-bg*(1-gvf)))**vg*(1+pl.exp(-bg))**vg
    fg_den1 = (1+pl.exp(-bg*(1-gvf)))**vg*(1-pl.exp(-bg))**vg
    fga  = 1-(fg_num1/fg_den1)
    #
    # pure drop alternative formula
    dp_oil = (1-fga)*(1-fla)*(rho_cal/cnst/(rho_o**2)*(mu_o/mu_cal)**y)**(1/x)
    dp_wat = (1-fga)*(fla)*(rho_cal/cnst/(rho_w**2)*(mu_w/mu_cal)**y)**(1/x)
    dp_gas = (fga)*(rho_cal/cnst/(rho_g**2)*(mu_g/mu_cal)**y)**(1/x)
    dp = 1/((dp_oil+dp_wat+dp_gas)**(x))*q**x
    return dp
