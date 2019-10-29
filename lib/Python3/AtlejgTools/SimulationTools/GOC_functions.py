from scipy.interpolate import interp1d
import AtlejgTools.Utils as UT
import numpy as np

'''
Some useful (simple?) Gas-Oil-Contact functions.
Note1:
 If a file of GOC-observations is expected, the format is like this:
    MD [m]   TVD [m]   +/- [m]
   2265.0    1550.8     0.3
   2985.0    1553.7     0.3
   6367.0    1552.0     0.3
'''

from pylab import *
from scipy.optimize import curve_fit

def goc_sinusoidal(md, goc0, wave_l, ampl, shift):
    return goc0 + ampl*sin(2*pi/wave_l*md - shift)

def goc_sin1(wave_l, ampl, fname):
    '''
    Creates a sinusoidal curve with given wavelength and amplitude
    that gives the best fit to observations.
    Note: This function does not use the uncertainty given.
    File of GOC observations expected. Format is given in Note1
    '''
    mds, gocs, dgcocs = loadtxt(fname, skiprows=1).T
    goc_func = lambda md, goc0, shift: goc_sinusoidal(md, goc0, wave_l, ampl, shift)
    goc0, shift = curve_fit(goc_func, mds, gocs)[0]
    return lambda md: goc_sinusoidal(md, goc0, wave_l, ampl, shift)

def goc_linear1(xlsfile, sheet, md_col, goc_col, skip=1, last_col=99999):
    '''
    returns a function that is a linear interpolation of the entries in the provided excel.
    input:
     xlsfile : file name
     sheet   : sheet name
     md_col  : name of column with measured depth (on Troll: typically 'A')
     goc_col : name of column with GOC            (on Troll: typically 'DF')
    '''
    mds = np.array([float(x) for x in UT.read_excel_column(sheet, md_col, skip+1, last_col, file=xlsfile) if is_numlike(x)])
    gocs = UT.read_excel_column(sheet, goc_col, skip+1, last_col, file=xlsfile)
    ixs = find(list(map(is_numlike, gocs)))
    if len(ixs) == 0: raise Exception('No GOC-values found')
    if len(ixs) == 1: return lambda md: gocs[ixs[0]]
    mds  = np.array(mds)[ixs]
    gocs = np.array(gocs)[ixs]
    return interp1d(mds, gocs, fill_value='extrapolate')

def goc_linear2(fname):
    '''
    returns a function that is a linear interpolation of points given in file
    input:
     fname : file name. format is given in Note1
    '''
    mds, gocs, dgcocs = loadtxt(fname, skiprows=1).T
    return interp1d(mds, gocs, fill_value='extrapolate')

if __name__ == '__main__':

    wave_l = float(sys.argv[1])
    ampl   = float(sys.argv[2])
    fname  = sys.argv[3]

    md_meas, goc_meas, dgcocs = loadtxt(fname, skiprows=1).T

    mds = linspace(md_meas[0], md_meas[-1], 100)

    goc_func = goc_sin1(wave_l, ampl)
    gocs     = goc_func(mds)

    figure()
    plot(mds, -gocs, '--')
    plot(md_meas, -goc_meas, 'o')
    grid(1)
    ylim(1.01*min(-gocs), .99*max(-gocs))

    show()
