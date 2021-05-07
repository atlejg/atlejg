import numpy as np
import pandas as pd
import AtlejgTools.Utils as UT
import py_wake 

SEP     = '='
TYP_CT  = 'ct'
TYP_PWR = 'pwr'

def read_input(fnm):
    '''
    reads a WindModeller input file into a dict.
    notes
     - booleans are given as proper python booleans (True/False)
     - wind-directions are always given as a list in inp['theta list']
    '''
    inp = {}
    for line in open(fnm).readlines():
        if not SEP in line: continue
        k, v = line.split(SEP)
        inp[k.strip()] = v.strip()
    #
    # use proper booleans
    for k, v in inp.items():
        if v.lower() == 'true' : inp[k] = True
        if v.lower() == 'false': inp[k] = False
    #
    # make sure we have wind-directions on a standard format, i.e. a 'theta list' of floats
    mode = inp['wind direction option'].lower()
    if   mode == 'list':
        inp['theta list'] = [float(theta) for theta in inp['theta list']]
    elif mode == 'start end increment':
        inp['theta list'] = np.arange(float(inp['theta start']), float(inp['theta end'])+1, float(inp['theta increment']))
    elif mode == 'number':
        inp['theta list'] = np.linspace(0, 360, int(inp['number of directions']), endpoint=False)
    #
    return inp

def _write_profile(typ, xs, ys, todir, fnm):
    #
    if typ == TYP_CT:
        nm1, nm2, unit = 'CThrust', 'cthrust', ' '
    else:
        nm1, nm2, unit = 'PowerC', 'power', 'W'
    txt = f'''
[Name],
{nm1},
[Spatial Fields],
speed,
[Data],
speed [m s^-1], {nm2} [{unit}]
'''.lstrip()
    for x,y in zip(xs, ys):
        txt += f'{x:5.1f}, {y:.3f}\n'
    f = open(f'{todir}/{fnm}', 'w')
    f.write(txt)
    f.close()
    print(f.name, 'was created')
    return f.name

def convert_wtg_file(fnm, todir='.'):
    '''
    write power & thrust profiles in WindModeller format.
    note: assumes units are W & 1
    - input
      fnm : a wtg-file for one turbine
    - returns
      a py_wake.wind_turbines.WindTurbines object
    '''
    wtg = py_wake.wind_turbines.WindTurbines.from_WAsP_wtg(fnm)
    ucp = wtg.upct_tables[0]   # assume it's only one
    base = UT.clean_string_for_fname(wtg.name())
    fnm_ct  = _write_profile(TYP_CT,  ucp[:,0], ucp[:,1], todir, f'{base}_ct.txt')
    fnm_pwr = _write_profile(TYP_PWR, ucp[:,0], ucp[:,2], todir, f'{base}_pwr.txt')
    return fnm_ct, fnm_pwr, wtg

def read_wt_positions(inp=None, fnm=None):
    '''
    reads wt-positions from a WindModeller file.
    one of inp or fnm must be given.
    inp is a dict as given by read_input
    '''
    fnm = inp['WT location file'] if inp else fnm
    return pd.read_csv(fnm, sep=' ', skiprows=0, skipinitialspace=True, names=['name', 'x', 'y', 'h'], index_col='name')

def relative_wt_positions(wts):
    '''
    when running offshore cases (without terrain-data) it seems we *MUST*
    have wtg-positions around origo
    input:
        wts: as read by read_wt_positions
    '''
    wts2 = wts.copy()
    wts2.x = wts.x - wts.x.mean()
    wts2.y = wts.y - wts.y.mean()
    return wts2

def write_wt_positions(wtgs, fnm):
    wtgs.to_csv(fnm, sep=' ', float_format='%.1f', header=False)

def read_results(inp, rtype=TYP_PWR, ws_ix=0):
    wakedir = 'wakes' if inp['wake model'] else 'no_wakes'
    fnm = f'{inp["target directory"]}/{wakedir}/'
    if rtype == TYP_PWR:
        fnm += 'WT_hub_power_'
    fnm += f'speed{ws_ix+1:d}.csv'
    return pd.read_csv(fnm, sep=',', skiprows=3, index_col=0, header=None)
