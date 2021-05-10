import os
import logging
import numpy as np
import pandas as pd
import xarray
import AtlejgTools.Utils as UT
import py_wake 
import matplotlib.pyplot as plt
import AtlejgTools.WindResourceAssessment.Utils as WU

SEP        = '='
TYP_CT     = 'ct'
TYP_PWR    = 'pwr'
TYP_WL     = 'wl'
ARROW_SIZE = 200.

logging.basicConfig(level=logging.INFO)

def read_input(casenm, file_ext='.wm'):
    '''
    reads a WindModeller input file into a dict.
    notes
     - booleans are given as proper python booleans (True/False)
     - wind-directions are always given as a list in inp['theta list'] as floats
     - wind-speeds are floats
     - adds 'casenm' for identification. this is the basename of the file-name
    '''
    inp = {}
    fnm = casenm + file_ext
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
        inp['theta list'] = [float(theta) for theta in inp['theta list'].split(',')]
    elif mode == 'start end increment':
        inp['theta list'] = np.arange(float(inp['theta start']), float(inp['theta end'])+1, float(inp['theta increment']))
    elif mode == 'number':
        inp['theta list'] = np.linspace(0, 360, int(inp['number of directions']), endpoint=False)
    #
    # want speed list as floats
    inp['speed list'] = [float(speed) for speed in inp['speed list'].split(',')]
    #
    # adds 'casenm' for identification
    if not 'casenm' in inp: inp['casenm'] = UT.basename(fnm)
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
    inp is a dict as given by read_input or casenm
    '''
    if not fnm:
        if type(inp) is str: inp = read_input(inp)
        fnm = inp['WT location file']
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

def get_resultsfiles(inp, rtype):
    wakedir = 'wakes' if inp['wake model'] else 'no_wakes'
    pattern = f'{inp["target directory"]}/{wakedir}/'
    if rtype == TYP_PWR:
        pattern += 'WT_hub_power_'
    else:
        raise Exception(f'have not implemented result-type {rtype}')
    #
    pattern += f'speed*.csv'
    fnms = UT.glob(pattern)
    fnms.sort(key=lambda x: int(x.split('speed')[1].split('.')[0]))  # sort on speed-index 1,2,3...
    return fnms

def read_results(inp, rtype=TYP_PWR):
    '''
    read WindModeller result-files into an xarray
    - input
        inp: case dictionary as read from read_input or casenm
    '''
    #
    if type(inp) is str: inp = read_input(inp)
    #
    # read all results into a list (one matrix for each velocity)
    data = []
    for fnm in get_resultsfiles(inp, rtype):
        logging.info(f'reading {fnm}')
        res = pd.read_csv(fnm, sep=',', skiprows=3, index_col=0, header=None)
        data.append(res.values)
    #
    # create xarray
    layout = read_wt_positions(inp)
    dims   = ["ws", "wt", "wd"]
    coords = dict(ws=inp['speed list'],
                  wt=range(len(layout)),
                  wd=inp['theta list'],
                  x=(["wt"], layout.x.values),
                  y=(["wt"], layout.y.values),
                  nm=(["wt"], res.index),
                 )
    return xarray.DataArray(data=data, dims=dims, coords=coords, attrs=dict(description=inp["casenm"]))

def res(inp, rtype=TYP_PWR):
    '''
    just an alias for read_results
    '''
    return read_results(inp, rtype)

def calc_wakelosses(net, gross=None, per_wd=False):
    #
    attrs = net.attrs     # lost in operations
    net = net.mean('ws')
    if not per_wd:
        net = net.mean('wd')
        if not gross is None: gross = gross.mean('wd')
    #
    m = gross if not gross is None else net.max()
    wl = (m-net)/m * 100
    wl.attrs = attrs
    return wl

def wtg_scatterplot(sim0, per_wd=False, descr='', fmt='.1f', scaler=1.):
    '''
    plots a scatter-plot with values per wtg.
    - input
        sim0  : simulation results. typically wake-loss or power
        per_wd: if True, it will average all directions
                else, it will make one plot per wind-dir
        descr: title prefix
    '''
    #
    descr += sim0.attrs["description"]
    if per_wd:
        wds = sim0.wd
        figsz = (10,7)
    else:
        if 'wd' in sim0.dims:
            sim = sim0.mean('wd')
        else:
            sim = sim0
        wds = [-1]
        figsz = (7,7)
    #
    for i, wd in enumerate(wds):
        if per_wd:
            sim = sim0.sel(wd=wd)
        if 'ws' in sim.dims: sim = sim.mean('ws')
        #
        plt.figure(figsize=figsz)
        plt.scatter(sim.x, sim.y, s=50, c='k')
        for j, val  in enumerate(sim.values/scaler):
            plt.text(sim.x[j]+30, sim.y[j]+30, f'{val:{fmt}}')
        #
        titl = descr
        if per_wd:
            # add wind-arrow. make extra space for it
            x = 1.2*ARROW_SIZE + sim.x.max()
            y = sim.y.max()-ARROW_SIZE
            dx, dy = WU.winddir_components(wd)
            plt.arrow(x, y, ARROW_SIZE*dx, ARROW_SIZE*dy, head_width=ARROW_SIZE/5)
            plt.xlim(right=x+2.5*ARROW_SIZE)
            titl += (r' $[%d^o]$'%wd)
        else:
            titl += ' [all directions]'
        plt.title(titl)
        plt.xticks([])
        plt.yticks([])
        plt.axis('equal')

def cp(old, new):
    '''
    copy an existing WindModeller input-file to a new one
    and also update the 'target directory'
     - input
         old: existing case without file extension
         new: new case without file extension
     - notes
         - file must be in current working directory
         - old file must have file extension wm
    '''
    if not os.path.exists(f'{old}.wm'):
        print(f'No such file {old}.wm. Returning')
        return
    if os.path.exists(f'{new}.wm'):
        print(f'File {new}.wm exists. Returning')
        return
    cmd = f'sed s/{old}/{new}/g {old}.wm > {new}.wm'
    UT.tcsh(cmd)
