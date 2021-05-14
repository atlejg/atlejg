'''
useful functions for analysing results from WindModeller.
uses xarray.DataArray for this.
'''

__author__  = 'atle j. gyllensten'
__version__ = '0.1'
__email__   = 'agy@equinor.com'
__date__    = '2021-05'

import os
import logging
import numpy as np
import pandas as pd
import xarray as xr
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

def read_case(casenm, file_ext='.wm'):
    '''
    reads a WindModeller case-file (typically a00.wm) into a dict.
    - input
      * casenm: this is the basename of the file-name
    - notes
     * booleans are given as proper python booleans (True/False)
     * wind-directions are always given as a list in case['theta list'] as np.array
     * wind-speeds is np.array
     * adds 'casenm' for identification. this is the basename of the file-name
    '''
    #
    # make sure casenm is without extension and fnm is with
    casenm = UT.basename(casenm)
    fnm = casenm + file_ext
    if not os.path.exists(fnm):
        raise Exception(f'No such file: {fnm}')
    #
    # read keys / values into dict
    case = {}
    for line in open(fnm).readlines():
        if not SEP in line: continue
        key, val = line.split(SEP)
        key, val = key.strip(), val.strip()
        # use proper booleans
        if   val.lower() == 'true' : val = True
        elif val.lower() == 'false': val = False
        case[key] = val
    #
    # make sure we have wind-directions on a standard format, i.e. a 'theta list' of floats
    mode = case['wind direction option'].lower()
    if mode == 'list':
        case['theta list'] = np.array([float(theta) for theta in case['theta list'].split(',')])
    elif mode == 'start end increment':
        case['theta list'] = np.arange(float(case['theta start']), float(case['theta end'])+1, float(case['theta increment']))
    elif mode == 'number':
        case['theta list'] = np.linspace(0, 360, int(case['number of directions']), endpoint=False)
    #
    # want speed list as floats
    case['speed list'] = np.array([float(speed) for speed in case['speed list'].split(',')])
    #
    # adds 'casenm' for identification
    if not 'casenm' in case: case['casenm'] = UT.basename(fnm)
    return case

def write_profile(typ, xs, ys, todir, fnm):
    '''
    writes a pair of vectors to file with a format accepted by WindModeller
    '''
    #
    if typ == TYP_CT:
        nm1, nm2, unit = 'CThrust', 'cthrust', ' '
    elif typ == TYP_PWR:
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
    note: assumes power unit is W (thrust is dimensionless always, i think)
    - input
      * fnm : a wtg-file for one turbine in WAsP-format (which is xml)
    - returns
      * a py_wake.wind_turbines.WindTurbines object
    '''
    wtg = py_wake.wind_turbines.WindTurbines.from_WAsP_wtg(fnm)
    ucp = wtg.upct_tables[0]   # assume it's only one
    base = UT.clean_string_for_fname(wtg.name())
    fnm_ct  = write_profile(TYP_CT,  ucp[:,0], ucp[:,1], todir, f'{base}_ct.txt')
    fnm_pwr = write_profile(TYP_PWR, ucp[:,0], ucp[:,2], todir, f'{base}_pwr.txt')
    return fnm_ct, fnm_pwr, wtg

def read_wt_positions(case=None, fnm=None):
    '''
    reads wt-positions from a WindModeller case
    - input:
      * one of case or fnm must be given.
      * case: case dictionary from read_case or casenm (str)
    - returns
      * pandas DataFrame of wt-positions
    '''
    if not fnm:
        if type(case) is str: case = read_case(case)
        fnm = case['WT location file']
    names = ['name', 'x', 'y', 'h', 'diam', 'ct_nm', 'pwr_nm']
    return pd.read_csv(fnm, sep=' ', skiprows=0, skipinitialspace=True, names=names, index_col='name')

def relative_wt_positions(wts):
    '''
    when running offshore cases (without terrain-data) it seems we *MUST*
    have wt-positions around origo
    input:
     * wts: as read by read_wt_positions
    '''
    wts = wts.copy()
    wts.x = wts.x - wts.x.mean()
    wts.y = wts.y - wts.y.mean()
    return wts

def write_wt_positions(wtgs, fnm):
    '''
    write wt-positions to a csv-file
    - input:
      * wtgs: pandas DataFrame of wt-positions
    '''
    wtgs.to_csv(fnm, sep=' ', float_format='%.1f', header=False)

def get_resultsfiles(case, rtype):
    wakedir = 'wakes' if case['wake model'] else 'no_wakes'
    pattern = f'{case["target directory"]}/{wakedir}/'
    if rtype == TYP_PWR:
        pattern += 'WT_hub_power_'
    else:
        raise Exception(f'have not implemented result-type {rtype}')
    #
    pattern += f'speed*.csv'
    fnms = UT.glob(pattern)
    fnms.sort(key=lambda x: int(x.split('speed')[1].split('.')[0]))  # sort on speed-index 1,2,3...
    return fnms

def read_results(case, rtype=TYP_PWR):
    '''
    read WindModeller result-files into a DataArray
    - input
      * case : case dictionary from read_case or casenm (str)
      * rtype: result-type. only power (TYP_PWR) is implemented. TODO!
    '''
    #
    if type(case) is str: case = read_case(case)
    #
    # read all results into a list (one matrix for each velocity)
    data = []
    for fnm in get_resultsfiles(case, rtype):
        logging.info(f'reading {fnm}')
        res = pd.read_csv(fnm, sep=',', skiprows=3, index_col=0, header=None)
        data.append(res.values)
    #
    # create DataArray
    layout = read_wt_positions(case)
    dims   = ["ws", "wt", "wd"]
    coords = dict(ws=case['speed list'],
                  wt=range(len(layout)),
                  wd=case['theta list'],
                  x=(["wt"], layout.x.values),
                  y=(["wt"], layout.y.values),
                  nm=(["wt"], res.index),
                 )
    attrs = dict(description=case["casenm"], rtype=rtype)
    res = xr.DataArray(data=data, dims=dims, coords=coords, attrs=attrs)
    logging.info(f'Result dimension: {dict(res.sizes)}. Result type: {rtype}')
    return res


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

def wtg_scatterplot(sim0, per_wd=False, descr='', fmt='.1f', scaler=1., fontsz=13., ms=75):
    '''
    plots a scatter-plot with values per wtg.
    - input
      * sim0  : simulation results as DataArray. typically wake-loss, power, or efficiency
      * per_wd: if True, it will average all directions
                else, it will make one plot per wind-dir
      * descr: title prefix
      * fmt   : format of text
      * scaler: scaling (by division). typically for converting W to MW etc.
      * fontsz: size of text
      * ms    : marker size for wt's
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
        if 'ws' in sim.dims: sim = sim.mean('ws')   # could have been done already
        #
        plt.figure(figsize=figsz)
        plt.scatter(sim.x, sim.y, s=ms, c='k')
        for j, val  in enumerate(sim.values/scaler):
            plt.text(sim.x[j]+30, sim.y[j]+30, f'{val:{fmt}}', fontsize=fontsz)
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

def cp(old, new, file_ext='.wm'):
    '''
    copy an existing WindModeller case-file to a new one
    and also update the 'target directory'
     - input
       * old: existing case without file extension
              must be in current working directory
       * new: new case without file extension
     - returns
       * name of new case-file. will be put in current working directory
    '''
    if not os.path.exists(f'{old}{file_ext}'):
        print(f'No such file {old}{file_ext}. Returning')
        return
    if os.path.exists(f'{new}{file_ext}'):
        print(f'File {new}{file_ext} exists. Returning')
        return
    cmd = f'sed s/{old}/{new}/g {old}{file_ext} > {new}{file_ext}'
    UT.tcsh(cmd)
    logging.info(cmd)
    return new + file_ext

def write_winddata_file(scd, fnm, std_lvl=0.1, stab='U'):
    '''
    writes scada data into a WindModeller file for 'Wind Data Files' as required when
    using the 'Offshore Array Efficiency and Effective TI' option
    - input
      * scd     : scada data. either as WU.Scada object or a DataFrame
      * std_lvl : standard deviation level
      * stab   : stability regime (U[nstable], W[ery]U[nstable], S[table], ..)
    '''
    if not isinstance(scd, pd.DataFrame):
        scd = scd.data
    if not 'stddev' in scd: scd['stddev'] = scd.WindSpeed * std_lvl
    if not 'stab'   in scd: scd['stab']   = stab
    f = open(fnm, 'w')
    f.write('WFU_data,,,,\nDirection,Mean,SD,TIMESTAMP,Stability\n')
    cols = ['WindDir', 'WindSpeed', 'stddev', 'time', 'stab']
    scd.to_csv(f, sep=',', columns=cols, index=False, header=False, float_format='%.2f', date_format='%d/%m/%Y %H:%M')
    logging.info(f'{fnm} is created')

def read_turbine_efficencies(case, raw=False, debug=False):
    '''
    for a WindModeller-run of type
      wind data transposition option = Offshore Array Efficiency and Effective TI
    a efficiency matrix for each turbine is written to file.
    this routine reads this file into a DataArray for ease of access
    - input
      * case : case dictionary from read_case or casenm (str)
      * raw  : in the csv-file, it reports wd's and ws's that are slightly off the requested values.
               here, these values are replaced by the requested values (found in the case-file)
               if 'raw' is True, it will instead use the values found in the csv-file.
      * debug: get more data (for debug...)
    '''
    if type(case) is str: case = read_case(case)
    fnm = f'{case["target directory"]}/array_eff_TI/results/TurbineEfficiency.csv'
    logging.info(f'reading {fnm}')
    #
    csvfile = open(fnm)
    m3 = []                                     # 3-dim matrix for all data
    nms = []
    for line in csvfile:
        line = line.strip()
        if (not line) or ('Upstream' in line) or ('Turbine efficiency' in line):
            continue
        #
        if 'Turbine:' in line:
            if len(nms) > 0:
                m3.append(m2)
            m2 = []                             # 2-dim matrix of data for this turbine
            nm = line.split('Turbine:,')[1].strip()
            nms.append(nm)
            continue
        #
        # 'default': read data into current matrix
        rec = [float(x) for x in line.split(',')]
        m2.append(rec)
    m3.append(m2)
    csvfile.close()
    #
    # build a DataArray
    da = xr.DataArray(m3)  # convinent
    if not raw:
        wd = case['theta list']
        ws = case['speed list']
    else:
        # a bit tricky ...
        wd = da[0,:,0].values[:len(case['theta list'])]
        ws = da[0,:,1].values[::len(case['theta list'])]
    vals = da[:,:,2].values.reshape([len(nms), len(wd), len(ws)])
    dims   = ['wt', 'wd', 'ws']
    wtgs = read_wt_positions(case)
    coords = dict(
               wt=range(len(nms)),
               wd=wd,
               ws=ws,
               x =(['wt'], wtgs.x.values),
               y =(['wt'], wtgs.y.values),
               nm=(['wt'], nms)
             )
    eff = xr.DataArray(data=vals, dims=dims, coords=coords, attrs=dict(description=case['casenm']))
    logging.info(f'Result dimension: {dict(eff.sizes)}')
    #
    #  make sure data is consistent
    assert np.all(wtgs.index == nms)
    #
    if debug:
        return eff, m3, m2, rec, nms
    else:
        return eff


# aliases
res = read_results

