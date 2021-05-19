'''
useful functions for analysing results from WindModeller.
uses xarray.DataArray "as much as possible"
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

SEP        = '='                    # key-value separator in WindModeller case-file (wm-file)
TYP_CT     = 'ct'                   # thrust coefficient for turbine [-]
TYP_PWR    = 'pwr'                  # power as reported by WindModeller [W]
TYP_VEL    = 'normvel'              # normalised velocities. use velocity from 'speed list' for scaling
TYP_WKL    = 'wkl'                  # wake-loss
TYP_EFF    = 'eff'                  # turbine efficency (= 100% - wakeloss%) as reported by WindModeller [-]
ARROW_SIZE = 200.
UNIT_PRCNT = '%'
UNIT_WATT  = 'W'
UNIT_NONE  = '-'

logging.basicConfig(level=logging.INFO)

def _read_case(casenm, file_ext='.wm'):
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

def get_case(case):
    if type(case) is str: case = _read_case(case)
    return case

def read_wt_positions(case=None, fnm=None):
    '''
    reads wt-positions from a WindModeller case
    - input:
      * one of case or fnm must be given.
      * case: case dictionary from get_case or casenm (str)
    - returns
      * pandas DataFrame of wt-positions
    '''
    if not fnm:
        case = get_case(case)
        fnm = case['WT location file']
    names = ['name', 'x', 'y', 'h', 'diam', 'ct_nm', 'pwr_nm']
    layout = pd.read_csv(fnm, sep=' ', skiprows=0, skipinitialspace=True, names=names, index_col='name')
    layout['name'] = layout.index    # for clarity
    return layout

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


def _get_resultsfiles(case, rtype):
    mode = 'wakes' if case['wake model'] else 'no_wakes'
    pattern = f'{case["target directory"]}/{mode}/'
    logging.info(' pattern= '+pattern)
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
    read WindModeller result-files into a DataArray with dimensions (ws, wt, wd)
    - input
      * case : case dictionary from get_case or casenm (str)
      * rtype: result-type. only power (TYP_PWR) is implemented. TODO!
    '''
    #
    case = get_case(case)
    #
    # read all results into a list (one matrix for each velocity)
    data = []
    fnms = _get_resultsfiles(case, rtype)
    assert len(fnms) > 0
    for fnm in fnms:
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
    if rtype == TYP_PWR:
        unit = UNIT_WATT
    else:
        unit = UNIT_NONE
    attrs = dict(description=case["casenm"], rtype=rtype, unit=unit)
    res = xr.DataArray(data=data, dims=dims, coords=coords, attrs=attrs)
    logging.info(f'Result dimension: {dict(res.sizes)}. Result type: {rtype}')
    return res

def _get_profilefiles(case, ptype):
    pattern = f'{case["target directory"]}/Report/Profiles/Exported_*_'
    if ptype == TYP_VEL:
        pattern += 'normvel*.csv'
    else:
        raise Exception(f'have not implemented result-type {ptype}')
    #
    fnms = UT.glob(pattern)
    fnms.sort(key=lambda x: int(x.split('speed')[1].split('.')[0]))  # sort on speed-index 1,2,3...
    return fnms

def read_profiles(case, ptype=TYP_VEL):
    '''
    read WindModeller profile output into a ????
    - input
      * case : case dictionary from get_case or casenm (str)
      * ptype: profile-type. only velocity (TYP_VEL) is implemented. TODO!
    '''
    #
    case = get_case(case)
    mode = 'wakes' if case['wake model'] else 'no_wakes'
    #
    layout = read_wt_positions(case)
    #
    # a bit tricky since i need to resample all profiles (to make indexing work)
    fnms = _get_profilefiles(case, ptype)
    hr = []       # keep profile with highest resolution
    profs = {}
    for fnm in fnms:
        logging.info(f' fnm= {fnm}')
        r = pd.read_csv(fnm, sep=',', skiprows=5, names=['vals','z'], index_col='z')
        profs[os.path.basename(fnm)] = r
        if len(r) > len(hr): hr = r
    #
    # now we can resample to highest resolution
    for fnm in profs.keys():
        profs[fnm] = r.reindex_like(hr).interpolate()
    #
    # now build the 4-dim matrix
    r4 = []
    for nm in layout.name:
        r3 = []
        for wd in case['theta list']:
            r2 = []
            for i, ws in enumerate(case['speed list']):
                fnm = f'Exported_{nm}_{ptype}_{mode}_speed{i+1}_{wd:.0f}.csv'
                assert fnm in profs
                r2.append(profs[fnm].vals.values)
            r3.append(r2)
        r4.append(r3)
    #
    # create DataArray
    layout = read_wt_positions(case)
    dims   = ["wt", "wd", "ws", "z"]
    coords = dict(wt=range(len(layout)),
                  wd=case['theta list'],
                  ws=case['speed list'],
                  z=hr.index.values,
                  x=(["wt"], layout.x.values),
                  y=(["wt"], layout.y.values),
                  nm=(["wt"], layout.name),
                 )
    attrs = dict(description=case["casenm"], ptype=ptype)
    res = xr.DataArray(data=r4, dims=dims, coords=coords, attrs=attrs)
    logging.info(f'Result dimension: {dict(res.sizes)}. Profile type: {ptype}')
    #
    return res, list(profs.values())

def calc_wakelosses(net, gross=None, wd_avr=True, ws_avr=True):
    '''
    calculating wakelosses without using the 'Offshore Array Efficiency and Effective TI' method.
    here, the net power from a wake-simulation is compared to a gross from a nowake-simulation.
    - input
      * net   : result from WindModeller case for wake-simulations.
                typically from read_results, must be DataArray
      * gross : result from WindModeller case for nowake-simulations.
                typically from read_results, must be DataArray
                if gross is None, net is compared to the max power from all turbines for given wd.
      * wd_avr: calculate per wind-direction (wd) or aggregate using mean
      * ws_avr: calculate per wind-speed (ws) or aggregate using mean
    '''
    #
    attrs = net.attrs.copy()   # lost in operations
    if ws_avr:
        net = net.mean('ws')
        if not gross is None: gross = gross.mean('ws')
    if wd_avr:
        net = net.mean('wd')
        if not gross is None: gross = gross.mean('wd')
    #
    m = gross if not gross is None else net.max()
    wl = (m-net)/m * 100
    wl.attrs = attrs
    wl.attrs['rtype'] = TYP_WKL
    wl.attrs['unit'] = UNIT_PRCNT
    return wl

def descr(res, meta=True, get_list=False, sep=' '):
    '''
    gives a short description-string of the given result.
    - input
      * res     : simulation results as DataArray
      * meta    : boolean. include some meta-info or not
      * get_list: boolean. if True: will give (meta, wt, wd, ws) as separate strings
      * sep     : string to join with
    '''
    #
    txt_meta = ''
    if meta:
        if 'description' in res.attrs:
            txt_meta += f"{res.attrs['description']}"
        if 'rtype' in res.attrs:
            txt_meta += f" {res.attrs['rtype']}"
    #
    txt_wt = f"#wt: {res.sizes['wt']}"
    #
    txt_wd = ''
    if 'wd' in res.dims:
        txt_wd += ', wd: '
        txt = res.wd.values.__str__()
        if len(res.wd) <= 5:
            txt_wd += txt
        else:
            txt_wd += ' '.join(txt.split()[:3]) +' ... ' + ' '.join(txt.split()[-1:])
    #
    txt_ws = ''
    if 'ws' in res.dims:
        txt_ws += ', ws: '
        txt = res.ws.values.__str__()
        if len(res.ws) <= 5:
            txt_ws += txt
        else:
            txt_ws += ' '.join(txt.split()[:3]) +' ... ' + ' '.join(txt.split()[-1:])
    #
    txts = (txt_meta, txt_wt, txt_wd, txt_ws)
    if get_list:
        return txts
    else:
        return sep.join(txts)

def scatterplot(res, per_wd=False, fmt='.1f', scaler=1., fontsz=13., ms=75, add_arrow=False):
    '''
    plots a scatter-plot with values per wt.
    if result (res) has multiple ws'es, they will be averaged
    - input
      * res   : simulation results as DataArray. typically wake-loss, power, or efficiency
      * per_wd: if True, it will average all directions. else, make one plot per wind-dir
      * fmt   : format of text
      * scaler: scaling (by division). typically for converting W to MW etc.
      * fontsz: size of text
      * ms    : marker size for wt's
      * add_arrow: draw an arrow showing wind-direction. only for per_wd
    '''
    #
    figsz = (7,7)
    if per_wd:
        wds = res.wd.values
        if add_arrow: figsz = (10,7)
    else:
        if 'wd' in res.dims:
            res_ = res.mean('wd', keep_attrs=True)
        else:
            res_ = res
        wds = [-1]
    #
    #
    for i, wd in enumerate(wds):
        if per_wd:
            res_ = res.sel(wd=wd)
        if 'ws' in res_.dims: res_ = res_.mean('ws', keep_attrs=True)   # could have been done already
        #
        plt.figure(figsize=figsz)
        plt.scatter(res_.x, res_.y, s=ms, c='k')
        for j, val  in enumerate(res_.values/scaler):
            plt.text(res_.x[j]+30, res_.y[j]+30, f'{val:{fmt}}', fontsize=fontsz)
        #
        if per_wd:
            txt_ws = descr(res, get_list=True)[3]
            plt.title(f"{descr(res_)} wd: {wd} {txt_ws}")
            if add_arrow:
                x = 1.2*ARROW_SIZE + res_.x.max()
                y = res_.y.max()-ARROW_SIZE
                dx, dy = WU.winddir_components(wd)
                plt.arrow(x, y, ARROW_SIZE*dx, ARROW_SIZE*dy, head_width=ARROW_SIZE/5)
                plt.xlim(right=x+2.5*ARROW_SIZE)
        else:
            plt.title(descr(res))
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
      'wind data transposition option' = Offshore Array Efficiency and Effective TI
    an efficiency matrix for each turbine is written to file.
    this routine reads this file into a DataArray for ease of access
    - input
      * case : case dictionary from get_case or casenm (str)
      * raw  : in the csv-file, it reports wd's and ws's that are slightly off the requested values.
               here, these values are replaced by the requested values (found in the case-file)
               if 'raw' is True, it will instead use the values found in the csv-file.
      * debug: get more data (for debug...)
    - returns
      * DataArray. unit is %
    '''
    case = get_case(case)
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
    wts = read_wt_positions(case)
    coords = dict(
               wt=range(len(nms)),
               wd=wd,
               ws=ws,
               x =(['wt'], wts.x.values),
               y =(['wt'], wts.y.values),
               nm=(['wt'], nms)
             )
    attrs = dict(description=case['casenm'], rtype=TYP_EFF, unit=UNIT_PRCNT)
    eff = xr.DataArray(data=vals, dims=dims, coords=coords, attrs=attrs)
    logging.info(f'Result dimension: {dict(eff.sizes)}')
    #
    #  make sure data is consistent
    assert np.all(wts.index == nms)
    #
    if debug:
        return eff, m3, m2, rec, nms
    else:
        return eff

def gross(net, eff):
    '''
    caclulate gross effect based on the reported net & turbine efficency.
    - input
      * net: net power as given by read_results()
      * eff: turbine efficency as given by read_turbine_efficencies()
    - returns
      * DataArray with gross values
    '''
    scaler = 100. if eff.attrs['unit'] == UNIT_PRCNT else 1.
    gross = net / (eff/scaler)
    gross.attrs = net.attrs.copy()
    return gross

def performance(case, full=True, plot_it=True, ms=4, ylim=None):
    '''
    get (and plot) performance-data for a given case.
    - input
      * case
      * full   : also get (and plot) conservation of variables?
      * plot_it: boolean
      * ms     : marker size
      * ylim   : ylim for imbalances (like [-40, 40]). note that imbalances are plotted in % (but data is not %)
    - returns
      * perf   : DataFrame with relevant data. column stime is solving-time
    '''
    case = get_case(case)
    tmpfile = 't'
    cmd = f"grep 'CFD Solver wall clock' {case['target directory']}/*.out > {tmpfile}"
    logging.info(cmd)
    os.system(cmd)
    perf = pd.read_csv(tmpfile, sep=' ', usecols=[0,6], header=None, names=['nm', 'stime'], converters=dict(nm=UT.basename))
    vals = [x.split('_')[-3:-1] for x in perf.nm]
    vs = [int(x[0][-1]) for x in vals]
    thetas = [int(x[1]) for x in vals]
    perf['theta'] = thetas
    perf['v'] = vs
    #
    if full:
        eqs = ['U-Mom', 'V-Mom', 'W-Mom', 'P-Mass']
        cmd = f"grep -a7 'Normalised Imbalance Summary' {case['target directory']}/*.out | grep"
        for eq in eqs: cmd += f" -e '{eq}'"
        cmd += f" > {tmpfile}"
        logging.info(cmd)
        os.system(cmd)
        df = pd.read_csv(tmpfile, sep='|', usecols=[0,1,3], header=None, converters={0:UT.basename, 1: lambda x: x.strip()})
        nms = df[df[1]==eqs[0]][0].values
        assert np.all(perf.nm.values==nms)
        #
        # populate perf-data
        for eq in eqs:
            d = df[df[1]==eq]
            perf[eq] = d[3].values
    #
    if plot_it:
        plt.figure()
        for i, v in enumerate(case['speed list']):
            p = perf[perf.v==i+1]
            plt.plot(p.theta.values, p.stime.values, 's', ms=ms, label=f'v = {v:.1f} m/s')
        plt.xlabel('Wind direction')
        plt.ylabel('Simulation runtime [s]')
        plt.title(f"Case {case['casenm']}")
        plt.legend(loc='best')
        if full:
            for eq in eqs:
                plt.figure()
                for i, v in enumerate(case['speed list']):
                    p = perf[perf.v==i+1]
                    plt.plot(p.theta.values, p[eq].values*100, 's', ms=ms, label=f'v = {v:.1f} m/s')
                plt.xlabel('Wind direction')
                plt.ylabel(f'{eq} imbalance [%]')
                plt.title(f"Case {case['casenm']}")
                plt.legend(loc='best')
                if ylim: plt.ylim(*ylim)
    #
    os.unlink(tmpfile)
    return perf

# aliases
res = read_results

