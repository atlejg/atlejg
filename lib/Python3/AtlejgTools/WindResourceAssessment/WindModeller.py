'''
useful functions for analysing results from WindModeller.
uses xarray.DataArray "as much as possible"
WindModellerCase is the main part.
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
from scipy.interpolate import interp1d

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

def read_params(casenm, file_ext='.wm'):
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
    # read keys / values into paramter-dict
    pm = {}
    for line in open(fnm).readlines():
        if not SEP in line: continue
        key, val = line.split(SEP)
        key, val = key.strip(), val.strip()
        if ',' in val:
            # split comma-separated values into list
            val = [rec.strip() for rec in val.split(',')]
        else:
            if   val.lower() == 'true' : val = True       # use proper booleans
            elif val.lower() == 'false': val = False
        #
        pm[key] = val
    # make sure we have wind-directions on a standard format, i.e. a 'theta list' of floats
    mode = pm['wind direction option'].lower()
    if mode == 'list':
        pm['theta list'] = np.array([float(theta) for theta in pm['theta list']])
    elif mode == 'start end increment':
        pm['theta list'] = np.arange(float(pm['theta start']), float(pm['theta end'])+1, float(pm['theta increment']))
    elif mode == 'number':
        pm['theta list'] = np.linspace(0, 360, int(pm['number of directions']), endpoint=False)
    #
    # want speed list as floats
    if type(pm['speed list']) is list:
        pm['speed list'] = np.array([float(speed) for speed in pm['speed list']])
    else:
        pm['speed list'] = np.array([float(pm['speed list'])])
    #
    # want file-lists as lists
    if not type(pm['power files']) is list:
        pm['power files'] = [pm['power files']] 
    if not type(pm['thrust coefficient files']) is list:
        pm['thrust coefficient files'] = [pm['thrust coefficient files']] 
    #
    # add 'casenm' for identification
    if not 'casenm' in pm: pm['casenm'] = casenm
    return pm

def write_curve(typ, xs, ys, todir, fnm):
    '''
    writes a pair of vectors to file with a format accepted as an input-curve for WindModeller
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
    fnm_ct  = write_curve(TYP_CT,  ucp[:,0], ucp[:,1], todir, f'{base}_ct.txt')
    fnm_pwr = write_curve(TYP_PWR, ucp[:,0], ucp[:,2], todir, f'{base}_pwr.txt')
    return fnm_ct, fnm_pwr, wtg

def read_curve(fnm, as_func):
    '''
    reads a WindModeller-format input-curve into a DataFrame or an interp1d-function
      - input:
        * fnm: file name
        * as_func: boolean. get interp1d-function or the raw data as a DataFrame
    - returns
      * curve (function or DataFrame)
      * name of the curve
      * unit of variables
    '''
    if not os.path.exists(fnm):
        raise Exception(f"cannot find file {fnm}")
    curve = pd.read_csv(fnm, sep=',', header=None, skiprows=6)
    if as_func:
        curve = interp1d(curve[0], curve[1])
    lines = open(fnm).readlines()
    nm = lines[1].strip().replace(',','')
    units = [rec.strip() for rec in lines[5].split(',')]
    return curve, nm, units

def relative_wt_positions(layout):
    '''
    when running offshore cases (without terrain-data) it seems we *MUST*
    have wt-positions around origo
    input:
     * layout: as read by read_layout
    '''
    layout = layout.copy()
    layout.x = layout.x - layout.x.mean()
    layout.y = layout.y - layout.y.mean()
    return layout

def write_wt_positions(layout, fnm):
    '''
    write wt-positions to a csv-file
    - input:
      * layout: pandas DataFrame of wt-positions
    '''
    layout.to_csv(fnm, sep=' ', float_format='%.1f', header=False)

def calc_wakelosses(net, gross=None, wd_avr=True, ws_avr=True):
    '''
    calculating wakelosses without using the 'Offshore Array Efficiency and Effective TI' method.
    - input
      * net   : result from WindModeller case for wake-simulations.
                typically from _read_results, must be DataArray
      * gross : result from WindModeller case for nowake-simulations.
                typically from _read_results, must be DataArray
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
    if 'wd' in res.coords:
        txt_wd += ', wd: '
        txt = res.wd.values.__str__()
        if 'wd' in res.dims and len(res.wd) > 5:
            txt_wd += ' '.join(txt.split()[:3]) +' ... ' + ' '.join(txt.split()[-1:])
        else:
            txt_wd += txt
    #
    txt_ws = ''
    if 'ws' in res.coords:
        txt_ws += ', ws: '
        txt = res.ws.values.__str__()
        if 'ws' in res.dims and len(res.ws) > 5:
            txt_ws += ' '.join(txt.split()[:3]) +' ... ' + ' '.join(txt.split()[-1:])
        else:
            txt_ws += txt
    #
    txts = (txt_meta, txt_wt, txt_wd, txt_ws)
    if get_list:
        return txts
    else:
        return sep.join(txts)

def scatterplot(res, per_wd=False, fmt='.1f', scaler=1., fontsz=13., ms=75, add_arrow=False, nm_only=False, tight=True):
    '''
    plots a scatter-plot with values per wt.
    if result (res) has multiple ws'es, they will be averaged
    - input
      * res   : simulation results as DataArray
      * per_wd: if True, it will average all directions. else, make one plot per wind-dir
      * fmt   : format of numbers to be displayed
      * scaler: scaling (by multiplying). typically for converting W to MW etc.
      * fontsz: size of text
      * ms    : marker size for wt's
      * add_arrow: draw an arrow showing wind-direction. only for per_wd
      * nm_only: boolean. show name of turbine only (not a number)
      * tight : boolean. do tight_layout?
    '''
    #
    figsz = (7,7)
    if per_wd:
        wds = res.wd.values
        if add_arrow: figsz = (10,7)
    else:
        if 'wd' in res.dims:
            res2 = res.mean('wd', keep_attrs=True)
        else:
            res2 = res
        wds = [-1]
    #
    #
    for i, wd in enumerate(wds):
        if per_wd:
            res2 = res.sel(wd=wd)
        if 'ws' in res2.dims: res2 = res2.mean('ws', keep_attrs=True)   # could have been done already
        #
        plt.figure(figsize=figsz)
        plt.scatter(res2.x, res2.y, s=ms, c='k')
        for j, val  in enumerate(res2.values):
            if nm_only:
                txt = res.nm.values[j]
            else:
                txt = f'{val*scaler:{fmt}}'
            plt.text(res2.x[j]+30, res2.y[j], txt, fontsize=fontsz)
        #
        if per_wd:
            txt_ws = descr(res, get_list=True)[3]
            plt.title(f"{descr(res2)} wd: {wd} {txt_ws}")
            if add_arrow:
                x = 1.2*ARROW_SIZE + res2.x.max()
                y = res2.y.max()-ARROW_SIZE
                dx, dy = WU.winddir_components(wd)
                plt.arrow(x, y, ARROW_SIZE*dx, ARROW_SIZE*dy, head_width=ARROW_SIZE/5)
                plt.xlim(right=x+2.5*ARROW_SIZE)
        else:
            plt.title(descr(res))
        plt.xticks([])
        plt.yticks([])
        plt.axis('equal')
        if tight: plt.tight_layout()

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

def plot(res, ws, title='', wts=None, kwargs={'ls':'-','lw':3}, ylim=None, onefig=True, unit=''):
    '''
    plot plot(s) for the given data. typically used for plotting wakeloss per direction.
    - input
      * res   : simulation results as DataArray
      * ws    : wind-speed to use
      * title : prepend this title
      * wts   : which wt's to plot for (indices)
      * kwargs: plotting key-words
      * ylim  : [ymin, ymax]
      * onefig: boolean. one common figure, or one per wind-turbine
      * unit  : unit for yticks. will use attrs['unit'] if unit is not given
    '''  
    res = res.sel(ws=ws)
    if wts is None:
        wts = res.wt.values
    if onefig: plt.figure(figsize=(7,7))
    for i in wts:
        if not onefig: plt.figure(figsize=(7,7))
        y = res.sel(wt=i).values.ravel()
        plt.plot(res.wd, y, **kwargs, label=res.nm.values[i])
        plt.title(f"{title} ws = {ws}m/s") 
        plt.tight_layout()
        plt.legend(loc='best')
        if not ylim is None: plt.ylim(ylim)
        if not unit and 'unit' in res.attrs:
            unit = res.attrs['unit']
        if unit:
            yt = plt.yticks()[0]
            plt.yticks(yt, [f"{y:.0f}{unit}" for y in yt])
    plt.show()

def polarplot(res, ws, title='', wts=None, kwargs={'ls':'-','lw':3}, ylim=None, onefig=True, unit=''):
    '''
    plot polar-plot(s) for the given data. typically used for plotting wakeloss per direction.
    - notes
      * note1: plt.polar() defines the angle counter-clockwise from east (positve x-axis), but
               we want clockwise from north (positve y-axis). this is handled, but maybe not so elegantly..
    - input
      * res   : simulation results as DataArray
      * ws    : wind-speed to use
      * title : prepend this title
      * wts   : which wt's to plot for (indices)
      * kwargs: plotting key-words
      * ylim  : [ymin, ymax]
      * onefig: boolean. one common figure, or one per wind-turbine
      * unit  : unit for yticks. will use attrs['unit'] if unit is not given
    '''  
    res = res.sel(ws=ws)
    if wts is None:
        wts = res.wt.values
    if onefig: plt.figure(figsize=(7,7))
    theta = -res.wd.values/180*np.pi + np.pi/2    # negative and shift 90deg to get it right (see note1)
    for i in wts:
        if not onefig: plt.figure(figsize=(7,7))
        y = res.sel(wt=i).values.ravel()
        plt.polar(theta, y, **kwargs, label=res.nm.values[i])
        plt.title(f"{title} ws = {ws}m/s") 
        plt.tight_layout()
        plt.legend(loc='best')
        if not ylim is None: plt.ylim(ylim)
        if not unit and 'unit' in res.attrs:
            unit = res.attrs['unit']
        if unit:
            yt = plt.yticks()[0]
            plt.yticks(yt, [f"{y:.0f}{unit}" for y in yt])
        plt.xticks(plt.xticks()[0], [])           # remove xticks (see note1)
    plt.show()

class WindModellerCase(object):
#
    def __init__(self, name):
        self.pm   = read_params(name)
        # convinent stuff
        self.name   = self.pm['casenm']
        self.layout = self.read_layout()
        self.wd     = self.pm['theta list']
        self.ws     = self.pm['speed list']
        self.wt     = self.layout.name.values
#
    def get_curve(self, typ, label='', as_func=True):
        '''
        get the power-curve or CT-curve used in the given case.
        if more than one curve exists, you will get the first curve or you need to give the name of the curve.
        - input
          * typ: TYP_PWR or TYP_CT
          * label : name of power-curve. if not given, the first curve will be returned
          * as_func: boolean. get interp1d-function or the raw data as a DataFrame
        - returns
          * pwr (function or DataFrame)
        '''
        #
        if typ == TYP_PWR:
            fnms = self.pm['power files']
        elif typ == TYP_CT:
            fnms = self.pm['thrust coefficient files']
        else:
            raise Exception(f'No such type of curve available: {typ}')
        #
        if not label:
            return read_curve(fnms[0], as_func)
        for fnm in fnms:
            curve, nm, units = read_curve(fnm, as_func)
            if nm == label: return curve, nm, units
        #
        # should not get here
        raise Exception(f"cannot find curve {label}")
#
    def read_layout(self):
        '''
        reads wt-positions from a this case
        - returns
          * pandas DataFrame of wt-positions
        '''
        names = ['name', 'x', 'y', 'h', 'diam', 'ct_nm', 'pwr_nm']
        fnm = self.pm['WT location file']
        layout = pd.read_csv(fnm, sep=' ', skiprows=0, skipinitialspace=True, names=names, index_col='name')
        layout['name'] = layout.index    # for clarity
        return layout
#
    def _get_resultsfiles(self, rtype):
        mode = 'wakes' if self.pm['wake model'] else 'no_wakes'
        pattern = f'{self.pm["target directory"]}/{mode}/'
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
#
    def _read_results(self, rtype=TYP_PWR):
        '''
        read WindModeller result-files into a DataArray with dimensions (ws, wt, wd)
        - input
          * rtype: result-type. only power (TYP_PWR) is implemented. TODO?
        '''
        #
        # read all results into a list (one matrix for each velocity)
        data = []
        fnms = self._get_resultsfiles(rtype)
        assert len(fnms) > 0
        for fnm in fnms:
            logging.info(f'reading {fnm}')
            res = pd.read_csv(fnm, sep=',', skiprows=3, index_col=0, header=None)
            data.append(res.values)
        #
        # create DataArray
        dims   = ["ws", "wt", "wd"]
        coords = dict(ws=self.ws,
                      wt=range(len(self.layout)),
                      wd=self.wd,
                      x=(["wt"], self.layout.x.values),
                      y=(["wt"], self.layout.y.values),
                      nm=(["wt"], res.index),
                     )
        if rtype == TYP_PWR:
            unit = UNIT_WATT
        else:
            unit = UNIT_NONE
        attrs = dict(description=self.name, rtype=rtype, unit=unit)
        res = xr.DataArray(data=data, dims=dims, coords=coords, attrs=attrs)
        logging.info(f'Result dimension: {dict(res.sizes)}. Result type: {rtype}')
        return res
#
    def _get_profilefiles(self, ptype):
        pattern = f'{self.pm["target directory"]}/Report/Profiles/Exported_*_'
        if ptype == TYP_VEL:
            pattern += 'normvel*.csv'
        else:
            raise Exception(f'have not implemented result-type {ptype}')
        #
        fnms = UT.glob(pattern)
        fnms.sort(key=lambda x: int(x.split('speed')[1].split('.')[0]))  # sort on speed-index 1,2,3...
        return fnms
#
    def read_profiles(self, ptype=TYP_VEL):
        '''
        read WindModeller profile output into a ????
        - input
          * ptype: profile-type. only velocity (TYP_VEL) is implemented. TODO!
        '''
        #
        mode = 'wakes' if self.pm['wake model'] else 'no_wakes'
        #
        #
        # a bit tricky since i need to resample all profiles (to make indexing work)
        fnms = self._get_profilefiles(ptype)
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
        for nm in self.wt:
            r3 = []
            for wd in self.wd:
                r2 = []
                for i, ws in enumerate(self.ws):
                    fnm = f'Exported_{nm}_{ptype}_{mode}_speed{i+1}_{wd:.0f}.csv'
                    assert fnm in profs
                    r2.append(profs[fnm].vals.values)
                r3.append(r2)
            r4.append(r3)
        #
        # create DataArray
        dims   = ["wt", "wd", "ws", "z"]
        coords = dict(wt=range(len(self.layout)),
                      wd=self.wd,
                      ws=self.ws,
                      z=hr.index.values,
                      x=(["wt"], self.layout.x.values),
                      y=(["wt"], self.layout.y.values),
                      nm=(["wt"], self.wt),
                     )
        attrs = dict(description=self.name, ptype=ptype)
        res = xr.DataArray(data=r4, dims=dims, coords=coords, attrs=attrs)
        logging.info(f'Result dimension: {dict(res.sizes)}. Profile type: {ptype}')
        #
        return res, list(profs.values())
#
    def performance(self, full=True, plot_it=True, ms=4, ylim=None):
        '''
        get (and optionally plot) performance-data for a given case (i.e. simulation time, residuals etc)
        - input
          * case
          * full   : boolean. also get (and plot?) conservation of variables?
          * plot_it: boolean. plot data?
          * ms     : marker size
          * ylim   : ylim for imbalances (like [-40, 40]). see note1
        - returns
          * perf   : DataFrame with relevant data. column duration is solving-time. also kept as attribute _perf
        - notes
          * note1: imbalances are plotted in % (but data is not %)
        '''
        tmpfile = '/tmp/grepped.txt'
        cmd = f"grep 'CFD Solver wall clock' {self.pm['target directory']}/*.out > {tmpfile}"
        logging.info(cmd)
        os.system(cmd)
        perf = pd.read_csv(tmpfile, sep=' ', usecols=[0,6], header=None, names=['nm', 'duration'], converters=dict(nm=UT.basename))
        vals = [x.split('_')[-3:-1] for x in perf.nm]
        vs = [int(x[0][-1]) for x in vals]
        thetas = [int(x[1]) for x in vals]
        perf['theta'] = thetas
        perf['v'] = vs
        #
        if full:
            eqs = ['U-Mom', 'V-Mom', 'W-Mom', 'P-Mass']
            cmd = f"grep -a7 'Normalised Imbalance Summary' {self.pm['target directory']}/*.out | grep"
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
            for i, v in enumerate(self.ws):
                p = perf[perf.v==i+1]
                plt.plot(p.theta.values, p.duration.values, 's', ms=ms, label=f'v = {v:.1f} m/s')
            plt.xlabel('Wind direction')
            plt.ylabel('Simulation runtime [s]')
            plt.title(f"Case {self.name}")
            plt.legend(loc='best')
            if full:
                for eq in eqs:
                    plt.figure()
                    for i, v in enumerate(self.ws):
                        p = perf[perf.v==i+1]
                        plt.plot(p.theta.values, p[eq].values*100, 's', ms=ms, label=f'v = {v:.1f} m/s')
                    plt.xlabel('Wind direction')
                    plt.ylabel(f'{eq} imbalance [%]')
                    plt.title(f"Case {self.name}")
                    plt.legend(loc='best')
                    if ylim: plt.ylim(*ylim)
        #
        os.unlink(tmpfile)
        self._perf = perf
        return perf
#
    def efficency(self):
        '''
        for a WindModeller-run of type
          'wind data transposition option' = 'Offshore Array Efficiency and Effective TI'
        an efficiency matrix for each turbine is written to file.
        this routine reads this file into a DataArray for ease of access
        - notes
          * in the csv-file, it reports wd's and ws's that are slightly off the requested values.
            here, these values are replaced by the requested values (found in the case-file)
        - returns
          * DataArray. unit is %
        '''
        if not self.pm['wind data transposition option'] == 'Offshore Array Efficiency and Effective TI':
            raise Exception("'wind data transposition option' must be 'Offshore Array Efficiency and Effective TI'")
        fnm = f'{self.pm["target directory"]}/array_eff_TI/results/TurbineEfficiency.csv'
        logging.info(f'reading {fnm}')
        #
        csvfile = open(fnm)
        m3 = []                                     # 3-dim matrix for all data
        nms = []
        vals = np.zeros((len(self.layout), len(self.ws), len(self.wd)))
        k = -1
        for line in csvfile:
            line = line.strip()
            if (not line) or ('Upstream' in line) or ('Turbine efficiency' in line):
                continue
            #
            if 'Turbine:' in line:
                k += 1
                i, j = -1, 0
                nm = line.split('Turbine:,')[1].strip()
                nms.append(nm)
                continue
            #
            # 'default': read data into matrix
            i += 1
            rec = [float(x) for x in line.split(',')]
            vals[k,j,i] = rec[-1]
            if (i+1) % len(self.wd) == 0:
                i = -1
                j += 1
        csvfile.close()
        #
        # build a DataArray
        dims   = ['wt', 'ws', 'wd']
        coords = dict(
                   wt=range(len(nms)),
                   ws=self.ws,
                   wd=self.wd,
                   x =(['wt'], self.layout.x.values),
                   y =(['wt'], self.layout.y.values),
                   nm=(['wt'], nms)
                 )
        attrs = dict(description=self.name, rtype=TYP_EFF, unit=UNIT_PRCNT)
        eff = xr.DataArray(data=vals, dims=dims, coords=coords, attrs=attrs)
        logging.info(f'Result dimension: {dict(eff.sizes)}')
        #
        #  make sure data is consistent
        assert np.all(self.wt == nms)
        #
        return eff
#
    def net(self, reread=False):
        '''
        read simulated power per wt, wd, ws.
        - input:
          * reread: force re-reading if this has been read already
        '''
        if reread or not hasattr(self, '_net'):
            self._net = self._read_results(rtype=TYP_PWR)
        return self._net
#
    def gross(self, reread=False):
        '''
        see docs for efficency()
        this one is just net / efficency
        - input:
          * reread: see net()
        '''
        if not reread and hasattr(self, '_gross'):
            return self._gross
        net = self.net(reread=reread)
        eff = self.efficency()
        scaler = 100. if eff.attrs['unit'] == UNIT_PRCNT else 1.
        gross = net / (eff/scaler)
        gross.attrs = net.attrs.copy()
        self._gross = gross
        self._eff   = eff
        return gross
#
    def wakeloss(self, ws_list=None, wd_list=None, average=True):
        '''
        calculating wakeloss for each turbine
        - input
          * ws_list: selection of wind-speeds. None means all..
          * wd_list: selection of wind-dirs. None means all..
          * average: boolean => averaging over wd and ws
        - output
          * DataArray. unit is %
        '''
        net   = self.net()
        gross = self.gross()
        if not ws_list is None:
            net   = net.sel(ws=ws_list)
            gross = gross.sel(ws=ws_list)
        if not wd_list is None:
            net   = net.sel(wd=wd_list)
            gross = gross.sel(wd=wd_list)
        if average:
            net   = net.mean('ws').mean('wd')
            gross = gross.mean('ws').mean('wd')
        self.wl = (1 - net/gross) * 100
        self.wl.attrs = dict(rtype=TYP_WKL, description=self.name, unit=UNIT_PRCNT)
        return self.wl

# TEST FUNCTIONS

def test_efficency(testcase):
    wmc = WindModellerCase(testcase)
    eff = wmc.efficency()
    assert np.all(eff.values == np.linspace(0.1, 1, 10))

