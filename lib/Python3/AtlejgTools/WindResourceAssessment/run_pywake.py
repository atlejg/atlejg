#!/usr/bin/env python3

'''
this is a demo to show that pywake could be used for the WindFlow project (David H)


there are still a few things to be worked out:
    - it seems weibull-parameters change from location to location. why?
      (i dont pick up this, i only use one set of weilbull-paramters)

NOTES

 - note1
    the yml-input file in the __main__ part should look like this:

        case_nm:                   # if empty, will use basename of *this* file
            Hywind Scotland
        layout_file:
            ../InputData/layout1.csv
        weibull_file:
            ../InputData/weibull1.csv
        wtg_file:
            ../InputData/SWT-6.0-154 DDG3.wtg
        wake_model:
            NOJ
        tp_A:
            !!float   0.60         # Only used for TurboPark. 'A' parameter in dDw/dx = A*I(x)
        noj_k:
            !!float   0.04         # Only used for Jensen. Wake expansion parameter
        turb_intens:
            !!float   0.053
        ws_min:
            !!float    3.
        ws_max:
            !!float   30.
        delta_winddir:
            !!float   1.0          # delta wind-dir for calculations
        delta_windspeed:
            !!float   0.50         # delta wind-vel for calculations
        outfile:
            results.txt

 - note2
    naming conventions follows typical pywake lingo:
        - wt: wind turbine
          wd: wind direction
          ws: wind speed
          wl: wake loss

'''

import py_wake
from py_wake.wind_turbines import WindTurbines
from py_wake.site._site import UniformWeibullSite
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import os, sys, time
import glob, re, zipfile, logging, yaml
from scipy.interpolate import interp1d
import AtlejgTools.WindResourceAssessment.Utils as WU
import AtlejgTools.Utils as UT

N_SECTORS = 12
EPS       = 1e-9               # small non-zero value



def read_layout(fnm, sheetnm=None, n_wtgs=-1):
    if fnm.endswith('.csv'):
        layout = pd.read_csv(fnm, sep=';')
    else:
        layout = pd.read_excel(fnm, sheet_name=sheetnm, header=None, names=['x','y','nm'])
        layout['parknm'] = 'park'
        layout.nm = [str(val) for val in layout.nm]
    if  n_wtgs > 0:
        layout = layout.iloc[:n_wtgs,:]
    return layout

def write_for_webviz(fnm, net, gross):
    '''
    reports at farm level
    '''
    wl    = (gross - net) / gross * 100                      # wakeloss in %
    #maxpwr = max(wtgs.power(ws))
    #cf = net / (len(layout) * maxpwr * 365*24 / 1e9)         # capacity factor based on energy [GWh/GWh]
    f = open(fnm, 'w')
    f.write('ZONE,REGION,GROSS,NET,WAKELOSS,EFFIENCY\n')
    f.write(f'park,park,{gross:.1f},{net:.1f},{wl:.2f},{100-wl:.2f}\n')
    f.close()
    logging.info(f'writing results for webviz: {fnm}')

def write_results(opts, wtgs, net, gross):
    '''
    reports at wtg level
    '''
    net     = net.sum('wd').sum('ws')
    gross   = gross.sum('wd').sum('ws')
    wtg     = os.path.basename(opts.wtg_file)
    weibull = os.path.basename(opts.weibull_file)
    res = pd.DataFrame()
    res['x']      = net.x
    res['y']      = net.y
    res['h']      = net.h
    res['diam']   = wtgs.diameter()
    res['net']    = net.values
    res['gross']  = gross.values
    res['wkmod']  = opts.wake_model
    res['A_scl']  = opts.A_scaler
    res['layout'] = opts.layout_sheet
    res['nWTGs']  = opts.n_wtgs
    res['wtg']    = wtg
    res['weib']   = weibull
    res.to_csv(opts.outfile, sep=',')
    logging.info(f'writing results to : {opts.outfile}')

def set_defaults(opts):
    if not hasattr(opts, 'layout_sheet'): opts.layout_sheet = ''
    if not hasattr(opts, 'n_wtgs'):       opts.n_wtgs       = -1
    if not hasattr(opts, 'A_scaler'):     opts.A_scaler     = 1.

def get_weibull(fnm, A_scaler=1.):
    '''
    - input
      * fnm: either a csv-file or a rasmus's xlsx-format
    '''
    if fnm.endswith('.csv'):
        weib = pd.read_csv(fnm, sep=';')
        weib.columns = ['dir', 'freq', 'A', 'K']
    elif fnm.endswith('.xlsx'):
        weib = pd.read_excel(fnm, skiprows=6, header=None, names=['A', 'K', 'freq'], usecols=[1,2,3], nrows=N_SECTORS)
        weib['dir'] = np.linspace(0, 360, N_SECTORS, endpoint=False)
    else:
        raise Exception(f'weibull-file {fnm} not supported')
    weib.A *= A_scaler
    return weib

def read_wtgs(wtg_file, gap=30.):
    '''
    - input
      * wtg_file: either a wtg-file or a generic Equinor WTG in xlsx-format, typically EQN-D250-15MW.xlsx
      * gap     : minimum distance from blade-tip to sea-surface
    '''
    if wtg_file.endswith('.wtg'):
        return WindTurbines.from_WAsP_wtg(wtg_file)
    elif wtg_file.endswith('.xlsx'):
        data = pd.read_excel(wtg_file, sheet_name='PowerCurve', header=0)
        scaler = 1000. if 'kW' in data.columns[1] else 1.
        data.set_axis(['ws', 'pwr', 'ct'], axis=1, inplace=True)
        nm = UT.basename(wtg_file)
        diam = float(re.search('-D(\d*)-', wtg_file).groups()[0]) # extract it from EQN-D250-15MW.xlsx
        height = diam/2 + gap
        ct_func = interp1d(data.ws, data.ct,  bounds_error=False, fill_value=(EPS,EPS))
        pwr_func = interp1d(data.ws, scaler*data.pwr,  bounds_error=False, fill_value=(EPS,EPS))
        return WindTurbines(names=[nm], diameters=[diam], hub_heights=[height], ct_funcs=[ct_func], power_funcs=[pwr_func], power_unit='W')
    else:
        raise Exception(f'wtg-file {wtg_file} not supported')


def main(yml_file):
    '''
    - input
      * yml_file
    - returns
      * sim, layout, opts, wtgs, site, wake_model
    '''
    #
    # house-keeping
    tic = time.perf_counter()
    logging.basicConfig(level=logging.INFO)
    #
    opts = UT.get_yaml(yml_file)
    set_defaults(opts)
    logging.info(f"Reading input file: {yml_file}")
    #
    # setup things
    layout = read_layout(opts.layout_file, opts.layout_sheet, opts.n_wtgs)
    wtgs   = read_wtgs(opts.wtg_file)
    weib   = get_weibull(opts.weibull_file)
    site   = UniformWeibullSite(weib.freq, weib.A, weib.K, opts.turb_intens)
    #
    # initialize the chosen wake model
    wake_model = opts.wake_model.upper()
    if opts.wake_model.upper() == 'TP' or 'TURBO' in opts.wake_model.upper():
        wake_model = py_wake.TP(site, wtgs, k=opts.tp_A)
    elif opts.wake_model.upper()== 'NOJ':
        wake_model = py_wake.NOJ(site, wtgs, k=opts.noj_k)
    else:
        raise Exception('TurboPark or the NOJ wake models are the only options available.')
    #
    # run wake model for all combinations of wd and ws
    wd = np.arange(0, 360, opts.delta_winddir)
    ws = np.arange(np.floor(opts.ws_min), np.ceil(opts.ws_max), opts.delta_windspeed)
    sim = wake_model(layout.x.values, layout.y.values, wd=wd, ws=ws)
    #
    # calculate and report
    net   = sim.aep(with_wake_loss=True)
    gross = sim.aep(with_wake_loss=False)
    write_results(opts, wtgs, net, gross)
    #
    if opts.webviz_output:
        write_for_webviz(opts.webviz_file, net.sum().values, gross.sum().values)
    #
    # house-keeping
    toc = time.perf_counter()
    logging.info(f"Total runtime: {toc-tic:0.1f} seconds")
    #
    return sim, layout, opts, wtgs, weib, site, wake_model

################################## -- MAIN LOGIC -- ###########################


if __name__ == '__main__':

    # run program
    yml_file = sys.argv[1]
    sim, layout, opts, wtgs, weib, site, wake_model = main(yml_file)

