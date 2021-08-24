#!/usr/bin/env python3

'''
this is a demo to show that pywake could be used for the WindFlow project (David H)


there are still a few things to be worked out:
    - it seems weibull-parameters change from location to location. why?
      (i dont pick up this, i only use one set of weilbull-paramters)

NOTES

 - note1
    assumes only one WTG-type.

 - note2
    assumes only one Weibull-distribution

 - note3
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
import os, sys, time
import re, logging
from scipy.interpolate import interp1d
import AtlejgTools.Utils as UT

N_SECTORS = 12
EPS       = 1e-9               # small non-zero value

def read_layout(fnm, sheetnm=None, n_wtgs=-1):
    '''
    supported:
      *.csv  : must have columns x, y, nm
      *.xlsx : must have columns x, y, nm
      *.txt  : format: WindModeller layout
    '''
    if fnm.endswith('.csv'):
        layout = pd.read_csv(fnm, sep=';')
    elif fnm.endswith('.xlsx'):
        layout = pd.read_excel(fnm, sheet_name=sheetnm, header=None, names=['x','y','nm'])
        layout['parknm'] = 'park'
        layout.nm = [str(val) for val in layout.nm]
    elif fnm.endswith('.txt'):
        layout = pd.read_csv(fnm, header=None, delim_whitespace=True, usecols=[0,1,2], names=['nm', 'x','y'])
        layout['parknm'] = 'park'
    else:
        raise Exception(f'File-type of file {fnm} is not supported')
    if  n_wtgs > 0:
        layout = layout.iloc[:n_wtgs,:]
    return layout

def write_for_webviz(fnm, net, gross):
    '''
    reports at farm level
    - input:
      * fnm    : file to write
      * net    : net AEP for the farm
      * gross  : gross AEP for the farm
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
    - input:
      * opts   : options from the yaml
      * wtgs   : WindTurbines 
      * net    : net AEP per wtg
      * gross  : gross AEP per wtg
    '''
    wtg     = os.path.basename(opts.wtg_file)
    weibull = os.path.basename(opts.weibull_file)
    #
    res = pd.DataFrame()
    res['wtg_id']      = net.wt
    res['x']           = net.x
    res['y']           = net.y
    res['h']           = net.h
    res['diam']        = wtgs.diameter()
    res['net_AEP']     = net.values
    res['gross_AEP']   = gross.values
    res['wakemodel']   = opts.wake_model
    res['A_scaler']    = opts.A_scaler
    res['layout']      = opts.layout_sheet
    res['nWTGs']       = opts.n_wtgs
    res['wtg']         = wtg
    res['weibull']     = weibull
    res['turb_intens'] = opts.turb_intens
    res['noj_k']       = opts.noj_k
    res['tp_A']        = opts.tp_A
    #
    res.to_csv(opts.outfile, sep=',')
    #
    logging.info(f'writing results to : {opts.outfile}')

def set_defaults(opts):
    if not hasattr(opts, 'layout_sheet'):    opts.layout_sheet    = ''
    if not hasattr(opts, 'n_wtgs'):          opts.n_wtgs          = -1
    if not hasattr(opts, 'A_scaler'):        opts.A_scaler        = 1.
    if not hasattr(opts, 'case_nm'):         opts.case_nm         = 'pywake'
    if not hasattr(opts, 'wake_model'):      opts.wake_model      = 'TurbOPark'
    if not hasattr(opts, 'tp_A'):            opts.tp_A            = 0.60
    if not hasattr(opts, 'noj_k'):           opts.noj_k           = 0.04
    if not hasattr(opts, 'turb_intens'):     opts.turb_intens     = 0.05
    if not hasattr(opts, 'ws_min'):          opts.ws_min          = 2.
    if not hasattr(opts, 'ws_max'):          opts.ws_max          = 30.
    if not hasattr(opts, 'delta_winddir'):   opts.delta_winddir   = 1.
    if not hasattr(opts, 'delta_windspeed'): opts.delta_windspeed = 0.5
    if not hasattr(opts, 'outfile'):         opts.outfile         = 'pywake_results.csv'
    if not hasattr(opts, 'webviz_output'):   opts.webviz_output   = False
    if not hasattr(opts, 'webviz_file'):     opts.webviz_file     = ''
    if not hasattr(opts, 'wtg_file'):        opts.wtg_file        = ''
    if not hasattr(opts, 'pwr_file'):        opts.pwr_file        = ''
    if not hasattr(opts, 'ct_file'):         opts.ct_file         = ''
    if not hasattr(opts, 'hub_height'):      opts.hub_height      = -1.

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

def read_wtgs(opts):
    '''
    - input
      * opts: either wtg_file or (pwr_file and ct_file) must be given. 
          - wtg_file:            either a wtg-file or a generic Equinor WTG in xlsx-format, typically EQN-D250-15MW.xlsx
          - pwr_file & ct_file:  WindModeller pwr & ct curves
    '''
    if opts.wtg_file:
        if opts.wtg_file.endswith('.wtg'):
            return WindTurbines.from_WAsP_wtg(opts.wtg_file)
        elif opts.wtg_file.endswith('.xlsx'):
            data = pd.read_excel(opts.wtg_file, sheet_name='PowerCurve', header=0)
            scaler = 1000. if 'kW' in data.columns[1] else 1.
            data.set_axis(['ws', 'pwr', 'ct'], axis=1, inplace=True)
            nm = UT.basename(opts.wtg_file)
            diam = float(re.search('-D(\d*)-', opts.wtg_file).groups()[0]) # extract it from EQN-D250-15MW.xlsx
            ct_func = interp1d(data.ws, data.ct,  bounds_error=False, fill_value=(EPS,EPS))
            pwr_func = interp1d(data.ws, scaler*data.pwr,  bounds_error=False, fill_value=(EPS,EPS))
            return WindTurbines(names=[nm], diameters=[diam], hub_heights=[opts.hub_height], ct_funcs=[ct_func], power_funcs=[pwr_func], power_unit='W')
        else:
            raise Exception(f'wtg-file {opts.wtg_file} not supported')
    else:
        pwr_func, _, attr = wm.read_curve(opts.pwr_file, True)
        ct_func           = wm.read_curve(opts.ct_file, True)[0]
        unit = re.search('\[(\w*)\]', un[-1]).group()[1:-1]
        return WindTurbines(names=['wtg'], diameters=[opts.diam], hub_heights=[opts.hub_height], ct_funcs=[ct_func], power_funcs=[pwr_func], power_unit=unit)




def main(yaml_file):
    '''
    - input
      * yaml_file like this:
            #
            case_nm:            !!str      Doggerbank-C Demo
            #
            layout_file:        !!str      ../../../../InputData/Layouts/<LAYOUTFILE>
            layout_sheet:       !!str      <LAYOUT>
            weibull_file:       !!str      ../../../../InputData/WindMaps/<WINDRESOURCE>
            A_scaler:           !!float    <A_SCALER>
            wtg_file:           !!str      ../../../../InputData/WTGs/<WTG>
            n_wtgs:             !!int      <N_WTGS>
            #
            wake_model:         !!str      <WAKE_MODEL>
            #
            tp_A:               !!float    0.60                                      # Only used for TurbOPark. 'A' parameter in dDw/dx = A*I(x)
            noj_k:              !!float    0.04                                      # Only used for Jensen. Wake expansion parameter
            turb_intens:        !!float    0.053
            #
            ws_min:             !!float    2.
            ws_max:             !!float    30.
            delta_winddir:      !!float    1.0                                       # delta wind-dir for calculations
            delta_windspeed:    !!float    0.50                                      # delta wind-vel for calculations
            #
            outfile:            !!str      pywake_results.csv
            webviz_output:      !!bool     true
            webviz_file:        !!str      share/results/volumes/webviz.csv
    - returns
      * sim, layout, opts, wtgs, weib, site, wake_model
    - notes
      * many attributes of the yaml-file has a default value; see set_defaults()
    '''
    #
    # house-keeping
    tic = time.perf_counter()
    logging.basicConfig(level=logging.INFO)
    #
    opts = UT.get_yaml(yaml_file)
    set_defaults(opts)
    logging.info(f"Reading input file: {yaml_file}")
    #
    # setup things
    layout = read_layout(opts.layout_file, opts.layout_sheet, opts.n_wtgs)
    wtgs   = read_wtgs(opts)
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
    write_results(opts, wtgs, net.sum('wd').sum('ws'), gross.sum('wd').sum('ws'))
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
    yaml_file = sys.argv[1]
    sim, layout, opts, wtgs, weib, site, wake_model = main(yaml_file)

