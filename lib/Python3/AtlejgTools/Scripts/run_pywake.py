#!/usr/bin/env python3

'''

Oct 2021

improved version of run_pywake.py
this one handles multiple parks & multiple wtg's
eventually, it should be merged with pywake_for_knowl.py


NOTES

 - note1
    assumes only one Weibull-distribution

 - note2
    naming conventions follows typical pywake lingo:
        - wt: wind turbine
          wd: wind direction
          ws: wind speed
          wl: wake loss

 - note3
   must have an old version of xlrd: pip install xlrd==1.2.0

 - note4
        only TP/TurbOPark, ETP (Equinor TP), ZGauss, NOJ, NOJLocal and Fuga
        wake models are available.
        for Fuga, only 1 wtg-type is allowed

'''

import py_wake
from py_wake.wind_turbines import WindTurbines
from py_wake.site._site import UniformWeibullSite
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, time
import re, logging, argparse
from scipy.interpolate import interp1d
import AtlejgTools.Utils as UT
import AtlejgTools.WindYieldAssessment.Utils as WU
import AtlejgTools.WindYieldAssessment.WindModeller as wm
import AtlejgTools.Scripts.pywake_for_knowl as pfk

N_SECTORS = 12

def read_setup(opts):
    '''
    - input
      * opts must have the following attributes
          - layout_files  : csv-file which must have columns x, y, nm, parknm, wtg
          - wtg_directory : where to find the wtg-files
          - default_wtg   : wtg to use if wtg-field is empty in the layout-file
          - default_parknm: parknm to use if parknm-field is empty in the layout-file
          - sep           : separator for csv-file
    - returns
      * case: a case-object similair to the one used in pywake_for_knowl.
              has an attribute 'ids' with unique id for each location. id follows
              the same pattern as knowl uses (i.e. 1001, 1002, ... 1054, 2001, 2002,...)
      * wtgs: WindTurbines object
    '''
    layout = pd.DataFrame([])
    for fnm in opts.layout_files:
        layout0 = pd.read_csv(fnm, sep=opts.csv_sep)
        # put in default values for empty ones
        if not 'wtg'    in layout0.columns: layout0['wtg']    = default_wtg
        if not 'parknm' in layout0.columns: layout0['parknm'] = default_parknm
        layout0.fillna({'wtg':opts.default_wtg, 'parknm':opts.default_parknm}, inplace=True)
        layout = pd.concat([layout,layout0])
    #
    case = UT.Struct()
    case._wtgs = layout.wtg
    #
    # we group parks according to park-name (not wtg-type)
    # and have 1 wtg per park. this means 'types' will also
    # follow park-name
    case.park_nms = list(layout.parknm.unique())
    #
    case.xs        = layout.x.values
    case.ys        = layout.y.values
    case.nms       = layout.nm.values
    case.types     = np.zeros(len(layout), dtype=int)
    case.ids       = np.zeros(len(layout), dtype=int)
    case._wtg_list = []
    case.park_list = []
    for i, parknm in enumerate(case.park_nms):
        ixs = (layout.parknm==parknm).values.nonzero()[0]
        #
        case.types[ixs] = i
        #
        n_wtgs = len(ixs)
        case.ids[ixs] = (i+1)*1000 + np.arange(n_wtgs) + 1
        #
        wtg = layout.iloc[ixs[0]].wtg
        case._wtg_list.append(wtg)
        #
        logging.info(f'{parknm} :: {wtg}')
        #
        # WU assumes the park_list has useful info
        park = UT.Struct()
        park.name = parknm
        park.size = n_wtgs
        park.ids  = case.ids[ixs]
        park.xs   = case.xs[ixs] 
        park.ys   = case.ys[ixs] 
        case.park_list.append(park)
    #
    wtgs = WindTurbines.from_WAsP_wtg([f'{opts.wtg_directory}/{wtg}' for wtg in case._wtg_list])
    #
    for i, park in enumerate(case.park_list):
        park.wtg = UT.Struct()
        park.wtg.hub_height = wtgs.hub_height(i)
        park.wtg.name = wtgs.name(i)
    #
    return case, wtgs

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

def write_results(case, opts, wtgs, net, gross, suffix='.csv'):
    '''
    reports at wtg level
    - input:
      * opts   : options from the yaml
      * wtgs   : WindTurbines 
      * net    : net AEP per wtg
      * gross  : gross AEP per wtg
    '''
    weibull = os.path.basename(opts.weibull_file)
    #
    res = pd.DataFrame()
    res['wtg_id']      = net.wt
    res['x']           = net.x
    res['y']           = net.y
    res['h']           = net.h
    res['diam']        = wtgs.diameter(case.types)
    res['net_AEP']     = net.values
    res['gross_AEP']   = gross.values
    res['wakemodel']   = opts.wake_model
    res['A_scaler']    = opts.A_scaler
    res['weibull']     = weibull
    res['turb_intens'] = opts.turb_intens
    res['noj_k']       = opts.noj_k
    res['tp_A']        = opts.tp_A
    res['wtg']         = wtgs.name(case.types)
    #
    fnm = opts.resultfile_prefix + suffix
    res.to_csv(fnm, sep=opts.csv_sep)
    #
    logging.info(f'writing results to : {fnm}')

def set_defaults(opts):
    if not hasattr(opts, 'A_scaler'):          opts.A_scaler          = 1.
    if not hasattr(opts, 'case_nm'):           opts.case_nm           = 'pywake'
    if not hasattr(opts, 'wake_model'):        opts.wake_model        = 'ETP'
    if not hasattr(opts, 'tp_A'):              opts.tp_A              = 0.60
    if not hasattr(opts, 'noj_k'):             opts.noj_k             = 0.04
    if not hasattr(opts, 'turb_intens'):       opts.turb_intens       = 0.05
    if not hasattr(opts, 'ws_min'):            opts.ws_min            = 2.
    if not hasattr(opts, 'ws_max'):            opts.ws_max            = 30.
    if not hasattr(opts, 'delta_winddir'):     opts.delta_winddir     = 1.
    if not hasattr(opts, 'delta_windspeed'):   opts.delta_windspeed   = 0.5
    if not hasattr(opts, 'webviz_output'):     opts.webviz_output     = False
    if not hasattr(opts, 'webviz_file'):       opts.webviz_file       = ''
    if not hasattr(opts, 'wtg_directory'):     opts.wtg_directory     = '.'
    if not hasattr(opts, 'default_wtg'):       opts.default_wtg       = ''
    if not hasattr(opts, 'weibull_file'):      opts.weibull_file      = ''
    if not hasattr(opts, 'csv_sep'):           opts.csv_sep           = ';'
    if not hasattr(opts, 'dump_results'):      opts.dump_results      = True
    if not hasattr(opts, 'resultfile_prefix'): opts.resultfile_prefix = 'pywake_results'

def get_weibull(fnm, A_scaler=1., is_percent=True):
    '''
    - input
      * fnm: either a csv-file or a rasmus's xlsx-format (which is almost the same as knowl uses)
    '''
    if fnm.endswith('.csv'):
        weib = pd.read_csv(fnm, sep=';')
        weib.columns = ['dirs', 'freqs', 'As', 'Ks']
    elif fnm.endswith('.xlsx'):
        weib = pd.read_excel(fnm, skiprows=6, header=None, names=['As', 'Ks', 'freqs'], usecols=[1,2,3], nrows=N_SECTORS)
        weib['dirs'] = np.linspace(0, 360, N_SECTORS, endpoint=False)
    else:
        raise Exception(f'weibull-file {fnm} not supported')
    weib.As *= A_scaler
    if is_percent: weib.freqs /= 100.
    return weib

def deficit_model(site, wtgs, opts):
    '''
    for the choice of wake_model, see note4
    '''
    wake_model = opts.wake_model.upper()
    if wake_model   == 'FUGA':
        assert wtgs.uniq_wtgs == 1
        wf_model = py_wake.Fuga(opts.lut_path, site, wtgs)
    elif wake_model in ['TP', 'TURBOPARK']:
        wf_model = py_wake.TP(site, wtgs, k=opts.tp_A)       # 'standard' TurbOPark
    elif wake_model == 'ETP':
        wf_model = py_wake.ETP(site, wtgs, k=opts.tp_A)      # 'Equinor' TurbOPark
    elif wake_model == 'NOJ':
        wf_model = py_wake.NOJ(site, wtgs, k=opts.noj_k)
    elif wake_model == 'NOJLOCAL':
        wf_model = py_wake.NOJLocal(site, wtgs)
    elif wake_model == 'ZGAUSS':
        wf_model = py_wake.deficit_models.gaussian.ZongGaussian(
                        site,
                        wtgs,
                        superpositionModel=py_wake.superposition_models.WeightedSum(),
                        turbulenceModel=py_wake.turbulence_models.crespo.CrespoHernandez()
                    )
    else:
        raise Exception('Only Fuga, TP/TurbOPark, ETP (Equinor TP), ZGauss, NOJ, NOJLocal wake models are available.')
    return wf_model

def write_AEP(case, sim, wtgs, weib, delta_winddir, outfile):
    # coarsen it by averaging
    n_sectors = len(weib.freqs)
    width  = 360. / n_sectors
    n_bins    = int(width / delta_winddir)         # number of bins per sector
    simc = sim.coarsen(wd=n_bins).mean()
    #
    aeps = WU.calc_AEP(simc, wtgs, weib, park_nms=case.park_nms, verbose=True)
    WU.create_output(aeps, case, outfile)
    return simc, aeps

def main(mode, yaml_file):
    '''
    - input
      * yaml_file like this:
            case_nm:            !!str      Arkona
            #
            layout_files:
                - ../../../../InputData/wikinger.csv
                - ../../../../InputData/arkona.csv
            wtg_directory:      !!str      ../../../../InputData
            default_wtg:        !!str      SWT-6.0-154 DDG3-r1.230.wtg
            default_parknm:     !!str      arkona
            #
            weibull_file:       !!str      ../../../../InputData/wind_resource.xlsx
            A_scaler:           !!float    1
            #
            wake_model:         !!str      ETP
            tp_A:               !!float    0.6           # Only used for TurbOPark. 'A' parameter in dDw/dx = A*I(x)
            noj_k:              !!float    0.04          # Only used for Jensen. Wake expansion parameter
            turb_intens:        !!float    0.05  
            #
            ws_min:             !!float    2.
            ws_max:             !!float    30.
            delta_winddir:      !!float    1.0              # delta wind-dir for calculations
            delta_windspeed:    !!float    0.50             # delta wind-vel for calculations
            #
            resultfile_prefix:  !!str      pywake_results
            webviz_output:      !!bool     false
            webviz_file:        !!str      none
            csv_sep:            !!str      ;
    - returns
      * sim, simc, aeps, case, opts, wtgs, site, wf_model, weib
        simc is a coarsened version of sim
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
    if mode == 'knowl':
        knowl = pfk.read_knowl_file(opts.knowl_file)
        case, wtgs = pfk.read_inventory(opts.inventory_file)
        weib = knowl.weibulls[opts.weibull_index-1]
        opts.turb_intens = knowl.turb_intens
    else:
        case, wtgs = read_setup(opts)
        weib = get_weibull(opts.weibull_file)
    #
    # initialize
    site = UniformWeibullSite(weib.freqs, weib.As, weib.Ks, opts.turb_intens)
    wf_model = deficit_model(site, wtgs, opts)
    #
    # run simulation
    logging.info(f'run wake model {opts.wake_model} for all combinations of wd and ws')
    wd = np.arange(0, 360, opts.delta_winddir)
    ws = np.arange(np.floor(opts.ws_min), np.ceil(opts.ws_max), opts.delta_windspeed)
    sim = wf_model(case.xs, case.ys, type=case.types, wd=wd, ws=ws)
    #
    # reporting AEP per sector like knowl/fuga
    fnm = opts.resultfile_prefix + '.res'
    simc, aeps = write_AEP(case, sim, wtgs, weib, opts.delta_winddir, fnm)
    #
    # report all data
    net   = sim.aep(with_wake_loss=True)
    gross = sim.aep(with_wake_loss=False)
    write_results(case, opts, wtgs, net.sum('wd').sum('ws'), gross.sum('wd').sum('ws'))
    #
    # house-keeping
    toc = time.perf_counter()
    logging.info(f"Total runtime: {toc-tic:0.1f} seconds")
    #
    vals =  sim, simc, aeps, case, opts, wtgs, site, wf_model, weib
    keys = ['sim', 'simc', 'aeps', 'case', 'opts', 'wtgs', 'site', 'wf_model', 'weib']
    res = dict(zip(keys, vals))
    if opts.dump_results: pfk.dump(res, f'{opts.resultfile_prefix}.pck')
    return res

################################## -- MAIN LOGIC -- ###########################


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(dest='yaml_file', help='yaml-file describing the case')
    parser.add_argument('--knowl', '-k', action='store_true', help='input from Knowl?')
    args = parser.parse_args()

    mode = 'knowl' if args.knowl else 'default'

    res = main(mode, args.yaml_file)

