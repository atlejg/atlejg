#!/usr/bin/env python3

'''
@AUTHOR Atle J. Gyllensten, agy@equinor.com


A pywake-wrapper to be used for running simulations in an effective way,
especially as a part of the WindWorks-initiativ.

Oct 2021
  - now handles multiple parks & multiple wtg's
  - --knowl option to replace pywake_for_knowl.py (which is now deprecated)
  - MERGING run_pywake.py into windworks.py (a script for running ert-like on windows)

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
    must have an old version of xlrd: pip install xlrd==1.2.0 for pd.read_excel to work.

 - note4
    only TP/TurbOPark, ETP (Equinor TP), ZGauss, NOJ, NOJLocal and Fuga
    wake models are available.
    for Fuga, only 1 wtg-type is allowed

 - note5
    must have the __name__ == "__main__" when using apply_async on windows

'''

import py_wake
from py_wake.wind_turbines import WindTurbines
from py_wake.site._site import UniformWeibullSite
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, time, os
import re, logging, argparse
from scipy.interpolate import interp1d
import AtlejgTools.Utils as UT
import AtlejgTools.WindYieldAssessment.Utils as WU
import AtlejgTools.WindYieldAssessment.KnowlUtils as KU
import AtlejgTools.WindYieldAssessment.WindModeller as wm
import multiprocessing as mp

N_SECTORS    = 12
TYPE_EXCEL   = 'excel'
TYPE_WTGFILE = 'wtg'
TEST_REPO    = '/project/RCP/active/wind_yield_assessment/agy/WindWorks/Testdata/'

def wtg_from_excel(fnm):
    '''
    reads excel file into a WindTurbines object
    - input
      * fnm: a generic Equinor WTG in xlsx-format, typically EQN-D250-15MW.xlsx
    - returns
      * a WindTurbines object
    '''
    data = pd.read_excel(fnm, sheet_name='PowerCurve', header=0)
    props = pd.read_excel(fnm, sheet_name='Properties', header=None)
    props = dict(zip(props[0].values, props[1].values))     # make a dict out of it
    scaler = 1000. if 'kW' in data.columns[1] else 1.
    data.set_axis(['ws', 'pwr', 'ct'], axis=1, inplace=True)
    nm = UT.basename(fnm)
    diam = props['Rotor Diameter']
    hub_height = diam/2. + props['Air Gap']
    ct_func = WU.interpolate_curve(data.ws, data.ct)
    pwr_func = WU.interpolate_curve(data.ws, scaler*data.pwr)
    return WindTurbines(names=[nm], diameters=[diam], hub_heights=[hub_height], ct_funcs=[ct_func], power_funcs=[pwr_func], power_unit='W')

def wtg_from_winmodeller(pwr_file, ct_file, diam, hub_height=-1, gap=26.):
    '''
    reads WindModeller input-files into a WindTurbines object
    - input
      * pwr_file: WindModeller profile file for power
      * ct_file: WindModeller profile file for CT
    - returns
      * a WindTurbines object
    '''
    pwr_func, _, attr = wm.read_curve(pwr_file, True)
    ct_func, _, _     = wm.read_curve(ct_file, True)
    unit = re.search('\[(\w*)\]', attr[-1]).group()[1:-1]
    if hub_height < 0:
        hub_height = diam/2. + gap
    return WindTurbines(names=['wtg'], diameters=[diam], hub_heights=[hub_height], ct_funcs=[ct_func], power_funcs=[pwr_func], power_unit=unit)

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
        if fnm.endswith('.csv'):
            layout0 = pd.read_csv(fnm, sep=opts.csv_sep)
        else:
            layout0 = pd.read_excel(fnm)
        # put in default values for empty ones
        if not 'nm'     in layout0.columns: layout0['nm']     = [str(i) for i in range(len(layout0))]
        if not 'wtg'    in layout0.columns: layout0['wtg']    = opts.default_wtg
        if not 'parknm' in layout0.columns: layout0['parknm'] = opts.default_parknm
        layout0.fillna({'wtg':opts.default_wtg, 'parknm':opts.default_parknm}, inplace=True)
        #
        # concat into one big structure
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
    wtg_files = [f'{opts.wtg_directory}/{wtg}' for wtg in case._wtg_list]
    if opts.wtg_filetype == TYPE_EXCEL:
        assert len(wtg_files) == 1, 'For excel-wtgs only 1 park is implemented'
        wtgs = wtg_from_excel(wtg_files[0])
    else:
        wtgs = WindTurbines.from_WAsP_wtg(wtg_files)
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

def test_run_pywake(testdir, yaml_file, orig_prefix, knowl_mode=False, wake_model=None):
    '''
    generic test-function
    - input
      * testdir:     where testdata is found, relative to TEST_REPO
      * yaml_file:   the input file
      * orig_prefix: prefix to csv-file and resultf-file (txt-file) that
                     we will compare new results to
      * knowl_mode:  boolean
      * wake_model:  to be used when running in knowl_mode
    - returns
      * a results-dict
    '''
    cwd = os.getcwd()
    os.chdir(TEST_REPO+testdir)
    logging.info(os.getcwd())
    #
    # run simulatiojn
    res =  main(yaml_file, knowl_mode=knowl_mode, wake_model=wake_model)
    opts = res['opts']
    #
    # compare csv-files
    orig = pd.read_csv(orig_prefix+'.csv', sep=opts.csv_sep)
    new  = pd.read_csv(opts.resultfile_prefix+'.csv', sep=opts.csv_sep)
    #
    logging.info(f'comparing new results to {orig_prefix}.csv')
    txt_cols = ['wakemodel', 'weibull', 'wtg']
    assert np.allclose(new.drop(txt_cols, axis=1), orig.drop(txt_cols, axis=1))
    #
    # compare result-files (txt-files)
    orig = WU.read_output_file(orig_prefix+'.txt')
    new  = WU.read_output_file(opts.resultfile_prefix+'.txt')
    #
    logging.info(f'comparing new results to {orig_prefix}.txt')
    assert np.allclose(new.net, orig.net)
    assert np.allclose(new.gross, orig.gross)
    assert np.allclose(new.wl, orig.wl)
    #
    logging.info(f' testing OK')
    os.chdir(cwd)
    return res

def test_DBC():
    '''
    test on a small Doggerbank-C case
    '''
    logging.info('\n* Testing Doggerbank-C case')
    return test_run_pywake('DBC', '1.yaml', 'orig')

def test_Arkona():
    '''
    test on a Arkona-case 
    '''
    logging.info('\n* Testing Arkona case')
    return test_run_pywake('Arkona/Test_1', '2.yaml', 'orig')

def test_Arkona_knowlinput():
    '''
    test on a Arkona-case with Knowl inputfiles
    '''
    logging.info('\n* Testing Arkona case (knowl-mode)')
    return test_run_pywake('Arkona/Test_Knowl/ARK_ext', '1.yaml', 'orig', knowl_mode=True)

def test_all():
    test_DBC()
    test_Arkona()
    test_Arkona_knowlinput()
    logging.info('\n* ALL TESTS OK')

def set_defaults(opts, knowl_mode):
    if not hasattr(opts, 'A_scaler'):          opts.A_scaler          = 1.
    if not hasattr(opts, 'case_nm'):           opts.case_nm           = 'pywake'
    if not hasattr(opts, 'wake_model'):        opts.wake_model        = 'ETP'
    if not hasattr(opts, 'tp_A'):              opts.tp_A              = 0.60
    if not hasattr(opts, 'noj_k'):             opts.noj_k             = 0.04
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
    if not hasattr(opts, 'dump_results'):      opts.dump_results      = False
    if knowl_mode:
        if not hasattr(opts, 'inventory_file'):    opts.inventory_file    = 'Inventory.xml'
        if not hasattr(opts, 'weibull_index'):     opts.weibull_index     = 1
        if not hasattr(opts, 'resultfile_prefix'): opts.resultfile_prefix = 'FugaOutput_1'
        if not hasattr(opts, 'selected'):          opts.selected          = []       # [] means all
        if not hasattr(opts, 'knowl_file'):
            opts.knowl_file = glob.glob('knowl_v*input.xlsx')[0]
    else:
        if not hasattr(opts, 'resultfile_prefix'): opts.resultfile_prefix = 'pywake_results'
        if not hasattr(opts, 'turb_intens'):       opts.turb_intens       = 0.05
        if not hasattr(opts, 'wtg_filetype'):      opts.wtg_filetype      = TYPE_WTGFILE

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

def write_AEP(case, pwr, wtgs, weib, outfile):
    '''
    write AEP
    '''
    n_sectors = len(weib.freqs)
    pwrc = WU.coarsen(pwr, n_sectors)
    #
    aeps = WU.calc_AEP(pwrc, wtgs, weib, park_nms=case.park_nms, verbose=True)
    WU.create_output(aeps, case, outfile)
    return pwrc, aeps

def write_AEP_old(case, pwr, wtgs, weib, outfile):
    '''
    write AEP
    '''
    n_sectors = len(weib.freqs)
    dwd = pwr.wd.values[1] - pwr.wd.values[0]
    width  = 360. / n_sectors
    n_bins    = int(width / dwd)         # number of bins per sector
    pwrc = pwr.coarsen(wd=n_bins).mean()
    #
    aeps = WU.calc_AEP(pwrc, wtgs, weib, park_nms=case.park_nms, verbose=True)
    WU.create_output(aeps, case, outfile)
    return pwrc, aeps

def main(yaml_file, wake_model=None, knowl_mode=False):
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
      * wake_model: 'ETP' etc (see note4) or None.
        if set, it will override the value from opts (yaml-file).
        this option is needed when being called from Knowl.
      * knowl_mode : True if you are using knowl-input files
    - returns
      * sim, pwrc, aeps, case, opts, wtgs, site, wf_model, weib
        pwrc is a wd-coarsened version of sim.Power. it is coarsened from delta_winddir to the sectors
        given by the weibull.
    - notes
      * many attributes of the yaml-file has a default value; see set_defaults()
    '''
    #
    # house-keeping
    tic = time.perf_counter()
    logging.basicConfig(level=logging.INFO)
    #
    if yaml_file:
        logging.info(f"Reading input file: {yaml_file}")
        opts = UT.get_yaml(yaml_file)
        if not hasattr(opts, 'resultfile_prefix'): opts.resultfile_prefix = UT.basename(yaml_file)
    else:
        opts = UT.Struct()
    set_defaults(opts, knowl_mode)
    if wake_model: opts.wake_model = wake_model
    #
    # setup things
    if not knowl_mode:
        case, wtgs = read_setup(opts)
        weib = get_weibull(opts.weibull_file)
    else:
        knowl = KU.read_knowl_file(opts.knowl_file)
        case, wtgs = KU.read_inventory(opts.inventory_file, selected=opts.selected)
        weib = knowl.weibulls[opts.weibull_index-1]
        opts.turb_intens = knowl.turb_intens
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
    fnm = opts.resultfile_prefix + '.txt'
    pwrc, aeps = write_AEP(case, sim.Power, wtgs, weib, fnm)
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
    vals =  sim, pwrc, aeps, case, opts, wtgs, site, wf_model, weib
    keys = ['sim', 'pwrc', 'aeps', 'case', 'opts', 'wtgs', 'site', 'wf_model', 'weib']
    res = dict(zip(keys, vals))
    if opts.dump_results: KU.dump(res, f'{opts.resultfile_prefix}.pck')
    return res

################################## -- MAIN LOGIC -- ###########################

def parse_ertconfig(ertfile, return_raw=False):
    '''
    reads ert config-file into Struct
    not too smart... does not substitute define'd values etc.
    tries to hide the ugly stuff
    '''
    raw = {}
    for line in open(ertfile):
        line = line.strip()
        if not line             : continue                                         # empty lines
        if line.startswith('--'): continue                                         # comments
        rec = line.split()
        key = rec[0]
        if key in ('DEFINE', 'SIMULATION_JOB', 'FORWARD_MODEL') and not key in raw:  # put multiple entries in list
            raw[key] = {}
        if key in raw: raw[key][rec[1]] = rec[2:]
        else         : raw[key]         = rec[1:]
    #
    # the UGLY stuff
    cfg = UT.Struct()
    cfg.design_file    = list(raw['FORWARD_MODEL'].keys())[0].split('/')[-1].replace(',', '').strip()
    cfg.templ_file     = '/'.join(list(raw['FORWARD_MODEL'].keys())[1].split('/')[-2:]).replace(',', '').strip()
    cfg.result_file    = list(raw['FORWARD_MODEL'].values())[1][0].split('/')[-1].replace(')', '').strip()
    cfg.runpath        = raw['RUNPATH'][0]
    cfg.n_realizations = int(raw['NUM_REALIZATIONS'][0])
    cfg.max_submit     = int(raw['MAX_SUBMIT'][0])
    #
    loglevel    = raw['LOG_LEVEL'][0]
    cfg.loglevel = logging.DEBUG if loglevel=='DEBUG' else logging.INFO
    #
    if return_raw: return cfg, raw
    return cfg

def run_simulation(simfile):
    cwd = os.getcwd()
    simdir, yaml_file = os.path.split(simfile) 
    os.chdir(simdir)
    logging.info(f'running {yaml_file} in {simdir}')
    #
    main(yaml_file)
    #
    os.chdir(cwd)
    return

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument('--ert', '-e', help='ert config-file (for running multiple cases)', default=None)
    parser.add_argument('--yaml', '-y', help='yaml-file describing the case', default=None)
    parser.add_argument('--wakemodel', '-w', help='overwrites input in yaml-file', default=None) # useful for knowl
    parser.add_argument('--knowl', '-k', action='store_true', help='input is from Knowl?')
    args = parser.parse_args()
    
    if args.ert:
        cfg = parse_ertconfig(args.ert)
        #
        logging.basicConfig(level=cfg.loglevel)
        #
        dm = pd.read_excel(cfg.design_file)    # design-matrix
        #
        keys = [f'<{k}>' for k in dm.columns[3:].values]
        #
        n_realizations = min(cfg.n_realizations, len(dm))
        logging.info(f'running {n_realizations} simulations')
        #
        parallel = (cfg.max_submit > 1)          # boolean
        if parallel: pool = mp.Pool(processes=cfg.max_submit)
        #
        for i in range(n_realizations):
            vals = dm.iloc[i, 3:].values
            simdir = cfg.runpath % (i, 0)
            os.makedirs(simdir, exist_ok=True)
            simfile = f'{simdir}/{cfg.result_file}'
            repl = list(zip(keys, vals))
            UT.replace_in_file(repl, cfg.templ_file, simfile)
            #
            logging.info(f'creating {simfile}')
            if parallel:
                pool.apply_async(run_simulation, args=[simfile])
            else:
                run_simulation(simfile)
        #
        if parallel:
            pool.close()
            pool.join()
        print('done')
    else:
        res = main(args.yaml, wake_model=args.wakemodel, knowl_mode=args.knowl)
