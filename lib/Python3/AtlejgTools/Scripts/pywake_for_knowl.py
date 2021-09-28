#!/usr/bin/env python3

'''
this is a demo to show that pywake could be used for simulation of knowl-cases.
the idea is that knowl is used for setting up the case with layout, wtg's etc, and in stead of calling
for example Fuga, it could call pywake.

the information from knowl is transferred using the Inventory.xml (un-zipped from FugaAdapted.wwh).
it will need to get the knowl data-directory where it will look for Inventory.xml. if it does
not find it, it will look for FugaAdapted.wwh and un-zip it.

i tried to use xml-reader to parse Inventory.xml:
import xml.etree.ElementTree as ET
r = ET.parse('Inventory.xml').getroot()
unfortunately, Inventory.xml is not "proper" xml as it does not close all it's tags (<tagname> ... </tagname>)
so the result is a strongly nested structure that is almost impossible to navigate.
(i made a cleanup-version Inventory2.xml that worked, but it was a dead end anyway).

so, I ended up parsing Inventory.xml myself, in a quite naiive way.

there are still a few things to be worked out:
    - it seems weibull-parameters change from location to location. why?
      (i dont pick up this, i only use one set of weilbull-paramters)
if pywake is to be called from knowl, i suggest that we make a more useful way of transferring the data needed.

NOTES

 - note1
    the yaml-input file in the __main__ part should look like this:

        case_nm: Doggerbank        # if empty, will use basename of *this* file
        inventory_file:            # if empty, will use Inventory.xml

        # model input / parameters
        #
        lut_path:        !!str    ../FugaLUTs/Z0=0.00012000Zi=00790Zeta0=2.00E-7    # path to look-up tables (Fuga only)
        tp_A:            !!float  0.60                # Only used for TurboPark. 'A' parameter in dDw/dx = A*I(x)
        noj_k:           !!float  0.04                # Only used for Jensen. Wake expansion parameter
        delta_winddir:   !!float  1.0                 # delta wind-dir for calculations
        delta_windspeed: !!float  0.50                # delta wind-vel for calculations

        # output file-names
        output_fnm1:     !!str    FugaOutput_1.txt    # the file needed when running 'Wake only' from knowl

        # options
        #
        plot_wakemap:    !!bool   false
        plot_layout:     !!bool   false
        plot_wind:       !!bool   false
        legend_scaler:   !!float  0.70                # lower limit of legend is nominal wind speed times this factor

    SHORT VERSION:
        inventory_file:  !!str     Inventory.xml
        knowl_file:      !!str     knowl_v4.5_input.xlsx

 - note2
    on weibull:
        Fuga.exe uses individual weibull for each wtg.
        i have not found an easy way to do this with pywake.
        we have so far used UniformWeibullSite, which obviously does not support this.
        the 'next level' is WaspGridSite, but this is made for complex onshore cases and
        seems a bit tricky.
        therefore, for version-1, we just stick to one weibull given in the knowl_file. 
        this is given by opts.weibull_index (which starts on 1).
        update: i have now implemented the use of a wrg-file

 - note3
    on fuga-model
        according to knut seim, fuga cannot be used in cases with more than one type of wtg
        so, it will abort if there is more than one (unique) wtg.
        also, it will need path to fuga look-up tables. (lut_path, see note1 above)

 - note4
    naming conventions follows typical pywake lingo:
        - wt: wind turbine
          wd: wind direction
          ws: wind speed
          wl: wake loss

 - note5
    about the different TurbOPark-versions, see mail from Knut Seim, 20/8-21

'''

import argparse
import py_wake
from py_wake.wind_turbines import WindTurbines
from py_wake.site._site import UniformWeibullSite
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import os, sys, time, pickle
import multiprocessing as mp
import glob, re, zipfile, logging, yaml
from scipy.interpolate import interp1d
import AtlejgTools.WindResourceAssessment.Utils as WU
import AtlejgTools.Utils as UT
import AtlejgTools.Scripts.run_pywake as run_pywake

N_SECTORS = 12
WWH_FILE  = 'FugaAdapted.wwh'  # the zipped inventory file (for Fuga)
SEP       = os.path.sep
EPS       = 1e-9               # small non-zero value


def set_default_opts(opts):
    '''
    see note1 for description of attributes
    '''
    if not hasattr(opts, 'case_nm'):            opts.case_nm            = ''
    if not hasattr(opts, 'inventory_file'):     opts.inventory_file     = 'Inventory.xml'
    if not hasattr(opts, 'output_fnm1'):        opts.output_fnm1        = 'FugaOutput_1.txt'
    if not hasattr(opts, 'tp_A'):               opts.tp_A               =  0.60
    if not hasattr(opts, 'noj_k'):              opts.noj_k              =  0.04
    if not hasattr(opts, 'legend_scaler'):      opts.legend_scaler      =  0.70
    if not hasattr(opts, 'plot_wakemap'):       opts.plot_wakemap       =  False
    if not hasattr(opts, 'plot_layout'):        opts.plot_layout        =  False
    if not hasattr(opts, 'plot_wind'):          opts.plot_wind          =  False
    if not hasattr(opts, 'delta_winddir'):      opts.delta_winddir      = 1.
    if not hasattr(opts, 'delta_windspeed'):    opts.delta_windspeed    = 0.5
    if not hasattr(opts, 'weibull_index'):      opts.weibull_index      = 1
    if not hasattr(opts, 'shift_wind_map'):     opts.shift_wind_map     = False      # mostly used for development
    if not hasattr(opts, 'report_aep'):         opts.report_aep         = True       # mostly used for development
    if not hasattr(opts, 'logfile') or not opts.logfile: opts.logfile   = None
    #
    if not hasattr(opts, 'knowl_file'): opts.knowl_file = glob.glob('knowl_v*input.xlsx')[0]
    #
    return opts


def  _nparks(sheet):
    a = sheet.iloc[7,1::3]
    return len(a[a.notna()].values)

def _get_weibulls(sheet):
    wbs = []
    for i in range(_nparks(sheet)):
        wb = UT.Struct()
        binsz = 360. / N_SECTORS
        wb.dirs = binsz*np.arange(N_SECTORS)
        wb.As = np.array(sheet.iloc[7:7+N_SECTORS, 3*i+1], dtype=float)
        wb.Ks = np.array(sheet.iloc[7:7+N_SECTORS, 3*i+2], dtype=float)
        wb.freqs = np.array(sheet.iloc[7:7+N_SECTORS, 3*i+3], dtype=float) / 100.
        wb.ix = np.argmax(wb.freqs)
        wb.main_dir = wb.dirs[wb.ix]
        wb.main_ws = wb.As[wb.ix]
        wbs.append(wb)
    return wbs

def read_knowl_file(knowl_file):
    logging.info(f'reading {knowl_file}')
    knowl = UT.Struct()
    sheet = pd.read_excel(knowl_file, sheet_name='WindResource')
    knowl.weibulls = _get_weibulls(sheet)
    # see knowl_v*.m for location of parameters (search for j0, j1, j2...)
    knowl.turb_intens = sheet.iloc[39,1]
    knowl.weiba_fnm = sheet.iloc[71,1]
    knowl.weibk_fnm = sheet.iloc[72,1]
    try:
        knowl.wrg_file  = knowl_dir+SEP+sheet.iloc[76,1]
        if not type(knowl.wrg_file) is str: knowl.wrg_file = None
    except:
        knowl.wrg_file  = None
    return knowl

def _clean(line):
    line = line.replace('=" ', '="')     # remove spaces that mess things up for some fields
    return line

def _location(line):
    line = _clean(line)
    x, y = [float(x.split('=')[1].replace('"','')) for x in line.split()[1:-1]]
    return x, y

def _wtg_data(line):
    line = _clean(line)
    ws, power, ct = [float(x.split('=')[1].replace('"','')) for x in line.split()[1::2]]
    return ws, power, ct

def _rose_sector(line):
    line = _clean(line)
    index, angle, freq = [float(x.split('=')[1].replace('"','')) for x in line.split()[1:-4]]
    return angle, freq

def _weibull(line):
    line = _clean(line)
    angle, _, freq, wa, wk = [float(x.split('=')[1].replace('"','')) for x in line.split()[2:-1]]
    return angle, (freq, wa, wk)

def _park_nm(line):
    line = _clean(line)
    return line.split('"')[5]

def _wtg_info(line):
    line = _clean(line)
    recs = line.split('"')
    return recs[3], float(recs[9])

def _hub_height(line):
    line = _clean(line)
    return float(line.split('>')[1].split('<')[0].strip())

def _wtgs(data, wtg_nms, hub_hs, diams):
    ws  = [x[0] for x in data]
    ixs = (np.diff(ws) < 0).nonzero()[0] + 1   # indices indicating start of new wtg
    wtgs = []
    ii = 0
    for i, ix in enumerate(np.concatenate((ixs, [len(data)]))):     # len(data) to capture last slice
        wtg = UT.Struct()
        wtg.name       = wtg_nms[i]
        wtg.hub_height = hub_hs[i]
        wtg.diameter   = diams[i]
        wtg.data       = data[ii:ix]
        wtg.ws         = np.array([x[0] for x in wtg.data])
        # use EPS to avoid hard zeros
        wtg.pwr        = np.maximum([x[1] for x in wtg.data], EPS)
        wtg.ct         = np.maximum([x[2] for x in wtg.data], EPS)
        wtg.ct_func    = interp1d(wtg.ws, wtg.ct,  bounds_error=False, fill_value=(EPS,EPS))
        wtg.pwr_func   = interp1d(wtg.ws, wtg.pwr, bounds_error=False, fill_value=(EPS,EPS))
        #
        wtgs.append(wtg)
        ii = ix
    return wtgs

def _create_weibull(weib):
    dirs  = []
    freqs = []
    As    = []
    Ks    = []
    for dir, data in weib.items():
        dirs.append(dir)
        freqs.append(data[0])
        As.append(data[1])
        Ks.append(data[2])
    wb = UT.Struct()
    wb.dirs  = np.array(dirs)
    wb.freqs = np.array(freqs)
    wb.As    = np.array(As)
    wb.Ks    = np.array(Ks)
    wb.ix = np.argmax(wb.freqs)
    wb.main_dir = wb.dirs[wb.ix]
    wb.main_ws = wb.As[wb.ix]
    return wb

def _get_windturbines(wtgs_raw):
    nms = [wtg.name       for wtg in wtgs_raw]
    ds  = [wtg.diameter   for wtg in wtgs_raw]
    hhs = [wtg.hub_height for wtg in wtgs_raw]
    cfs = [wtg.ct_func    for wtg in wtgs_raw]
    pfs = [wtg.pwr_func   for wtg in wtgs_raw]
    #
    wtgs =  WindTurbines(names=nms, diameters=ds, hub_heights=hhs, ct_funcs=cfs, power_funcs=pfs, power_unit='W')
    #
    # add min and max windspeed of all wtgs
    assert(not hasattr(wtgs, 'ws_min'))    # just in cases ..
    assert(not hasattr(wtgs, 'ws_max'))
    assert(not hasattr(wtgs, 'uniq_wtgs'))
    wtgs.ws_min = min([wtg.ws[0] for wtg in wtgs_raw])
    wtgs.ws_max = max([wtg.ws[-1] for wtg in wtgs_raw])
    wtgs.uniq_wtgs = len(set(nms))
    return wtgs

def read_inventory(fnm, selected=[]):
    '''
    reads knowl-inventory (xml) file for setting up windfarm simulation using pywake.
    #
    typical use then:
    case, wtgs = read_inventory(inventory_file)
    #
    input:
        - fnm      : name of knowl inventory file
    #
    output:
        - case     : all parks collected into one object.
                     includes lists of individual WTG's and parks as read from the inventory (wtg_list and park_list)
        - parks    : list of individual parks
        - wtgs     : an WindTurbines object for all WTG used
    '''
    #
    # making sure inventory_file is available
    if not os.path.exists(fnm):
        assert(os.path.exists(WWH_FILE))
        _unzip(WWH_FILE, '.')
    assert(os.path.exists(fnm))
    #
    park_nms = []
    ids      = []
    locs     = []
    wtg_data = []
    #rose     = set()     # seems to be given reduntantly. dont really know what this is anyway... TODO!
    weib     = {}
    weibs    = []
    hub_hs   = []
    diams    = []
    wtg_nms  = []
    #
    logging.info(f'reading {fnm}')
    lines = open(fnm).readlines()
    for line in lines:
        if 'Turbine site group' in line:
            park_nms.append(_park_nm(line))
            continue
        if '<Location x-Location' in line:
            locs.append(_location(line))
            continue
        if '<DataPoint' in line:
            wtg_data.append(_wtg_data(line))
            continue
        if '<WindTurbineGenerator' in line:
            nm, diam = _wtg_info(line)
            wtg_nms.append(nm)
            diams.append(diam)
            continue
        if '<Height' in line:
            hub_hs.append(_hub_height(line))
            continue
        #if '<ProductionRoseSector' in line:
        #    rose.add(_rose_sector(line))
        #    continue
        if '<WeibullWind' in line:
            angl, data = _weibull(line)
            weib[angl] = data
            continue
        if '</RveaWeibullWindRose' in line:
            weibs.append(_create_weibull(weib))
            weib = {}
        m = re.search('"Turbine site (\d\d\d\d)"', line)
        if m:
            ids.append(int(m.groups()[0]))
            continue
    #
    park_nms = [x for x in park_nms if x != 'AllProjects']
    #
    wtgs_raw = _wtgs(wtg_data, wtg_nms, hub_hs, diams)
    #
    ixs = (np.diff(ids) > 1).nonzero()[0] + 1    # indices indicating start of new park-area. .. 1035, 2001, ..
    #
    wtgs_ = []  # used (selected) wtgs
    locs_ = []  # used (selected) locs
    nms_  = []  # used (selected) names
    ids_  = []  # used (selected) ids
    #
    # each park is a slice of the lists collected
    # => make a list of parks
    parks = []
    ii = 0
    type_i = -1      # need to make sure type's are 0,0,0, ... ,1,1,1, ...
    for i, ix in enumerate(np.concatenate((ixs, [len(ids)]))):     # len(ids) to capture last slice
        #
        if selected and park_nms[i] not in selected:
            ii = ix
            continue
        #
        type_i += 1
        p = UT.Struct()
        p.locs  = locs[ii:ix]
        p.size  = len(p.locs)
        p.xs    = np.array([loc[0] for loc in p.locs])
        p.ys    = np.array([loc[1] for loc in p.locs])
        p.ids   = ids[ii:ix]
        p.wtg   = wtgs_raw[i]
        p.types = type_i * np.ones(p.size, dtype=int)
        p.name  = park_nms[i]
        parks.append(p)
        #
        # make sure we use only the selected data below
        wtgs_.append(p.wtg)
        locs_.extend(p.locs)
        nms_.append(p.name)
        ids_.append(p.ids)
        # housekeeping
        logging.info(f'Park: {p.name}')
        logging.info(f'WTG: {p.wtg.name}')
        ii = ix
    #
    # finally, collect all parks into one big case / park
    case = UT.Struct()
    case.locs  = locs_
    case.size  = len(locs_)
    case.xs    = np.array([loc[0] for loc in locs_])
    case.ys    = np.array([loc[1] for loc in locs_])
    case.ids   = ids_
    case.wtgs  = []
    case.parks = []
    case.names = []
    case.weibs = weibs
    types = []
    for p in parks:
        types.extend(p.types)
        case.wtgs.extend([p.wtg]*p.size)
        case.parks.extend([p]*p.size)
        case.names.extend([p.name]*p.size)
    case.wtg_names = [wtg.name for wtg in case.wtgs]
    case.hub_heights = [wtg.hub_height for wtg in case.wtgs]
    case.types = np.array(types)
    #
    # some useful stuff
    case.wtg_list  = wtgs_
    case.park_list = parks
    case.n_parks   = len(parks)
    case.pwr_funcs = [wtg.pwr_func for wtg in case.wtg_list]
    case.park_nms  = nms_
    #
    wtgs = _get_windturbines(wtgs_)
    return case, wtgs

def _read_header(fnm):
    hdr = pd.read_csv(fnm, nrows=6, header=None, delim_whitespace=True)
    keys = hdr.T.iloc[0,:].values
    vals = hdr.T.iloc[1,:].values
    hd = dict(zip(keys, vals))
    hd['ncols'] = int(hd['ncols'])
    hd['nrows'] = int(hd['nrows'])
    return hd

def get_site(weiba_fnm, weibk_fnm, mwb):
    '''
    this function should solve note2
        what i've learned:
         - neet to use WaspGridSite for variable Weibulls
           * requires regular grid of the weibulls
           * this is what is found in wind-maps like vortex.a.140.asc
        what remains:
         - not sure how the vortex-files handles wind-direction or frequencies.
    - input
      * weiba_fnm: name of weibull A-parameter file (found in knowl-excel)
      * weibk_fnm: name of weibull k-parameter file (found in knowl-excel)
      * mwb      : main weibull (fron knowl excel input)
    '''
    h = float(re.search('(\d+)', weiba_fnm).group())   # extract h from file-name (typically vortex.a.140.asc)
    hd = _read_header(weiba_fnm)
    a = pd.read_csv(weiba_fnm, skiprows=6, header=None, delim_whitespace=True)
    k = pd.read_csv(weibk_fnm, skiprows=6, header=None, delim_whitespace=True)
    x = hd['xllcorner'] + hd['cellsize']*np.arange(hd['ncols'])
    y = hd['yllcorner'] + hd['cellsize']*np.arange(hd['nrows'])
    binsz = 360. / N_SECTORS
    wds = binsz*np.arange(N_SECTORS)
    dims   = ["h", "wd", "x", "y"]
    coords = dict(h=[h],
                  wd=mwb.dirs,
                  x=x,
                  y=y,
                 )
    a_map = np.tile(a, [len(mwb.dirs), 1])
    wba = xr.DataArray(data=[a_map], dims=dims, coords=coords)
    stop

def _unzip(wwh_file, todir):
    logging.info(f'unzipping {wwh_file}')
    z = zipfile.ZipFile(wwh_file)
    z.extractall(path=todir)

def get_yaml(fnm):
    '''
    read yaml-file into a Struct for easy access.
    '''
    yml = yaml.load(open(fnm), Loader=yaml.SafeLoader)
    s = UT.Struct()
    for k,v in yml.items():
        s.__dict__[k] = v
    return s

def get_opts(yaml_file):
    opts = get_yaml(yaml_file) if yaml_file else UT.Struct()
    set_default_opts(opts)
    return opts

def _dump(res, fnm):
    f = open(fnm, 'wb')
    pickle.dump(res, f)
    f.close()
    logging.info(f'dumped results to {fnm}')

def load(fnm):
    '''
    loading pickled results (from earlier simulation)
    '''
    f = open(fnm, 'rb')
    res = pickle.load(f)
    logging.info(f'loaded results from {fnm}')
    return res

def plot_flowmap(sim_res, ws0, ws_min, ws_max, wd0, wake_model, case_nm, plot_wt):
    logging.info('Plotting wakemap')
    plt.figure()
    levels = np.linspace(ws_min, ws_max, int((ws_max-ws_min)*10)+1)
    grid =  py_wake.HorizontalGrid(resolution=1500, extend=0.1) # grid=None defaults to HorizontalGrid(resolution=500, extend=0.2)
    flow_map = sim_res.flow_map(grid=grid, wd=wd0, ws=ws0)
    flow_map.plot_wake_map(levels=levels, plot_windturbines=plot_wt)    # , plot_ixs=False)
    plt.title(f'{case_nm} :: {wake_model} :: {ws0:.0f} m/s :: {wd0:.0f} deg')
    plt.show()


def main(wake_model, yaml_file=None, selected=[], dump_results=True):
    '''
    pick up knowl case description and run PyWake simulation.
    if wake_model is None, it *must* be given in the yaml_file
    - input
      * wake_model: NOJ, NOJLOCAL, ETP ...
      * yaml_file:  see note1
      * selected:   Which parks to use. Default is [], which gives all. This is useful since it means
                    only one Inventory.xml is needed (covering all parks in the knowl-project)
    - returns
      * sim_res
      * sim
      * aeps
      * sim
      * case
      * knowl
      * opts
      * wtgs
      * site
      * wf_model
    '''
    #
    opts = get_opts(yaml_file)
    #
    knowl = read_knowl_file(opts.knowl_file)
    case, wtgs = read_inventory(opts.inventory_file, selected=selected)
    #
    logging.basicConfig(level=logging.INFO, filename=opts.logfile)
    tic = time.perf_counter()
    #
    # setup things
    if not knowl.wrg_file:
        weib = knowl.weibulls[opts.weibull_index-1]            # we just use one. see note2. TODO!
        site = UniformWeibullSite(weib.freqs, weib.As, weib.Ks, knowl.turb_intens)
    else:
        wrg = WU.read_wrg(knowl.wrg_file, tke=knowl.turb_intens)
        if opts.shift_wind_map:                                # for testing/development only
            wrg = wrg.assign_coords(z=case.hub_heights[:1])
        wrg = wrg.rename(A='Weibull_A', k='Weibull_k', f='Sector_frequency', sec='wd', z='h', tke='TI')
        wrg = wrg.assign_coords(wd=(wrg.wd-1)*30)
        site = py_wake.site.XRSite(wrg)
    #
    # pick and initialize the chosen wake model
    if not wake_model: wake_model = opts.wake_model
    wake_model = wake_model.upper()
    wf_model = run_pywake.deficit_model(wake_model, opts, site, wtgs)
    #
    # run simulations
    logging.info(f'run wake model {wake_model} for all combinations of wd and ws')
    wd = np.arange(0, 360, opts.delta_winddir)
    ws = np.arange(np.floor(wtgs.ws_min/2), np.ceil(wtgs.ws_max)+1, opts.delta_windspeed)
    sim_res = wf_model(case.xs, case.ys, type=case.types, wd=wd, ws=ws)
    assert(np.all(np.equal(sim_res.x.values, case.xs)))
    assert(np.all(np.equal(sim_res.y.values, case.ys)))
    #
    if opts.report_aep:
        #
        # coarsen it by averaging
        n_sectors = len(weib.freqs)
        width  = 360. / n_sectors
        n_bins    = int(width / opts.delta_winddir)         # number of bins per sector
        sim = sim_res.coarsen(wd=n_bins).mean()
        #
        # reporting AEP etc
        if not hasattr(opts, 'output_fnm1') or not opts.output_fnm1:
            fnm = wake_model + '.txt'
        else:
            fnm = opts.output_fnm1
        #
        aeps = WU.calc_AEP(sim, wtgs, weib, park_nms=case.park_nms, verbose=True)
        WU.create_output(aeps, case, fnm)
    else:
        aeps = None
        sim = sim_res
    #
    # optional stuff
    if opts.plot_wakemap:
        ws_min, ws_max = np.floor(opts.legend_scaler*opts.plot_ws), opts.plot_ws+1
        plot_flowmap(sim_res, opts.plot_ws, ws_min, ws_max, opts.plot_wd, wake_model, opts.case_nm, opts.plot_wt)
    #
    if opts.plot_layout:
        logging.info('Plotting layout')
        plt.figure()
        wtgs.plot(case.xs, case.ys)
        plt.show()
    #
    if opts.plot_wind:
        logging.info('Plotting wind-distribution')
        plt.figure()
        site.plot_wd_distribution(n_wd=12, ws_bins=[0,5,10,15,20,25])
        plt.show()
    #
    toc = time.perf_counter()
    logging.info(f"Total runtime: {toc-tic:0.1f} seconds")
    #
    res =  sim_res, sim, aeps, case, knowl, opts, wtgs, site, wf_model
    if dump_results: _dump(res, f'{wake_model}.pck')
    return res

def run_single(directory, wake_model, yaml_file, selected=[]):
    cwd = os.getcwd()
    if not os.path.exists(directory):
        print('creating', directory)
        os.mkdir(directory)
    os.chdir(directory)
    print(directory, wake_model, yaml_file)
    main(wake_model, yaml_file=yaml_file, selected=selected)
    os.chdir(cwd)

def read_casefile(csvfile, sep='&', comment='#'):
    '''
    see docs for run_multiple
    '''
    cases = pd.read_csv(csvfile, sep=';', comment='#')
    parklist = []
    for i, case in cases.iterrows():
        if case.parks is np.NaN:
            parks = []    # get all parks
        else:
            parks = [park.strip() for park in case.parks.split('&')]
        parklist.append(parks)
    cases.parks = parklist
    #
    return cases

def run_multiple(csvfile, max_cpus=3):
    '''
    useful for running a set of simulations.
    the csv-file should look like something like this (to be opened in excel!):
        directory;parks;wake_model;yaml
        UKE_int;DEP NW & DEP SE & DEP & SEP;ETP;../../2.yaml
        UKE_ext;DEP NW & DEP SE & DEP & SEP & DOW & ShS & RaceBank & Triton Knoll;ETP;../../2.yaml
        UKE_fut;DEP NW & DEP SE & DEP & SEP & DOW & ShS & RaceBank & Triton Knoll & R4;ETP;../../2.yaml
        UKE_int;DEP NW & DEP SE & DEP & SEP;NOJ;../../2.yaml
        UKE_ext;DEP NW & DEP SE & DEP & SEP & DOW & ShS & RaceBank & Triton Knoll;NOJ;../../2.yaml
        UKE_fut;DEP NW & DEP SE & DEP & SEP & DOW & ShS & RaceBank & Triton Knoll & R4;NOJ;../../2.yaml
        UKE_int;DEP NW & DEP SE & DEP & SEP;NOJLOCAL;../../2.yaml
        UKE_ext;DEP NW & DEP SE & DEP & SEP & DOW & ShS & RaceBank & Triton Knoll;NOJLOCAL;../../2.yaml
        UKE_fut;DEP NW & DEP SE & DEP & SEP & DOW & ShS & RaceBank & Triton Knoll & R4;NOJLOCAL;../../2.yaml
    '''
    cases = read_casefile(csvfile)
    pool = mp.Pool(max_cpus)
    results = []
    #
    for i, case in cases.iterrows():
        pool.apply_async(run_single, args=(case.directory, case.wake_model, case.yaml, case.parks))
    #
    pool.close()
    pool.join()
    print('done')

def run_multiple_old(csvfile, max_cpus=3):
    '''
    useful for running a set of simulations.
    the csv-file should look like something like this:
        directory                wake_model   yaml
        #DOW_int                  ETP       ../1.yaml
        #DOW_ext                  ETP       ../1.yaml
        #DOW_fut                  ETP       ../1.yaml
        DOW_int                  NOJ       ../1.yaml
        DOW_ext                  NOJ       ../1.yaml
        DOW_fut                  NOJ       ../1.yaml
        DOW_int                  NOJLOCAL  ../1.yaml
        DOW_ext                  NOJLOCAL  ../1.yaml
        DOW_fut                  NOJLOCAL  ../1.yaml
    '''
    cases = pd.read_csv(csvfile, delim_whitespace=True, comment='#')
    pool = mp.Pool(max_cpus)
    results = []
    #
    for i, case in cases.iterrows():
        pool.apply_async(run_single, args=(case.directory, case.wake_model, case.yaml))
    #
    pool.close()
    pool.join()
    print('done')


################################## -- MAIN LOGIC -- ###########################

if __name__ == '__main__':

    #
    # get necessary input
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--wake_model", default='NOJ',
                        help="wake-model: Fuga, TP, ETP, NOJ, NOJLOCAL. (TP=TurboPark, ETP=Equinor-TP, NOJ=NO-Jensen). Required")
    parser.add_argument("-y", "--yaml_file",   default=None,
                        help="name of yaml-file of options for this program. Optional")
    args = parser.parse_args()
    #
    # run program
    sim_res, sim, aeps, case, knowl, opts, wtgs, site, wf_model = \
        main(args.wake_model, yaml_file=args.yaml_file)

