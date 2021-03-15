#!/usr/bin/env python

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

        case_nm:                   # if empty, will use basename of *this* file
            Doggerbank

        # model input / parameters
        #
        lut_path:
            ../FugaLUTs/Z0=0.00012000Zi=00790Zeta0=2.00E-7
        tp_A:
            !!float   0.60         # Only used for TurboPark. 'A' parameter in dDw/dx = A*I(x)
        noj_k:
            !!float   0.04         # Only used for Jensen. Wake expansion parameter


        # options
        #
        legend_scaler:
            !!float   0.70         # lower limit of legend is nominal wind speed times this factor
        plot_layout:
            !!bool    false
        plot_wind:
            !!bool    false


 - note2
    on weibull:
        Fuga.exe uses individual weibull for each wtg.
        i have not found an easy way to do this with pywake.
        we have so far used UniformWeibullSite, which obviously does not support this.
        the 'next level' is WaspGridSite, but this is made for complex onshore cases and
        seems a bit complext.
        therefore, for version-1, we just stick to the first weibull given in the knwl_file

 - note3
    on fuga-model
        according to knut seim, fuga cannot be used in cases with more than one type of wtg
        so, it will abort if there is more than one (unique) wtg.
        also, it will need path to fuga look-up tables. (lut_path, see note1 above)

'''

import py_wake
from py_wake.wind_turbines import WindTurbines
from py_wake.site._site import UniformWeibullSite
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys
import glob, re, zipfile, logging, yaml

N_SECTORS = 12
INV_FILE  = 'Inventory.xml'
WWH_FILE  = 'FugaAdapted.wwh'
SEP       = os.path.sep

class Struct(object):
    '''
    just a lazy way of introducing objects
    '''
    pass

def get_default_opts():
    '''
    see note1 for description of attributes
    '''
    opts = Struct()
    opts.case_nm = ''
    opts.inventory_file = 'Inventory.xml'
    opts.output_fnm1    = 'FugaOutput_1.txt'
    opts.output_fnm2    = 'FugaOutput_2.txt'
    opts.tp_A           =  0.60
    opts.noj_k          =  0.04
    opts.legend_scaler  =  0.70
    opts.plot_wakemap   =  False
    opts.plot_layout    =  False
    opts.plot_wind      =  False
    #
    return opts


def  _nparks(sheet):
    a = sheet.iloc[7,1::3]
    return len(a[a.notna()].values)

def _get_weibulls(sheet, nparks):
    wbs = []
    for i in range(nparks):
        wb = Struct()
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

def read_knowl_input(fnm):
    logging.info(f'reading {fnm}')
    knwl = Struct()
    sheet = pd.read_excel(fnm, sheet_name='WindResource')
    np = _nparks(sheet)
    knwl.weibulls = _get_weibulls(sheet, np)
    knwl.turb_intens = sheet.iloc[27+N_SECTORS,1]
    return knwl

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
        wtg = Struct()
        wtg.name       = wtg_nms[i]
        wtg.hub_height = hub_hs[i]
        wtg.diameter   = diams[i]
        wtg.data       = data[ii:ix]
        wtg.ws         = [x[0] for x in wtg.data]
        wtg.pwr        = [x[1] for x in wtg.data]
        wtg.ct         = [x[2] for x in wtg.data]
        wtg.ct_func    = lambda u: np.interp(u, wtg.ws, wtg.ct)
        wtg.pwr_func   = lambda u: np.interp(u, wtg.ws, wtg.pwr)
        #
        wtgs.append(wtg)
        ii = ix
    return wtgs

def _windrose(rose):
    rose = sorted(list(rose))
    wr = Struct()
    wr.dirs = np.array([x[0] for x in rose])
    wr.freqs = np.array([x[1] for x in rose])
    wr.ix = np.argmax(wr.freqs)
    wr.main_dir = wr.dirs[wr.ix]
    return wr

def _get_weibull(weib):
    dirs  = []
    freqs = []
    As    = []
    Ks    = []
    for dir, data in weib.items():
        dirs.append(dir)
        freqs.append(data[0])
        As.append(data[1])
        Ks.append(data[2])
    wb = Struct()
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
    assert('ws_min' not in wtgs.__dict__.keys())    # just in cases ..
    assert('ws_max' not in wtgs.__dict__.keys())
    assert('uniq_wtgs' not in wtgs.__dict__.keys())
    wtgs.ws_min = min([wtg.ws[0] for wtg in wtgs_raw])
    wtgs.ws_max = max([wtg.ws[-1] for wtg in wtgs_raw])
    wtgs.uniq_wtgs = len(set(nms))
    return wtgs

def read_inventory(fnm):
    '''
    reads knowl-inventory (xml) file for setting up windfarm simulation using pywake.
    #
    typical use then:
    case, parks, weibull, wtgs, wtgs_raw  = read_inventory(inventory_file)
    #
    input:
        - fnm      : name of knowl inventory file
    #
    output:
        - case     : all parks collected into one object
        - parks    : list of individual parks
        - wtgs     : an WindTurbines object for all WTG used
        - wtgs_raw : list of individual WTG's as read from the inventory
    '''
    park_nms = []
    ids      = []
    locs     = []
    wtg_data = []
    rose     = set()     # seems to be given reduntantly. dont really know what this is anyway... TODO!
    #weib     = {}        # seems to be given almost reduntantly. TODO!. Not in use - use knowl-input instead
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
        if '<ProductionRoseSector' in line:
            rose.add(_rose_sector(line))
            continue
        #if '<WeibullWind' in line:
            #angl, data = _weibull(line)
            #weib[angl] = data
            #continue
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
    # make a list of each park
    parks = []
    ii = 0
    for i, ix in enumerate(np.concatenate((ixs, [len(ids)]))):     # len(ids) to capture last slice
        p = Struct()
        p.locs  = locs[ii:ix]
        p.size  = len(p.locs)
        p.xs    = np.array([loc[0] for loc in p.locs])
        p.ys    = np.array([loc[1] for loc in p.locs])
        p.ids   = ids[ii:ix]
        p.wtg   = wtgs_raw[i]
        p.types = i * np.ones(p.size, dtype=int)
        p.name  = park_nms[i]
        parks.append(p)
        ii = ix
    #
    # finally, collect all parks into one big case / park
    case = Struct()
    case.locs  = locs
    case.size  = len(locs)
    case.xs    = np.array([loc[0] for loc in locs])
    case.ys    = np.array([loc[1] for loc in locs])
    case.ids   = ids
    case.types = []
    case.wtgs  = []
    case.parks = []
    case.names = []
    for p in parks:
        case.types.extend(p.types)
        case.wtgs.extend([p.wtg]*p.size)
        case.parks.extend([p]*p.size)
        case.names.extend([p.name]*p.size)
    #
    wtgs = _get_windturbines(wtgs_raw)
    #
    #weibull = _get_weibull(weib)
    #
    return case, wtgs, parks, wtgs_raw

def create_output1(net, gross, case, fnm=None):
    '''
    create ouput for knowl
    '''
    unit = net.Description.split('[')[1][:-1]
    assert unit == 'GWh'
    sectors = []
    for l, wd in enumerate(weibull.dirs):
        sector = []
        for i in range(case.size):
            n = net[i][l]
            g = gross[i][l]
            sector.append([case.ids[i], n.x.values, n.y.values, n.h.values, g.sum().values, n.sum().values, case.wtgs[i].name])
        sectors.append(pd.DataFrame(sector))
    # 
    if fnm:
        fms = {0:lambda x: f'Turbine site {x:d}', 1:int, 2:int, 3:int, 4:lambda x: f'{x:1.3f}', 5:lambda x: f'{x:1.3f}'}
        header = f'''FUGA Annual Energy production estimates
z0[m]   99999
zi[m]   99999
zeta0   99999
'''
        f = open(fnm, 'w')
        f.write(header)
        for i, wdr in enumerate(sectors):
            f.write(f'\nSector\t{i}\n')
            f.write(f'      Xpos   Ypos   Height   GrAEP   NetAEP\n')
            f.write(f' Site [m]    [m]    [m]      [GWh]   [GWh]\n')
            f.write(wdr.to_string(index=False, header=False, formatters=fms))
            #
            # write a 'summary' for each park and the total
            n = wdr.iloc[:,4].sum()
            g = wdr.iloc[:,5].sum()
            f.write(f'\nAllProjects                     {n:.4f}   {g:.4f}\n')
            for typeno in sorted(np.unique(case.types)):
                ixs = (case.types==typeno).nonzero()[0]
                p = case.parks[ixs[0]]
                n = wdr.iloc[ixs].iloc[:,4].sum()
                g = wdr.iloc[ixs].iloc[:,5].sum()
                f.write(f'{p.name:20s}             {n:.4f}   {g:.4f}\n')
        f.close()
        logging.info(f'{fnm} was created')
    #
    return sectors

def create_output2(power, case, fnm=None):
    unit = power.Description.split('[')[1][:-1]
    assert unit == 'W'
    vels = []
    for k, wd in enumerate(power.ws):
        vel = []
        for i in range(case.size):
            p = power[i,:,k].values / 1000.       # W => kW
            vel.append([case.ids[i]] + list(p))
        vels.append(pd.DataFrame(vel))
    # 
    if fnm:
        real = lambda x: f'{x:.4f}'
        fms = {0:lambda x: f'Turbine site {x:d}', 1:real, 2:real, 3:real, 4:real, 5:real,
               6:real, 7:real, 8:real, 9:real, 10:real, 11:real}
        binsz = 360. / len(power.wd)
        header = f'''AllProjects
py_wake results
Directional averaging: Simple average ±{binsz/2:.2f} deg
Thrust evaluated by speed at turbine hub
z0[m]   99999
zi[m]   99999
zeta0   99999
Results given as power [kW]
'''
        f = open(fnm, 'w')
        f.write(header)
        dirs = ' '.join(f'{x:.2f}' for x in power.wd.values)
        for i, vel in enumerate(vels):
            f.write(f'\n\nu [m/s] {power.ws.values[i]:.2f}\n')
            f.write(f'Dir [deg]        {dirs}\n')
            f.write(vel.to_string(index=False, header=False, formatters=fms))
        f.write('\n')
        f.close()
        logging.info(f'{fnm} was created')
    #
    return vels

def _unzip(wwh_file, knwl_dir):
    logging.info(f'unzipping {wwh_file}')
    z = zipfile.ZipFile(wwh_file)
    z.extractall(path=knwl_dir)

def get_yaml(fnm):
    '''
    read yaml-file into a Struct for easy access.
    and also avoid the open()
    '''
    yml = yaml.load(open(fnm))
    s = Struct()
    for k,v in yml.items():
        s.__dict__[k] = v
    return s

def get_input(knwl_dir, yml_file):
    opts = get_yaml(yml_file) if yml_file else get_default_opts()
    #
    knwl = read_knowl_input(glob.glob(knwl_dir+SEP+'knowl*.xlsx')[0])
    #
    # read inventory_file
    inv_file = knwl_dir + SEP + INV_FILE
    if not os.path.exists(inv_file):
        wwh_file = knwl_dir + SEP + WWH_FILE
        assert(os.path.exists(wwh_file))
        _unzip(wwh_file, knwl_dir)
    assert(os.path.exists(inv_file))
    return knwl, opts, inv_file 


################################## -- MAIN LOGIC -- ###########################


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    '''
    get necessary input
    '''
    wake_model = sys.argv[1]                                 # Fuga, TP / *Turbo*, or NOJ. TP / *Turbo* = TurboPark, NOJ = Jensen
    knwl_dir   = sys.argv[2] if len(sys.argv) > 2 else '.'   # where to find the knowl-data
    yml_file   = sys.argv[3] if len(sys.argv) > 3 else None  # yaml input file. see note1 above.
    # 
    knwl, opts, inv_file = get_input(knwl_dir, yml_file)
    #
    case, wtgs, _, _  = read_inventory(inv_file)
    #
    weibull = knwl.weibulls[0]                               # for now, we just use the first one. see note2. TODO!

    '''
    pick and initialize the chosen wake model
    '''
    site = UniformWeibullSite(weibull.freqs, weibull.As, weibull.Ks, knwl.turb_intens)
    #
    wake_model = wake_model.upper()                          # convinient
    if wake_model == 'FUGA':
        assert wtgs.uniq_wtgs == 1                           # see note3
        wf_model = py_wake.Fuga(opts.lut_path, site, wtgs)
    elif wake_model == 'TP' or 'TURBO' in wake_model:
        wf_model = py_wake.TP(site, wtgs, k=opts.tp_A)
    elif wake_model == 'NOJ':
        wf_model = py_wake.NOJ(site, wtgs, k=opts.noj_k)
    else:
        raise Exception('The Fuga, TP, or the NOJ wake models are the only options available.')

    '''
    run wake model for all combinations of wd and ws
    '''
    ws = np.arange(np.floor(wtgs.ws_min), np.ceil(wtgs.ws_max)+1)
    sim0 = wf_model(case.xs, case.ys, type=case.types, wd=weibull.dirs, ws=ws)

    '''
    reporting AEP etc
    '''
    net   = sim0.aep(with_wake_loss=True)
    assert(np.all(np.equal(net.x.values, case.xs)))
    assert(np.all(np.equal(net.y.values, case.ys)))
    #
    gross = sim0.aep(with_wake_loss=False)
    g     = gross.values.sum()
    n     = net.values.sum()
    wloss = (g-n) / g * 100                                 # wake loss in %
    #
    print(f'Gross = {g:.0f}')
    print(f'Net   = {n:.0f}')
    print(f'Loss  = {wloss:.1f} %')
    #
    output_fnm1 = knwl_dir + SEP + opts.output_fnm1
    create_output1(net, gross, case, fnm=output_fnm1)
    #
    power = sim0.Power
    #
    assert(np.all(np.equal(power.x.values, case.xs)))
    assert(np.all(np.equal(power.y.values, case.ys)))
    #
    output_fnm2 = knwl_dir + SEP + opts.output_fnm2
    create_output2(power, case, fnm=output_fnm2)

    '''
    optional stuff
    '''
    if opts.plot_wakemap:
        # plot flow map for the dominant wind direction / wind-speed
        # grid=None defaults to HorizontalGrid(resolution=500, extend=0.2)
        plt.figure()
        ws_plot = int(weibull.main_ws)
        ws1, ws2 = np.floor(opts.legend_scaler*ws_plot), ws_plot+1
        levels = np.linspace(ws1, ws2, int((ws2-ws1)*10)+1)
        flow_map = sim0.flow_map(grid=None, wd=weibull.main_dir, ws=ws_plot)
        flow_map.plot_wake_map(levels=levels, plot_ixs=False)
        plt.title(f'{opts.case_nm} :: {opts.wake_model} :: {ws_plot:d} m/s :: {weibull.main_dir:.0f} deg')
        plt.show()
    #
    if opts.plot_layout:
        plt.figure()
        wtgs.plot(case.xs, case.ys)
        plt.show()
    #
    if opts.plot_wind:
        plt.figure()
        site.plot_wd_distribution(n_wd=12, ws_bins=[0,5,10,15,20,25])
        plt.show()

