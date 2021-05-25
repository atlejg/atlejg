import numpy as np
import xarray as xr
import pandas as pd
import logging
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import dateutil

REAL = lambda x: f'{x:.4f}'
INT  = lambda x: f'{x:.0f}'


def _binArray(data, axis, binstep, binsize, func=np.nanmean):
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis]= argdims[axis], argdims[0]
    data = data.transpose(argdims)
    data = [func(np.take(data,np.arange(int(i*binstep),int(i*binstep+binsize)),0),0) for i in np.arange(dims[axis]//binstep)]
    data = np.array(data).transpose(argdims)
    return data

def G_func(alpha, k):
    # 1/k times the incomplete gamma function of the two arguments 1/k and alpha^k
    # Note, the scipy incomplete gamma function, gammainc, must be multiplied with gamma(k) to match the
    # the G function used in the European Wind Atlas
    import scipy.special as sc
    return 1 / k * sc.gamma(1 / k) * sc.gammainc(1 / k, alpha**k)

def _weigthed_pwr(u, power, wb_a, wb_k):
    """Calculate the weibull weighted power
    Parameters
    ----------
    Power : xarray DataArray
        Power
    Returns
    -------
    y : array_like
    """
    # see https://orbit.dtu.dk/en/publications/european-wind-atlas, page 95
    #
    u0, u1 = u[:-1], u[1:]
    #alpha0, alpha1 = (u0 / wb_a)+np.random.rand()*1e-10, (u1 / wb_a)+np.random.rand()*1e-10      # knut s. used this. i try without
    alpha0, alpha1 = (u0 / wb_a), (u1 / wb_a)
    p0, p1 = power[:, :-1], power[:, 1:]
    #
    wpw = (    p0 * np.exp(-alpha0**wb_k)                                                                          # eq 6.5, p0 * cdf(u0:)
             + (p1 - p0) / (alpha1 - alpha0) * (G_func(alpha1, wb_k) - G_func(alpha0, wb_k))                     # eq 6.4 linear change p0 to p1
             - p1 * np.exp(-alpha1**wb_k)                                                                          # eq 6.5, - p1 * cdf(u1:)
    )
    #
    return wpw

def calc_power(sim, wtgs, weib, park_nms=[]):
    '''
    calculates weibull-averaged net and gross power for each sector.
    # 
    loosely based on scriptified version of knut s. seim's notebook:
    /cygdrive/c/Appl/atlejg/lib/Python3/AtlejgTools/WindResourceAssessment/KnutSeim_stuff/py_wake_demo2.py
    #
    notes: 
        n1: assumes one WTG-type for each park
        n2: assumes that each park has it's own type-number
        n3: assumes one weibull for the whole area
    # 
    input:
        sim  : simulation result from pywake. have typically been coarsen'ed
        wtgs : wtgs
        weib : weibull-struct for the entire area (n3)
    output:
        net & gross power per park in GW (per WTG & sector & wind-speed)
    '''
    #
    wb_K = interp1d(weib.dirs, weib.Ks, kind='nearest', fill_value='extrapolate')
    wb_A = interp1d(weib.dirs, weib.As, kind='nearest', fill_value='extrapolate')
    #
    pwrs = []
    for i in np.unique(sim.type.values):                                                                             # loop each park
        #
        pwr = sim.where(sim.type==i, drop=True).Power
        nm = park_nms[i] if park_nms else ''
        #
        gross  = np.tile(wtgs.power(pwr.ws, pwr.type), [pwr.sizes['wt'], pwr.sizes['wd'], 1])                         # gross power - from WTG power curve
        #
        # weighting according to European Wind Atlas using the weibull
        net_w   = np.zeros(pwr.values.shape)
        gross_w = np.zeros(pwr.values.shape)
        for n, wd in enumerate(pwr.wd.values):
            net_w[:,n,:-1]   = _weigthed_pwr(pwr.ws.values, pwr.values[:,n,:], wb_A(wd), wb_K(wd))
            gross_w[:,n,:-1] = _weigthed_pwr(pwr.ws.values, gross[:,n,:],      wb_A(wd), wb_K(wd))
        #
        na = xr.DataArray(data=net_w,   dims=pwr.dims, coords=pwr.coords, attrs=dict(description=f'Net power {nm}'))
        ga = xr.DataArray(data=gross_w, dims=pwr.dims, coords=pwr.coords, attrs=dict(description=f'Gross power {nm}'))
        pwrs.append([na, ga])
    #
    return pwrs

def calc_AEP(sim, wtgs, weib, park_nms=[], verbose=False):
    '''
    calculates weibull-averaged net and gross AEP for each sector.
    # 
    loosely based on scriptified version of knut s. seim's notebook:
    /cygdrive/c/Appl/atlejg/lib/Python3/AtlejgTools/WindResourceAssessment/KnutSeim_stuff/py_wake_demo2.py
    #
    notes: 
        n1: assumes one WTG-type for each park
        n2: assumes that each park has it's own type-number
        n3: assumes one weibull for the whole area
    # 
    input:
        sim  : simulation result from pywake
        wtgs : wtgs
        weib : weibull-struct for the entire area (n3)
    output:
        net & gross AEP per park in GWh (per WTG & sector)
    '''
    #
    wb_f = interp1d(weib.dirs, weib.freqs, kind='nearest', fill_value='extrapolate')
    cf = 24*365*1e-9                                                                                               # conversion factor => GWh. 365.24?? TODO
    #
    pwrs = calc_power(sim, wtgs, weib, park_nms=park_nms)
    aeps = []
    if not park_nms: park_nms = [''] * len(aeps)
    #
    for nm, (net, gross) in zip(park_nms, pwrs):                                                                   # loop each park
        naep, gaep = net.sum('ws'), gross.sum('ws')
        for wd in sim.wd.values:
            naep.sel(wd=wd).values *= wb_f(wd) * cf
            gaep.sel(wd=wd).values *= wb_f(wd) * cf
        aeps.append([naep, gaep])
        #
        if verbose:
            n = naep.sum().values
            g = gaep.sum().values
            wloss = (g-n)/g * 100
            #
            print(f'{nm} Gross AEP (GWh) {g:.1f}')
            print(f'{nm} Net AEP (GWh) {n:.1f}')
            print(f'{nm} Wake loss (%) {wloss:.2f}')
        #
    return aeps

def calc_AEP_old(sim0, pwr_funcs, weibs, dwd=1., park_nms=[], verbose=False, return_pwr=False):
    '''
    calculates averaged AEP for each sector.
    # 
    based on scriptifie version of knut s. seim's notebook:
    /cygdrive/c/Appl/atlejg/lib/Python3/AtlejgTools/WindResourceAssessment/KnutSeim_stuff/py_wake_demo2.py
    #
    notes: 
        n1: assumes one WTG-type & weibull for each park
        n2: only regular sectors supported
    # 
    input:
        sim0     : simulation result from pywake
        pwr_funcs: power function for each WTG (see n1)
        weibs    : weibull per sector per park. must have attributes weib.Ks, weib.As, weib.freqs (see n2)
                   note that sector width is given by number of entries in weibull (see n3)
        dwd      : delta wd. [deg]. gives number of bins.
    output:
        net & gross AEP per park in GWh (per WTG & sector)
    '''
    #
    aeps = []
    pwrs = []
    #
    for i, (pwr_func, weib) in enumerate(zip(pwr_funcs, weibs)):                                                       # loop each park
        #
        nm = park_nms[i] if park_nms else ''
        #
        sim = sim0.where(sim0.type==i, drop=True)
        #
        n_sectors = len(weib.freqs)
        sector_w  = 360. / n_sectors
        n_bins    = int(sector_w / dwd)                                                                                # number of bins per sector
        #
        # power per WTG and wind-direction (wt, wd)
        npwrb = _binArray(sim.Power.values,1, n_bins, n_bins)                                                          # net power from simulation. binned
        gpwr  = np.tile(pwr_func(sim.ws), [sim.Power.shape[0], sim.Power.shape[1], 1])                                 # gross power - from WTG power curve
        gpwrb = _binArray(gpwr,1, n_bins, n_bins)                                                                      # gross power - from WTG power curve. binned
        #
        # weighting of bins
        nwpw = np.ndarray((len(sim.wt), n_sectors, len(sim.ws)-1))
        gwpw = np.ndarray((len(sim.wt), n_sectors, len(sim.ws)-1))
        for n in range(n_sectors):
            nwpw[:,n,:] = _weigthed_pwr(sim.ws.values, npwrb[:,n,:], weib.As[n], weib.Ks[n])
            gwpw[:,n,:] = _weigthed_pwr(sim.ws.values, gpwrb[:,n,:], weib.As[n], weib.Ks[n])
        #
        # calculate AEP per sector & wtg (using weibull freqs)
        cf = 24*365*1e-9                                                                                               # conversion factor => GWh. 365.24?? TODO
        naep = cf * weib.freqs * nwpw.sum(axis=2)                                                                      # summing over all velocities
        gaep = cf * weib.freqs * gwpw.sum(axis=2)
        #
        if verbose:
            n = naep.sum()
            g = gaep.sum()
            wloss = (g-n)/g * 100
            #
            print(f'{nm} Gross AEP (GWh) {g:.1f}')
            print(f'{nm} Net AEP (GWh) {n:.1f}')
            print(f'{nm} Wake loss (%) {wloss:.2f}')
        #
        # save as DataArray (not SimulationResult - too complex...)
        #
        xa = sim.isel({'wd':range(n_sectors), 'ws':0}).drop('ws')         # just a template with correct dimensions
        nx = xr.DataArray(data=naep.copy().T, dims=xa.dims, coords=xa.coords, attrs=dict(description=f'Net production {nm}'))
        gx = xr.DataArray(data=gaep.copy().T, dims=xa.dims, coords=xa.coords, attrs=dict(description=f'Gross production {nm}'))
        aeps.append([nx, gx])
        if return_pwr:
            nx = xr.DataArray(data=nwpw.mean(axis=2).copy().T, dims=xa.dims, coords=xa.coords, attrs=dict(description=f'Net power {nm}'))
            gx = xr.DataArray(data=gwpw.mean(axis=2).copy().T, dims=xa.dims, coords=xa.coords, attrs=dict(description=f'Gross power {nm}'))
            pwrs.append([nx, gx])
    #
    if return_pwr: return aeps, pwrs
    else         : return aeps

def _write_park_aep(net, gross, park, f):
    # formatting
    fms = {0:lambda x: f'Turbine site {x:.0f}', 1:INT, 2:INT, 3:INT, 4:REAL, 5:REAL, 6:str}
    #
    f.write(f'      Xpos   Ypos   Height   GrAEP   NetAEP\n')
    f.write(f' Site [m]    [m]    [m]      [GWh]   [GWh]\n')
    hub_heights = [park.wtg.hub_height]*park.size
    df = pd.DataFrame([park.ids, park.xs, park.ys, hub_heights, gross, net], dtype=np.float64).T  # needs to be float64
    df[6] = [park.wtg.name]*park.size
    f.write(df.to_string(index=False, header=False, formatters=fms))


def create_output(aeps, case, fnm):
    '''
    create ouput ala Fuga (for knowl)
    only the 'For wake only' option is supported
    note1: prefix 'n' is for net, 'g' is for gross
    '''
    #
    header = f'''py_wake Annual Energy production estimates
z0[m]   99999
zi[m]   99999
zeta0   99999
'''
    f = open(fnm, 'w')
    f.write(header)
    #
    # write sector by sector
    wds = aeps[0][0].wd.values
    for i, wd in enumerate(wds):
        f.write(f'\nSector\t{i}\n')
        for park, (net, gr) in zip(case.park_list, aeps):
            _write_park_aep(net.sel(wd=wd), gr.sel(wd=wd), park, f)
        # write a 'summary' for the total and each park
        ns, gs = [aep[0].sel(wd=wd).sum() for aep in aeps], [aep[1].sel(wd=wd).sum() for aep in aeps]
        f.write(f'\nAllProjects                      {sum(gs).values:.4f}   {sum(ns).values:.4f}\n')
        for park, n, g in zip(case.park_list, ns, gs):
            f.write(f'{park.name:20s}             {g.values:.4f}   {n.values:.4f}\n')
    #
    # write all sectors
    f.write(f'\nAll sectors\n')
    for j, (park, aep) in enumerate(zip(case.park_list, aeps)):
        _write_park_aep(aep[0].sum(axis=1), aep[1].sum(axis=1), park, f)
    #
    # write a 'summary' for the total and each park
    ns, gs = [aep[0].sum(axis=1) for aep in aeps], [aep[1].sum(axis=1) for aep in aeps]
    #
    n, g = sum([sum(x) for x in ns]), sum([sum(x) for x in gs])
    f.write(f'\nAllProjects                      {g.values:.4f}   {n.values:.4f}\n')
    for park, n, g in zip(case.park_list, ns, gs):
        f.write(f'{park.name:20s}             {sum(g).values:.4f}   {sum(n).values:.4f}\n')
    #
    f.close()
    logging.info(f'{fnm} was created')

def read_output_file(fnm):
    '''
    reads a typical output file (from create_output_aep or from Fuga) into an xarray.
    useful for comparing results.
    '''
    lines = open(fnm).readlines()
    tmpfile = 'tmpfile'
    ns = []
    gs = []
    sector = []
    #
    for line in lines:
        #if 'All sectors' in line: break
        if 'AllProjects' in line:
            f = open(tmpfile, 'w')
            f.writelines(sector)
            f.close()
            df = pd.read_csv(tmpfile, sep='\s+', header=None, usecols=range(2,9))
            ns.append(df[7].values)
            gs.append(df[6].values)
            sector = []
            continue
        if 'Turbine site' in line:
            sector.append(line)
            continue
    os.unlink(tmpfile)
    #
    nsecs = len(ns) - 1
    wds = 360/nsecs * np.arange(nsecs+1)
    wds[-1] = -1                     # this is the 'all sectors'
    #
    coords = dict(wd=wds, wt=df[2].values, x=(["wt"], df[3].values), y=(["wt"], df[4].values))
    dims   = ["wd", "wt"]
    #
    net   = xr.DataArray(data=ns, dims=dims, coords=coords, attrs=dict(description=f'Net production'))
    gross = xr.DataArray(data=gs, dims=dims, coords=coords, attrs=dict(description=f'Gross production'))
    return net, gross

def compare_outputs(fnm1, fnm2, lbl1, lbl2, ms=60):
    '''
    reads two output files and compares them (by plotting)
    useful for comparing results from different sources.
    typical usage: compare_outputs('FugaOutput_1.txt', 'pywake1.txt', 'FUGA', 'TurboPark')
    ms: marker-size. set to None for default
    '''
    net1, gr1 = read_output_file(fnm1)
    net2, gr2 = read_output_file(fnm2)
    #
    for wd in net1.wd.values:
        n1, g1 = net1.sel(wd=wd), gr1.sel(wd=wd)
        wl1 = (g1-n1) / g1 *100
        n2, g2 = net2.sel(wd=wd), gr2.sel(wd=wd)
        wl2 = (g2-n2) / g2 *100
        min_val = min(wl1.min(), wl2.min())
        max_val = max(wl1.max(), wl2.max())
        #
        wd_txt = f'wd = {wd:.0f} deg' if wd >= 0 else 'wd = ALL'
        #
        plt.figure(figsize=(18,5))
        #
        plt.subplot(131)
        plt.scatter(wl1.x.values, wl1.y.values, c=wl1.values, s=ms)
        plt.axis('equal')
        cb = plt.colorbar()
        cb.set_label('Wake loss [%]')
        plt.clim(min_val, max_val)
        plt.title(f'{lbl1} ({wd_txt})')
        plt.xticks([]), plt.yticks([])
        #
        plt.subplot(132)
        plt.scatter(wl2.x.values, wl2.y.values, c=wl2.values, s=ms)
        plt.axis('equal')
        cb = plt.colorbar()
        cb.set_label('Wake loss [%]')
        plt.clim(min_val, max_val)
        plt.title(f'{lbl2} ({wd_txt})')
        plt.xticks([]), plt.yticks([])
        #
        wld = wl1 - wl2
        plt.subplot(133)
        plt.scatter(wld.x.values, wld.y.values, c=wld.values, s=ms)
        plt.axis('equal')
        plt.colorbar()
        cb.set_label('Wake loss [%]')
        plt.title('Wake loss difference')
        plt.xticks([]), plt.yticks([])
        plt.tight_layout()
    return net1, gr1, net2, gr2

def create_output_power(power, case, fnm=None):
    unit = power.Description.split('[')[1][:-1]
    assert unit == 'W'
    vels = []
    for k in range(len(power.ws)):
        vel = []
        for i in range(case.size):
            p = power[i,:,k].values / 1000.       # W => kW
            vel.append([case.ids[i]] + list(p))
        vels.append(pd.DataFrame(vel))
    # 
    if fnm:
        fms = {0:lambda x: f'Turbine site {x:d}', 1:REAL, 2:REAL, 3:REAL, 4:REAL, 5:REAL, 6:REAL, 7:REAL, 8:REAL, 9:REAL, 10:REAL, 11:REAL}
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
            f.write(f'\nTotal    ')
            tot = pd.DataFrame(vel.sum()[1:]).T
            f.write(tot.to_string(index=False, header=False, formatters=fms))
        f.write('\n')
        f.close()
        logging.info(f'{fnm} was created')
    #
    return vels

def winddir_components(wind_dir):
    '''
    wind_dir is wind direction, i.e. angle in degrees from north
    '''
    dx = -np.sin(wind_dir/180*np.pi)
    dy = -np.cos(wind_dir/180*np.pi)
    return dx, dy

class Scada(object):
#
    def __init__(self, fnm='', stype='camille'):
        '''
        read scada csv-file (from camille) into a DataFrame
        - input
            fnm: file name (csv-file, could be gzip'ed)
            stype: scada format type. for now only 'camille' is implemented
        - notes
            DataArray less useful for scada-data since wind-speeds and wind-directions are not regular
        '''
        if stype == 'camille':
            if fnm:
                self.data = pd.read_csv(fnm, parse_dates=['time'], converters={'time':dateutil.parser.parse})
            else:
                self.data = pd.DataFrame([])
            self.fnm = fnm
        else:
            raise Exception(f'SCADA type {stype} not implemented')
#
    def select(self, wd_min, wd_max, ws_min, ws_max, nms=[]):
        '''
        selecting part of the data.
        - input
          * wd_min : min wind-dir
          * wd_max : max wind-dir
          * ws_min : min wind-speed
          * ws_max : max wind-speed
          * nms    : wind turbine names (list). if [], all are included
        - returns
          * DataFrame with the selected data. also available as attribute sel
        '''
        sel = self.data              # pointer to all
        ixs = np.logical_and(sel.WindDir>wd_min, sel.WindDir<=wd_max)
        sel = sel[ixs]               # now it's a copy
        ixs = np.logical_and(sel.WindSpeed>ws_min, sel.WindSpeed<=ws_max)
        sel = sel[ixs]
        if nms:
            ixs = [False]*len(sel)   # init
            for nm in nms:
                ixs = np.logical_or(ixs, (sel.Turbine==nm).values)
            sel = sel[ixs]
        self.sel = sel               # potentially useful
        return sel
#
    def net(self, wd_min=-1., wd_max=361., ws_min=-1., ws_max=9999., nms=[]):
        '''
        return net power production for selected conditions.
        use defaults if you want all data
        - input
          * wd_min : min wind-dir
          * wd_max : max wind-dir
          * ws_min : min wind-speed
          * ws_max : max wind-speed
          * nms    : wind turbine names (list). if [], all are included
        - returns
          * Series with the selected data
        '''
        sel = self.select(wd_min, wd_max, ws_min, ws_max, nms)
        return sel.groupby('Turbine').mean().ActivePower
#
    def gross(self, pwr, wd_min=-1., wd_max=361., ws_min=-1., ws_max=9999., nms=[]):
        '''
        return gross power production for selected conditions.
        use defaults if you want all data
        - input
          * pwr    : power-function. assumes one turbine-type only. TODO! 
          * wd_min : min wind-dir
          * wd_max : max wind-dir
          * ws_min : min wind-speed
          * ws_max : max wind-speed
          * nms    : wind turbine names (list). if [], all are included
        - returns
          * Series with the selected data
        '''
        sel = self.select(wd_min, wd_max, ws_min, ws_max, nms)
        sel['gross'] = pwr(sel.WindSpeed)
        return sel.groupby('Turbine').mean().gross


