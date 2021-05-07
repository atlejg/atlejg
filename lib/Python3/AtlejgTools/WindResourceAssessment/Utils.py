import numpy as np

def _binArray(data, axis, binstep, binsize, func=np.nanmean):
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis]= argdims[axis], argdims[0]
    data = data.transpose(argdims)
    data = [func(np.take(data,np.arange(int(i*binstep),int(i*binstep+binsize)),0),0) for i in np.arange(dims[axis]//binstep)]
    data = np.array(data).transpose(argdims)
    return data

def _g_func(alpha, k):
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
    alpha0, alpha1 = (u0 / wb_a), (u1 / wb_a)
    alpha0, alpha1 = (u0 / wb_a)+np.random.rand()*1e-10, (u1 / wb_a)+np.random.rand()*1e-10
    p0, p1 = power[:, :-1], power[:, 1:]
    #
    wpw = (    p0 * np.exp(-alpha0**wb_k)                                                                          # eq 6.5, p0 * cdf(u0:)
             + (p1 - p0) / (alpha1 - alpha0) * (_g_func(alpha1, wb_k) - _g_func(alpha0, wb_k))                     # eq 6.4 linear change p0 to p1
             - p1 * np.exp(-alpha1**wb_k)                                                                          # eq 6.5, - p1 * cdf(u1:)
    )
    #
    return wpw

def calc_AEP_old(sim, pwr_func, weib, dwd=1., verbose=False):
    '''
    calculates averaged AEP for each sector.
    # 
    based on scriptifie version of knut s. seim's notebook:
    /cygdrive/c/Appl/atlejg/lib/Python3/AtlejgTools/WindResourceAssessment/KnutSeim_stuff/py_wake_demo2.py
    #
    notes: 
        n1: only single WTG supported
        n2: only single weibull supported
        n3: only regular sectors supported
    # 
    input:
        sim      : simulation result from pywake
        pwr_func : power function for WTG (see n1)
        weib     : weibull per sector. must have attributes weib.Ks, weib.As, weib.freqs (see n2)
                   note that sector width is given by number of entries in weibull (see n3)
        dwd      : delta wd. [deg]. gives number of bins.
    output:
        net & gross AEP in GWh per WTG & sector
    '''
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
    gwpw = np.ndarray((len(sim.wt), n_sectors, len(sim.ws)-1))
    nwpw = np.ndarray((len(sim.wt), n_sectors, len(sim.ws)-1))
    for n in range(n_sectors):
        gwpw[:,n,:] = _weigthed_pwr(sim.ws.values, gpwrb[:,n,:], weib.As[n], weib.Ks[n])
        nwpw[:,n,:] = _weigthed_pwr(sim.ws.values, npwrb[:,n,:], weib.As[n], weib.Ks[n])
    #
    # calculate AEP per sector & wtg (using weibull freqs)
    cf = 24*365*1e-9                                                                                               # conversion factor => GWh. 365.24?? TODO
    naep = cf * weib.freqs * nwpw.sum(axis=2)                                                                      # summing over all velocities
    gaep = cf * weib.freqs * gwpw.sum(axis=2)
    #
    if verbose:
        g = gaep.sum()
        n = naep.sum()
        wloss = (g-n)/g * 100
        #
        print(f'Gross AEP (GWh) {g:.1f}\nNet AEP (GWh) {n:.1f}\nWake loss (%) {wloss:.2f}\n')
    #
    return gaep, naep

def calc_AEP(sim0, pwr_funcs, weibs, dwd=1., verbose=False):
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
    #
    for i, (pwr_func, weib) in enumerate(zip(pwr_funcs, weibs)):
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
            print(f'Gross AEP (GWh) {g:.1f}\nNet AEP (GWh) {n:.1f}\nWake loss (%) {wloss:.2f}\n')
        #
        aeps.append([naep, gaep])
    return aeps

def winddir_components(wind_dir):
    '''
    wind_dir is wind direction, i.e. angle in degrees from north
    '''
    dx = -np.sin(wind_dir/180*np.pi)
    dy = -np.cos(wind_dir/180*np.pi)
    return dx, dy

