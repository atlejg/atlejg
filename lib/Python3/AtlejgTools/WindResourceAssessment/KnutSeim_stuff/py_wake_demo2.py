#!/usr/bin/env python
# coding: utf-8

'''
this is a scriptified version of knut s. seim's notebook.
'''

import py_wake
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


import Dudgeon2

wake_model = sys.argv[1]         # 'Fuga', 'TP', or 'NOJ'
n_old      = int(sys.argv[2])    #  number of "old" WTGs

wt_x, wt_y = Dudgeon2.wt_x, Dudgeon2.wt_y
types = np.zeros(len(wt_x)); types[n_old:] = 1

site = Dudgeon2.DudgeonSite()
wtgs = Dudgeon2.TwoWTGs(n_old, len(wt_x)-n_old, Dudgeon2.D, Dudgeon2.D, Dudgeon2.hh, Dudgeon2.hh)

plot_layout = False
plot_wind   = True
do_more     = False

# ### Layout

if plot_layout:
    fig,ax=plt.subplots(figsize=(10,6))
    wtgs.plot(wt_x, wt_y, ax=ax)


# ### Site wind rose

if plot_wind:
    fig,ax=plt.subplots(figsize=(10,6))
    site.plot_wd_distribution(n_wd=12, ws_bins=[0,5,10,15,20,25], ax=ax)


# ### Wake model

# Wake model input
#
A = 0.6 # Only used if wake_model='TP', A parameter in dDw/dx = A*I(x)
k = 0.04 # Only used if wake_model='NOJ', wake expansion parameter
dwd = 0.5  # resolution in wind direction
dws = 0.25 # resolution in wind speed


if wake_model.upper()   =='FUGA':
    lut_path = 'FugaLUTs/Z0=0.00012000Zi=00790Zeta0=2.00E-7/'
    wf_model = py_wake.Fuga(lut_path, site, wtgs)
elif wake_model.upper() =='TP':
    wf_model = py_wake.TP(site, wtgs, k=k)
elif wake_model.upper() =='NOJ':
    wf_model = py_wake.NOJ(site, wtgs, k=A)
else:
    print('The Fuga, TP, or the NOJ wake models are the only options in this notebook. For other wake models, see the PyWake documentation.')


# Run wake model for all combinations of wd and ws
sim_res = wf_model(wt_x, wt_y, type=types, wd=np.arange(0,360,dwd), ws=np.arange(3,30+dws,dws)) 

gross = sim_res.aep(with_wake_loss=False).values.sum()
net = sim_res.aep().values.sum()
wl = ((sim_res.aep(with_wake_loss=False).sum()-
       sim_res.aep().sum())/sim_res.aep(with_wake_loss=False).sum()).values*100


pd.DataFrame({'Wake model: {}'.format(wake_model): np.round([gross, net, wl], 1)},
             index=['Gross AEP (GWh)', 'Net AEP (GWh)', 'Wake loss (%)'])


# ### Flow map

flow_map = sim_res.flow_map(grid=None, # defaults to HorizontalGrid(resolution=500, extend=0.2)
                            wd=230,
                            ws=10.)


fig,ax = plt.subplots(figsize=(8,8))
#flow_map.plot_wake_map(levels=linspace(4, 10, 60), ax=ax);
flow_map.plot_wake_map(levels=np.arange(2,12), ax=ax);
#flow_map.plot_wake_map(ax=ax);
ax.set_title(f'Wake model = {wake_model}. n_old = {n_old}')



if do_more:
    # ### Weibull weighted AEP

    def binArray(data, axis, binstep, binsize, func=np.nanmean):
        data = np.array(data)
        dims = np.array(data.shape)
        argdims = np.arange(data.ndim)
        argdims[0], argdims[axis]= argdims[axis], argdims[0]
        data = data.transpose(argdims)
        data = [func(np.take(data,np.arange(int(i*binstep),int(i*binstep+binsize)),0),0) for i in np.arange(dims[axis]//binstep)]
        data = np.array(data).transpose(argdims)
        return data

    def WeightedPower(u, power, A, k):
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
        def G(alpha, k):
            # 1/k times the incomplete gamma function of the two arguments 1/k and alpha^k
            # Note, the scipy incomplete gamma function, gammainc, must be multiplied with gamma(k) to match the
            # the G function used in the European Wind Atlas
            import scipy.special as sc
            return 1 / k * sc.gamma(1 / k) * sc.gammainc(1 / k, alpha**k)

        u0, u1 = u[:-1], u[1:]
        alpha0, alpha1 = (u0 / A), (u1 / A)
        alpha0, alpha1 = (u0 / A)+np.random.rand()*1e-10, (u1 / A)+np.random.rand()*1e-10
        p0, p1 = power[..., :-1], power[..., 1:]

        res = (p0 * np.exp(-alpha0**k) +  # eq 6.5, p0 * cdf(u0:)
               (p1 - p0) / (alpha1 - alpha0) * (G(alpha1, k) - G(alpha0, k)) -  # eq 6.4 linear change p0 to p1
               p1 * np.exp(-alpha1**k)
               )  # eq 6.5, - p1 * cdf(u1:)

        return res


    gpw = np.tile(wf_model.windTurbines.power(sim_res.ws),(sim_res.Power.shape[0],sim_res.Power.shape[1],1))


    tmp=binArray(sim_res.Power.values,1,60,60)
    gtmp=binArray(gpw,1,60,60)


    a=site.ds.Weibull_A[0:-1].values
    k=site.ds.Weibull_k[0:-1].values
    f=site.ds.Sector_frequency[0:-1].values


    wagpw=np.ndarray((gtmp.shape[0],gtmp.shape[1],gtmp.shape[2]-1))
    wapw=np.ndarray((tmp.shape[0],tmp.shape[1],tmp.shape[2]-1))
    for n in range(12):
        wagpw[:,n,:]=WeightedPower(sim_res.ws.values, gtmp[:,n,:], a[n], k[n])
        wapw[:,n,:]=WeightedPower(sim_res.ws.values, tmp[:,n,:], a[n], k[n])


    pw=0
    gpw=0
    for n in range(wapw.shape[1]):
        gpw+=f[n]*wagpw[:,n,:].sum()
        pw+=f[n]*wapw[:,n,:].sum()


    wwgross = gpw*24*365*1e-9
    wwnet = pw*24*365*1e-9
    wwwl = ((wwgross-wwnet)/wwgross)*100


    pd.DataFrame({'Wake model: {}'.format(wake_model): np.round([wwgross, wwnet, wwwl], 1)},
                 index=['Gross AEP (GWh)', 'Net AEP (GWh)', 'Wake loss (%)'])


show()

