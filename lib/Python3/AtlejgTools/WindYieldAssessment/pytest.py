import matplotlib.pyplot as plt
import AtlejgTools.Scripts.pywake_for_knowl as pywake_for_knowl
import AtlejgTools.Scripts.run_pywake as rp
import AtlejgTools.WindYieldAssessment.Utils as WU
import AtlejgTools.Utils as UT
import os, logging, sys
import numpy as np
import pandas as pd
import xarray as xr

ACCURACY = 1.1e-4

logging.basicConfig(level=logging.INFO)

def _ae(res1, res2, accur=ACCURACY):
    '''
    ae: almost equal
    we are comparing content of two files with values on format %.4f
    we could (should?) accept difference of last decimal place...
    '''
    err = abs(res2 - res1).values.max()
    logging.info(f' max error: {err:.2e}')
    return err <= accur

def zero_results(nwt, nwd, nws):
    '''
    create a zero-valued DataArray in a convinent way
    '''
    dims   = ["wt", "wd", "ws"]
    coords = dict(wt=range(nwt),
                  wd=np.linspace(0, 360, nwd, endpoint=False),
                  ws=np.linspace(4, 15, nws, endpoint=True),
                 )
    zeros = np.zeros((nwt, nwd, nws))
    return xr.DataArray(data=zeros, dims=dims, coords=coords)

def test_shifting(width=30.):
    '''
    when doing averaging into sectors, we usually want the first sector
    to be centered around wd=0 (north). this means we need to include results
    from wd<0 (i.e. 345-359 deg). for this purpose we need to shift all the
    data half a sector-width to the right.
    '''
    pwr = zero_results(1, 360, 1).sel(wt=0, ws=4.)
    pwr.values = 1e6*np.cos(pwr.wd/180.*np.pi) + 2e6
    #
    dwd = pwr.wd.values[1] - pwr.wd.values[0]
    #
    pwr_s = WU.shift_it(pwr, width)
    #
    plt.figure()
    plt.plot(pwr.wd, pwr.values, 'k-', label='Not shifted')
    plt.plot(pwr_s.wd, pwr_s.values, 'r-', label='Shifted')
    #
    n_sectors = int(360/width)
    n_bins    = int(width/dwd)
    pwrc1 = pwr.coarsen(wd=n_bins).mean()
    pwrc2 = WU.coarsen(pwr, n_sectors)
    plt.plot(pwrc1.wd, pwrc1.values, 'ko', ms=12, label='Averaged')
    plt.plot(pwrc2.wd, pwrc2.values, 'ro', ms=12, label='Averaged')
    plt.legend(loc='best')
    plt.show()


def test_synt_scada(n=1000, ws_max=15., cleanup=True):
    '''
    create syntetic scada data, write to file, read from file,
    plot heatmap, and exit
    the heatmap should look like some exotic stairs.
    - input
      * n      : number of test-data entries
      * ws_max : max wind-speed
      * cleanup: remove file after use
    '''
    fnm = 'synt_scada.csv'
    wd = np.linspace(0,(n-1)/2, n) % 360.
    ws = ws_max*wd/360
    data = pd.DataFrame(np.array([wd, ws]).T, columns=['WindDir', 'WindSpeed'])
    data['Turbine']     = 'WT-1'
    data['ActivePower'] = 1.
    data['time']        = '2018-01-01 00:00:00+00:00'
    data.to_csv(fnm, sep=',')
    synt = WU.Scada(fnm)
    synt.heatmap(cbar=True, cticks=True)
    if cleanup:
        os.unlink(fnm)
    return synt

