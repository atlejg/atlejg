import AtlejgTools.Scripts.pywake_for_knowl as pywake_for_knowl
import AtlejgTools.Scripts.run_pywake as run_pywake
import AtlejgTools.WindYieldAssessment.Utils as WU
import AtlejgTools.Utils as UT
import os, logging, sys
import numpy as np
import pandas as pd

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

def pywake_dbc_test(testdir='/project/RCP/active/wind_yield_assessment/agy/WindWorks/Testdata/DBC', yaml_file='1.yml'):
    '''
    a doggerbank test-case
    '''
    cwd = os.getcwd()
    os.chdir(testdir)
    _, _, opts, _, _, _, _ =  run_pywake.main(yaml_file)
    orig = pd.read_csv('orig.csv')
    new  = pd.read_csv(opts.outfile)
    assert np.all(new==orig)
    logging.info(f' testing OK')
    os.chdir(cwd)
