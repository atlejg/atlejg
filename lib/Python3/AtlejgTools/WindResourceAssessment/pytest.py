import AtlejgTools.Scripts.pywake_for_knowl as pywake_for_knowl
import AtlejgTools.Scripts.run_pywake as run_pywake
import AtlejgTools.WindResourceAssessment.Utils as WU
import AtlejgTools.Utils as UT
import os, logging, sys
import numpy as np
import pandas as pd

ACCURACY = 1.1e-4

if 'linux' in sys.platform.lower():
    KNOWLDIR  = r'/project/RCP/active/wind_resource_assessment/agy/Knowl/Testdata'
else:
    KNOWLDIR  = r'D:\OneDrive - Equinor\ACTIVE\Wind\knowl\Testdata'

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

def _knowl_cases(pattern, noj_only, accur=ACCURACY):
    cwd = os.getcwd()
    os.chdir(KNOWLDIR)
    case_dirs = UT.glob(pattern); case_dirs.sort()
    #
    for case_dir in case_dirs:
        logging.info('')
        logging.info(f'------ case_dir : {case_dir} ------')
        os.chdir(case_dir)
        #
        logging.info('------ NOJ ------')
        pywake_for_knowl.main('NOJ', dump_results=False)
        r1 = WU.read_output_file('NOJ.txt')           # existing file
        r2 = WU.read_output_file('FugaOutput_1.txt')  # created now
        #
        assert _ae(r1.net, r2.net, accur)
        assert _ae(r1.gross, r2.gross, accur)
        #
        if not noj_only:
            #
            logging.info(' ------ TurbOPark ------')
            pywake_for_knowl.main('ETP', dump_results=False)
            r1 = WU.read_output_file('ETP.txt')            # existing file
            r2 = WU.read_output_file('FugaOutput_1.txt')  # created now
            #
            assert _ae(r1.net, r2.net, accur)
            assert _ae(r1.gross, r2.gross, accur)
        os.chdir(KNOWLDIR)
    os.chdir(cwd)
    logging.info(f' testing OK')

def knowl_large_cases(noj_only=False):
    _knowl_cases('Phase?_*/', noj_only)

def knowl_small_case(noj_only=False):
    _knowl_cases('WestermostRough', noj_only)

def pywake_dbc_test(testdir='/project/RCP/active/wind_resource_assessment/WindFlow/Testdata/DBC', yaml_file='1.yml'):
    '''
    a doggerbank test-case
    '''
    cwd = os.getcwd()
    os.chdir(testdir)
    _, _, opts, _, _, _, _ =  run_pywake.main(yaml_file, dump_results=False)
    orig = pd.read_csv('orig.csv')
    new  = pd.read_csv(opts.outfile)
    assert all(new==orig)
    logging.info(f' testing OK')
    os.chdir(cwd)

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


