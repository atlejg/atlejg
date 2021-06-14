import AtlejgTools.WindResourceAssessment.pywake_for_knowl as pywake_for_knowl
import AtlejgTools.WindResourceAssessment.run_pywake as run_pywake
import AtlejgTools.WindResourceAssessment.Utils as WU
import AtlejgTools.Utils as UT
import os, logging
import numpy as np
import pandas as pd

ACCURACY = 1.1e-4
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
        logging.info(f'----- case_dir : {case_dir} ----- ')
        os.chdir(case_dir)
        #
        logging.info(' NOJ')
        pywake_for_knowl.main('NOJ')
        net1, gr1 = WU.read_output_file('NOJ.txt')
        net2, gr2 = WU.read_output_file('FugaOutput_1.txt')
        #
        assert _ae(net1, net2, accur)
        assert _ae(gr1, gr2, accur)
        #
        if not noj_only:
            #
            logging.info(' TurbOPark')
            pywake_for_knowl.main('TurbOPark')
            net1, gr1 = WU.read_output_file('TP.txt')
            net2, gr2 = WU.read_output_file('FugaOutput_1.txt')
            #
            assert _ae(net1, net2, accur)
            assert _ae(gr1, gr2, accur)
        os.chdir(KNOWLDIR)
    os.chdir(cwd)
    logging.info(f' testing OK')

def knowl_large_cases(noj_only=False):
    _knowl_cases('Phase?_*/', noj_only)

def knowl_small_case(noj_only=False):
    _knowl_cases('WestermostRough', noj_only)

def pywake_dbc_test(testdir='/project/RCP/active/wind_resource_assessment/WindFlow/Testdata/DBC', yaml_file='1.yml'):
    cwd = os.getcwd()
    os.chdir(testdir)
    _, _, opts, _, _, _, _ =  run_pywake.main(yaml_file)
    orig = pd.read_csv('orig.csv')
    new  = pd.read_csv(opts.outfile)
    assert all(new==orig)
    logging.info(f' testing OK')
    os.chdir(cwd)


