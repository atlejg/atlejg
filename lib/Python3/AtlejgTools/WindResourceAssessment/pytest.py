import AtlejgTools.WindResourceAssessment.pywake_for_knowl as pywake_for_knowl
import AtlejgTools.WindResourceAssessment.Utils as WU
import AtlejgTools.Utils as UT
import os, logging
import numpy as np

logging.basicConfig(level=logging.INFO)

testdir = 'D:\\OneDrive - Equinor\\ACTIVE\\Wind\\knowl\\Testdata'


def _test_cases(pattern):
    cwd = os.getcwd()
    os.chdir(testdir)
    case_dirs = UT.glob(pattern); case_dirs.sort()
    #
    for case_dir in case_dirs:
        logging.info(f'case_dir= {case_dir}')
        os.chdir(case_dir)
        #
        pywake_for_knowl.main('NOJ')
        net1, gr1 = WU.read_output_file('NOJ.txt')
        net2, gr2 = WU.read_output_file('FugaOutput_1.txt')
        #
        assert np.all(net1 == net2)
        #
        pywake_for_knowl.main('TurbOPark')
        net1, gr1 = WU.read_output_file('TP.txt')
        net2, gr2 = WU.read_output_file('FugaOutput_1.txt')
        #
        assert np.all(net1 == net2)
    os.chdir(cwd)
    logging.info(f'testing OK')

def test_large_cases():
    _test_cases('Phase?_*/')

def test_small_case():
    _test_cases('WestermostRough')

