import AtlejgTools.Utils as UT
import xlrd
import pylab as py

'''
OBSOLETE: Replaced by PICT.py

pict: passive inflow control technology

this was made to plot performance curves for RCP/ICD valves tested in porsgrunn.

atle j. gyllensten
agy@statoil.com
2011
'''

class PictExperiment(object):
    '''
    tested only for heavy oil tests: water, gas, oil, oil+water (not oil + gas etc)
    '''

    def __init__(self, fname, sheetnm, dpt_col, dpa_col, orat_col, wrat_col, grat_col, start_row, end_row, valve, fluid, rho_col=None):
        '''
        reads test data from excel sheet. column names must be supplied.
        dpt_col: dp_target column name/letter
        dpa_col: dp_actual column name/letter
        ...
        '''
        self.sheet = xlrd.open_workbook(fname).sheet_by_name(sheetnm)
        dp_t_ = UT.read_excel_column(self.sheet, dpt_col, start_row, end_row)
        dp_a_ = UT.read_excel_column(self.sheet, dpa_col, start_row, end_row)
        orat_ = UT.read_excel_column(self.sheet, orat_col, start_row, end_row)
        grat_ = UT.read_excel_column(self.sheet, grat_col, start_row, end_row)
        wrat_ = UT.read_excel_column(self.sheet, wrat_col, start_row, end_row)
        if rho_col:
            rho_ = UT.read_excel_column(self.sheet, rho_col, start_row, end_row)
        grat = []; orat = []; wrat = []; dp_t = []; dp_a = []; rho = []
        # need to handle cases where values are missing
        i = -1
        for i in range(len(dp_a_)):
            i += 1
            try:
                dp_a.append(float(dp_a_[i]))
                try:    dp_t.append(float(dp_t_[i]))
                except: dp_t.append(0.)
                try:    grat.append(float(grat_[i]))
                except: grat.append(0.)
                try:    orat.append(float(orat_[i]))
                except: orat.append(0.)
                try:    wrat.append(float(wrat_[i]))
                except: wrat.append(0.)
                if rho_col:
                    try:    rho.append(float(rho_[i]))
                    except: rho.append(0.)
            except: pass
        # now we should be fine
        self.dp_a  = py.array(dp_a)
        self.dp_t  = py.array(dp_t)
        self.orat  = py.array(orat)
        self.grat  = py.array(grat)
        self.wrat  = py.array(wrat)
        self.lrat  = self.orat + self.wrat
        self.flowr = self.lrat + self.grat
        if rho_col:
            self.rho  = py.array(rho)
        self.wc    = self.wrat / self.lrat
        self.gvf   = self.grat / self.flowr
        self.valve = valve
        self.fluid = fluid
        self.fname = fname

    def get(self, gvf=0, wc=1, up_only=True, down_only=False):
        if down_only: up_only = False # down_only overrides up_only
        indx1 = (abs(self.gvf - gvf) <= gvf/5.) # find where gvf is close to the requested value
        indx2 = (abs(self.wc - wc) <= wc/5.)    # find where wc is close to the requested value
        indx = indx1 * indx2
        q = self.flowr[indx]
        dp= self.dp_a[indx]
        if up_only:
            m = (dp == max(dp)).nonzero()[0][0]
            dp = dp[0:m+1]
            q  = q[0:m+1]
        if down_only:
            m = (dp == max(dp)).nonzero()[0][0]
            dp = dp[m:]
            q  = q[m:]
            # remove repeaters also..
            m = (dp == min(dp)).nonzero()[0][0]
            dp = dp[0:m+1]
            q  = q[0:m+1]
        return (q, dp)

    def plot(self, wc_list, gvf=0, up_only=True, down_only=False, **kwargs):
        py.figure()
        for wc in wc_list:
            q, dp = self.get(gvf=gvf, wc=wc, up_only=up_only, down_only=down_only)
            lbl = '%s wc=%.2f' % (self.valve, wc)
            py.plot(q, dp, '-*', label=lbl, **kwargs)
        py.xlabel('flowrate [m3/h]')
        py.ylabel('dp [bar]')
        py.legend(loc='best')
        py.title('%s. GVF = %.2f' % (self.fluid, gvf))
        py.grid(True)
        py.show()

if __name__ == '__main__':

    inp = UT.InputValues()
    inp.add_variable('start_row', 9)
    inp.add_variable('end_row', 100)
    inp.add_variable('target_dp_column', 'B')
    inp.add_variable('actual_dp_column', 'G')
    inp.add_variable('oil_rate_column', 'H')
    inp.add_variable('wat_rate_column', 'O')
    inp.add_variable('gas_rate_column', 'K')
    inp.add_variable('fluid', 'Peregrino')
    inp.add_variable('valve', 'AR7')
    inp.add_variable('pattern', 'P-S7-1i')

    if len(sys.argv) == 3:
        inp.create_input_file(sys.argv[1])
        sys.exit(0)

    inp.read_input_file(sys.argv[1])

    fname = UT.glob('*%s*.xls'%inp.pattern)[0]

    r = PictExperiment(
           fname,
           inp.pattern+'-#',
           inp.target_dp_column,
           inp.actual_dp_column,
           inp.oil_rate_column,
           inp.wat_rate_column,
           inp.gas_rate_column,
           inp.start_row,
           inp.end_row,
           inp.valve,
           inp.fluid)

    r.plot(gvf=0., wc_list=[1.])
