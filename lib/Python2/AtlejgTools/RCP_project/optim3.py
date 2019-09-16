import scipy.optimize
from scipy.optimize import leastsq
import xlrd
from fnmatch import fnmatch

'''
atle j. gyllensten
agy@statoil.com
feb 2014

reads experimental data from Excel spreadsheet as provided from inflow control tests performed in
HP AICD rig in Porsgrunn.
useful for doing curvefitting, plotting etc.

see documentation under load_experimental_data for more details.

nomenclature:
dp : pressure drop
mu  : viscosity
rho : density
vf  : volume fraction
q   : flowrate
kw  : key words
'''


RHO_CAL = 1000.
MU_CAL  = 1.

def _dp_calc(params, mu, rho, q):
   cnst, x, y = params
   return VC.rcp_coeff(rho, mu, RHO_CAL, MU_CAL, cnst, y) * q**x

def _dp_calc2(params, mu, rho, q):
   cnst, x, y, mu0, m = params
   return cnst * rho**m * (mu+mu0)**(-y) * q**x

def _func(params, dp, mu, rho, q):
   return dp - _dp_calc(params, mu, rho, q)

def _get_val2(val, pattern, v_o, v_w, v_g):
   if val == pattern + '_o': return v_o
   if val == pattern + '_w': return v_w
   if val == pattern + '_g': return v_g
   if type(val) is str:     return val     # it's a column letter
   return eval(val)                        # it's a number

class InflowTest(object):
   def __init__(self, fluid_nm='', valve_nm=''):
      self.fluid_nm = fluid_nm
      self.valve_nm = valve_nm
      self.fp    = UT.Struct()  # fluid props
#
   def set_fluid_props(self, mu_o,mu_w,mu_g, rho_o,rho_w,rho_g):
      self.fp.mu_o = mu_o
      self.fp.mu_w = mu_w
      self.fp.mu_g = mu_g
      self.fp.rho_o = rho_o
      self.fp.rho_w = rho_w
      self.fp.rho_g = rho_g
#
   def read(self, column, line_from, line_to): # convinient
      return UT.read_excel_column(self.sheet, column, line_from, line_to, return_array=True)
#
   def get_vals(self, val, line_from, line_to):
      nlines = line_to - line_from + 1
      try:
         return float(val)*ones(nlines)              # it's a hardcoded number
      except: pass
      keys = self.fp.__dict__.keys()
      if val in keys:
         return self.fp.__dict__[val]*ones(nlines)   # it's one of the fluid parameters (rho_o etc)
      column = val
      return self.read(column, line_from, line_to)   # will read values from spreadsheet
#
   def reset(self):
      self.data = {}
#
   def update(self, kw):   # kw = keywords
      if not self.data.has_key(kw.label):
         self.data[kw.label] = {'dp':zeros(0), 'mu':zeros(0), 'rho':zeros(0), 'q':zeros(0)}
      dp   = self.get_vals(kw.dp, kw.line_from, kw.line_to)
      mu1  = self.get_vals(kw.mu1, kw.line_from, kw.line_to)
      rho1 = self.get_vals(kw.rho1, kw.line_from, kw.line_to)
      q1   = self.get_vals(kw.q1, kw.line_from, kw.line_to)
      vf   = self.get_vals(kw.vf, kw.line_from, kw.line_to)
      if str(kw.vf) != '100':
         mu2  = self.get_vals(kw.mu2, kw.line_from, kw.line_to)
         rho2 = self.get_vals(kw.rho2, kw.line_from, kw.line_to)
         q2   = self.get_vals(kw.q2, kw.line_from, kw.line_to)
      else:
         mu2  = 0.
         rho2 = 0.
         q2   = 0.
      vf /= 100. # convert from %
      q1 *= 24.  # convert from m3/h to m3/d
      q2 *= 24.  # convert from m3/h to m3/d
      # append data
      self.data[kw.label]['dp']  = concatenate((self.data[kw.label]['dp'], dp))
      self.data[kw.label]['q']   = concatenate((self.data[kw.label]['q'], q1+q2))
      self.data[kw.label]['mu']  = concatenate((self.data[kw.label]['mu'],  vf*mu1  + (1-vf)*mu2))
      self.data[kw.label]['rho'] = concatenate((self.data[kw.label]['rho'], vf*rho1 + (1-vf)*rho2))
#
   def get(self, *args):
      '''
      returns list of dp, mu, rho, q for given labels
      uses wildcard matching ala unix (i.e. '1-ph *', 'gvf*' etc)
      '''
      dp, mu, rho, q = ([], [], [], [])
      labels = self.data.keys()
      for pattern in args:
         for label in labels:
            if fnmatch(label, pattern):
               print 'getting', label
               dp.extend(self.data[label]['dp'])
               mu.extend(self.data[label]['mu'])
               rho.extend(self.data[label]['rho'])
               q.extend(self.data[label]['q'])
      return (array(dp), array(mu), array(rho), array(q))
      #if not self.data.has_key(label): return None
#
   def get_all(self):
      '''
      returns list of dp, mu, rho, q for all labels
      '''
      dp, mu, rho, q = ([], [], [], [])
      for v in self.data.values():
         dp.extend(v['dp'])
         mu.extend(v['mu'])
         rho.extend(v['rho'])
         q.extend(v['q'])
      return (array(dp), array(mu), array(rho), array(q))
#
   def load_experimental_data(self, inpfile):
      '''
      reads input file to obtain experimental data.
      input file must have the following structure:

# file info
fname            = Copy of NK Tendeka 2.5mm TR7 2013.xlsx
sheet            = testmatrix

# fluid props
rho_o            = 920
rho_w            = 1030
rho_g            = 80.
mu_o             = 115.
mu_w             = 0.5
mu_g             = 0.012

# one phase oil
label            = oil 
start_row        = 29
end_row          = 36
dp               = I
vol_frac         = 100
flowrate1        = J
rho1             = rho_o
mu1              = mu_o

#
# two phase gas/oil
#

label            = gvf1 
start_row        = 157
end_row          = 160
vol_frac         = D
flowrate1        = O
flowrate2        = J
rho1             = rho_g
rho2             = rho_o
mu1              = mu_g
mu2              = mu_o

label            = gvf2 
start_row        = 157
end_row          = 160

      the variables may be a number, a column letter (referring to the spreadsheet) or
      a 'global' value (rho_g, rho_o ... given under fluid props).
      we could have any number of sections where data is read - typically one section
      for oil, one for for gas, one for water, one for 10% gvf etc etc. 
      the volume fraction is given in % (so vol_frac = 100 for single phases)
      each new section is started with 'label'. start_row and end_row must also always be present.
      for each new each section, if something is missing, it uses value from previous section.
      three phase experiments are not handled.

      '''
      #
      # first: read globals
      f = open(inpfile)
      for line in f:
         line = line.strip()
         if line.startswith('fname') : self.fname   = line.split('=')[1].strip(); continue
         if line.startswith('sheet') : self.sheetnm = line.split('=')[1].strip(); continue
         if line.startswith('rho_o') : rho_o = float(line.split('=')[1]);         continue
         if line.startswith('rho_w') : rho_w = float(line.split('=')[1]);         continue
         if line.startswith('rho_g') : rho_g = float(line.split('=')[1]);         continue
         if line.startswith('mu_o' ) : mu_o  = float(line.split('=')[1]);         continue
         if line.startswith('mu_w' ) : mu_w  = float(line.split('=')[1]);         continue
         if line.startswith('mu_g' ) : mu_g  = float(line.split('=')[1]);         continue
         if line.startswith('label' ): break
      self.set_fluid_props(mu_o,mu_w,mu_g, rho_o,rho_w,rho_g)
      self.sheet = xlrd.open_workbook(self.fname).sheet_by_name(self.sheetnm)
      #
      # now read sections
      f.seek(0) # reset
      self.reset()
      kw = UT.Struct()  # keywords
      nsections = 0
      for line in f:
         if 'label' in line:
            if nsections > 0: self.update(kw)
            kw.label = line.split('=')[1].strip()
            nsections += 1
            dp, mu1, mu2, rho1, rho2, q1, q2, vf = [0]*8
            continue
         if 'start_row' in line: kw.line_from = int(line.split('=')[1].strip());   continue
         if 'end_row'   in line: kw.line_to   = int(line.split('=')[1].strip());   continue
         if 'vol_frac'  in line: kw.vf        = line.split('=')[1].strip();        continue
         if 'dp'        in line: kw.dp        = line.split('=')[1].strip();        continue
         if 'flowrate1' in line: kw.q1        = line.split('=')[1].strip();        continue
         if 'flowrate2' in line: kw.q2        = line.split('=')[1].strip();        continue
         if 'rho1'      in line: kw.rho1      = line.split('=')[1].strip();        continue
         if 'rho2'      in line: kw.rho2      = line.split('=')[1].strip();        continue
         if 'mu1'       in line: kw.mu1       = line.split('=')[1].strip();        continue
         if 'mu2'       in line: kw.mu2       = line.split('=')[1].strip();        continue
      f.close()
      self.update(kw) # remember last one
#
   def plot(self, label, params, newfig=True, colour=None):
      if newfig: figure()
      d = self.data[label]
      q = d['q']
      dp_calc =  _dp_calc(params, d['mu'], d['rho'], q)
      q_plot  = q * 1000 / 24.  # m3/d -> l/h
      plot(q_plot, d['dp'], '*')
      plot(q_plot, dp_calc, 'o--')
      xlabel('flow rate [Al/h]')
      ylabel('dp [bar]')
      title(label)
      grid(True)
      legend(('exp', 'model'), loc='best')

it = InflowTest()
excelfile = sys.argv[1]
#fname = 'Copy of NK Tendeka 2.5mm TR7 2013.xlsx'
sheet = xlrd.open_workbook(excelfile).sheet_by_name('metadata')
v1=UT.read_excel_column(sheet, 'A', 1, 1000)
v2=UT.read_excel_column(sheet, 'B', 1, 1000)

it.load_experimental_data(sys.argv[1])
#rcp_params0 = (1e-3, 3, .6, 1., 2)

rcp_params0 = (1e-4, 3, .6)  # initial guess
dp,mu,rho,q = it.get_all()

rcp_params1 = leastsq(_func, rcp_params0, args=(dp,mu,rho,q))[0]
dp_calc1 = _dp_calc(rcp_params1, mu, rho, q)
rcp_params2 = leastsq(_func, rcp_params0, args=it.get('oil', 'water', 'gas'))[0]
dp_calc2 = _dp_calc(rcp_params2, mu, rho, q)
rcp_params3 = leastsq(_func, rcp_params0, args=it.get('gvf*'))[0]
dp_calc3 = _dp_calc(rcp_params3, mu, rho, q)

if True:
   figure()
   plot(q, dp, 'r*')
   plot(q, dp_calc1, 'k.')
   plot(q, dp_calc2, 'g.')
   plot(q, dp_calc3, 'b.')
   legend(('exp', 'all', 'singlephase', 'gvf'), loc='best')

if False:
   for label in it.data.keys(): it.plot(label, rcp_params1)
   UT.hardcopies(12, 'Figs', 'NK_all_')
   close('all')
   for label in it.data.keys(): it.plot(label, rcp_params2)
   UT.hardcopies(12, 'Figs', 'NK_singlephase_')
   close('all')
   for label in it.data.keys(): it.plot(label, rcp_params3)
   UT.hardcopies(12, 'Figs', 'NK_gvf_')
   close('all')

show()
