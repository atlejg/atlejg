'''
New version for experiments in Porsgrunn rig (replaces PictExperiment.py)
Nov 2015
Atle J. Gyllensten
agy@statoil.com

Nomenclature:
d = density
v = viscosity
q = volume flow
m = mass flow
dp = pressure
o = oil
w = water
g = gas
l = liquid (oil + water)
t = total  (oil + water + gas)

The valve opening (h_div_h0) and the reynolds number (re) is based on Reidar Schuller's report
"Final Report - New AICD model based on experimental data.docx" from Oct 2015.
'''

import AtlejgTools.RCP_project.PictExperiment        as PE
import AtlejgTools.RCP_project.Valve_characteristics as VC
import AtlejgTools.SimulationTools.UnitConversion    as U
import AtlejgTools.Utils                             as UT
import sys, os, xlrd, re
from pylab import *

def read_testmatrix(fname, sheetnm='testmatrix'):
   return xlrd.open_workbook(fname).sheet_by_name(sheetnm)

# globals. useful to cache data
CM = UT.CacheManager(read_testmatrix, verbose=True)
DATA = {}

def read_data(c, ps):
   '''
   c is a struct with 'constants'
    - column-names for dp, qo etc.
   ps is a parameter-set (a struct) with parameters for one data-serie
    - incl_pattern: a reg-exp pattern for selecting data (typically OW, OG etc)
      note: must match from first character
    - excl_pattern: a pattern for de-selecting data (simple python-mathcing)
    - indx1 and indx2: indices for selecting data (indx1 <= indx <= indx2)
      indx is 13 for test OW-13
    - excl_rowno : indices for de-selecting data (this is spreadsheet row-number)
   '''
   sheet = CM.get(ps.fname)
   s = UT.Struct()
   # read data
   if not ps.quickcurve:
      tags1 = ['p','T','dp','qo','qg','qw','do','dw']
      tags2 = tags1 + ['vo','vg','vw']
      names = UT.read_excel_column(sheet, c.nm, 1, sheet.nrows)
      s.nm = []
      for tag in tags2 + ['vo','vg','vw']: s.__dict__[tag] = []
      for i in range(sheet.nrows):
         if ps.excl_pattern in names[i]: continue
         if i+1 in ps.excl_rowno       : continue
         t2 = UT.read_excel_column(sheet, c.t2, i+1, i+1)[0]
         m = re.match('%s(\d+)'%ps.incl_pattern, names[i]) # note: must match from first character
         if not (m and t2): continue
         indx = int(m.group(1))
         if not ps.indx1 <= indx <= ps.indx2: continue
         s.nm.append(names[i])
         for tag in tags1:
            s.__dict__[tag].append(UT.read_excel_column(sheet, c.__dict__[tag], i+1, i+1)[0])
         print names[i], s.dp[-1]
      # clean it
      for tag in tags1:
         for i in range(len(s.nm)):
            if not s.__dict__[tag][i]: s.__dict__[tag][i] = 0
            else                     : s.__dict__[tag][i] = float(s.__dict__[tag][i])
      # convert to arrays
      for tag in tags1: s.__dict__[tag] = array(s.__dict__[tag]) 
      # derived data
      s.dg = ones(s.dp.shape)*UT.read_excel_column(sheet, c.dg[0],c.dg[1],c.dg[1])[0]
      s.vo = c.vo_func(s.p, s.T)
      s.vw = c.vw_func(s.p, s.T)
      s.vg = c.vg_func(s.p, s.T)
      s.ql = s.qo + s.qw
      s.qt = s.ql + s.qg
      if sum(s.ql) > 0:
         s.wc = s.qw/(s.ql)
      else:
         s.wc = zeros(s.dp.shape)
      s.gvf = s.qg/s.qt * 100   # in percent
      s.mo = s.qo * s.do
      s.mw = s.qw * s.dw
      s.mg = s.qg * s.dg
      s.ml = s.mo + s.mw
      s.mt = s.ml + s.mg
      s.wcp = 100 * s.wc  # wc in percent
      # Reidars model.
      # - mix properties
      s.dl = s.do - (s.do-s.dw)*s.wc
      s.dt = s.dl - (s.dl-s.dg)*s.gvf
      s.vl = s.vo - (s.vo-s.vw)*s.wc**c.v_exp_o
      s.vt = s.vl - (s.vl-s.vg)*s.gvf**c.v_exp_g
      # - reynolds-number and dimensionless numbers
      s.re = s.mt/3600./s.vt/pi/c.radius
      s.y1  = s.mt/3600./pi/2./c.radius/c.h0/(2*s.dt*s.dp*1e5)**0.5   # this is Cd*h/h0 - see eq.6, p.9
      s.h_div_h0  = s.y1/c.Cd                                         # this is h/h0
   else:
      # handle quickcurve. very limited data available (for now...)
      sh = sheet.book.sheet_by_name(ps.sheetnm)
      dp = array(UT.read_excel_column(sh, ps.dp_column, ps.indx1, ps.indx2))[::ps.decimate]
      qt = array(UT.read_excel_column(sh, ps.q_column, ps.indx1, ps.indx2))[::ps.decimate]
      if ps.npoints > 1:
         s.qt = linspace(0, max(qt), ps.npoints)
         s.dp = interp(s.qt, qt, dp)
      else:
         s.qt = qt
         s.dp = dp
   # convinient
   s.npoints = len(s.dp)
   s.sheet = sheet
   s.ps = ps       # for debugging etc.
   DATA[ps.label] = s
   return s

def read_data2(fnm1, fnm2):
   '''
   just a 'wrapper' for easy access.
   fnm1: name of input file for constants
   fnm2: name of input file for datasets
   '''
   c = UT.InputValues(fnm1)
   psets = get_parametersets(fnm2)
   data = []
   for ps in psets:
      data.append(read_data(c, ps))
   return data

def get_data(ident):
   '''
   easy access to cached data.
   if ident is an integer, it will return the data found in this position (sorted on labels,
   which is the default when you show the content of DATA)
   else, it will assume ident is a label to look for
   '''
   if type(ident) == int:
      keys = DATA.keys()
      keys.sort()
      return DATA[keys[ident]]
   if not ident in DATA.keys():
      print 'No such label found:', ident
      return None
   return DATA[ident]

def get_datas(valveconfig, ident, excl_phrase=None):
   '''
   easy access to cached data.
   valveconfig is typically 1v, overunder etc.
   ident is some phrase to look for - like GVF etc.
   excl_phrase is some phrase we *dont* want
   '''
   res = []
   keys = DATA.keys();
   keys.sort()
   for key in keys:
      includeit = False
      if valveconfig in key and ident in key: includeit = True
      if excl_phrase and excl_phrase in key : includeit = False
      if includeit: res.append(DATA[key])
   return res

def _set_defaults(ps):
   # some initial values
   if not 'xscaler'    in ps.__dict__.keys(): ps.xscaler    = 1.
   if not 'yscaler'    in ps.__dict__.keys(): ps.yscaler    = 1.
   if not 'unit_q'     in ps.__dict__.keys(): ps.unit_q     = 'm3/d'  # flow rate
   if not 'unit_dp'    in ps.__dict__.keys(): ps.unit_dp    = 'bar'   # pressure
   if not 'unit_v'     in ps.__dict__.keys(): ps.unit_v     = 'Pa*s'  # viscosity
   if not 'annotate'   in ps.__dict__.keys(): ps.annotate   = None
   if not 'active'     in ps.__dict__.keys(): ps.active     = True
   if not 'quickcurve' in ps.__dict__.keys(): ps.quickcurve = False
   if ps.quickcurve and not 'decimate' in ps.__dict__.keys(): ps.decimate   = 1

def get_parametersets(fname, tmpfile='t'):
   lines = open(fname).readlines()
   f = open(tmpfile, 'w')
   nsets = 0
   psets = []
   for line in lines:
      if 'fname' in line:
         nsets += 1
         if nsets > 1:
            f.close()
            psets.append(UT.InputValues(f.name))
            f = open(tmpfile, 'w')
      f.write(line)
   f.close()
   psets.append(UT.InputValues(f.name))
   os.unlink(tmpfile)
   for ps in psets: _set_defaults(ps)
   return psets

def reidars_model(re):
   y = []
   for re_ in re:
      if re_ < 200:
         y.append(0.13*re_**0.38)
      elif 200 <= re_ < 3000:
         y.append(300*re_**-1)
      else:
         y.append(1.1*re_**-0.2)
   return array(y)

def plotit(xvar, yvar, s, c, ps, figno=None):
   x = s.__dict__[xvar] * ps.xscaler
   y = s.__dict__[yvar] * ps.yscaler
   if figno: figure(figno)
   plot_ = loglog if c.logplot else plot
   plot_(x, y, ps.fmt, label=ps.label, ms=c.markersize)
   if ps.annotate:
      for i in range(s.npoints):
         v = s.__dict__[ps.annotate][i]
         text(x[i],y[i], '%i'%v)

def wsegaicd_dp(par, s):
   '''
   par is WSEGAICD params : (rho_cal, mu_cal, cnst, x, y)
   s is a PICT-struct (should i make it a class now ???)
   '''
   dp = zeros(s.npoints)
   cfq = 1000. if s.ps.unit_q == 'm3/d' else 1.   # conversion factor flow rate
   cfv = 1000. if s.ps.unit_v == 'Pa*s' else 1.   # conversion factor viscosity
   for i in range(s.npoints):
      rhol = s.dw[i]*s.wc[i] + s.do[i]*(1-s.wc[i])
      rho  = s.dg[i]*s.gvf[i]/100. + rhol*(1-s.gvf[i]/100.)  # gvf in is % ...
      mul  = s.vw[i]*s.wc[i] + s.vo[i]*(1-s.wc[i])
      mu   = s.vg[i]*s.gvf[i]/100. + mul*(1-s.gvf[i]/100.)  # gvf in is % ...
      dp[i] = VC.rcp_dp3(par, rho, cfv*mu, cfq*s.qt[i])
   return dp
   
if __name__ == '__main__':

   c     = UT.InputValues(sys.argv[1])
   psets = get_parametersets(sys.argv[2])
   xvar  = sys.argv[3]
   yvar  = sys.argv[4]

   f1 = figure().number
   if c.reidarsplot: f2 = figure().number

   for ps in psets:
      if not ps.active: continue
      s = read_data(c, ps)
      plotit(xvar, yvar, s, c, ps, f1)
      if c.reidarsplot:
         figure(f2)
         loglog(s.re, s.h_div_h0, ps.fmt, label=ps.label, ms=c.markersize)

   figure(f1)
   xlabel('flow rate')
   ylabel('dP [bar]')
   grid(True)
   legend(loc='best')
   title(UT.basename(sys.argv[2]))

   if c.reidarsplot:
      figure(f2)
      re = r_[1,199] 
      loglog(re, reidars_model(re), '--k', label='Reidars model')
      re = r_[201, 2999]
      loglog(re, reidars_model(re), '--k')
      re = r_[3001, 1e6]
      loglog(re, reidars_model(re), '--k')
      xlabel('Re [-]')
      ylabel('h/h0')
      grid(True)
      legend(loc='best')
      title(UT.basename(sys.argv[2]))
      ylim(0.1, 1)

   show()
