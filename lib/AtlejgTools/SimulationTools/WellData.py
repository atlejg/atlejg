'''
Feb 2017
Atle J. Gyllensten
agy@statoil.com

Generic module to hold typical well data. See WellData
Assumes that all variables has the same sampling as the time vector. (This is *not* usually the case)

'''

import pylab as pl
import time, datetime
import AtlejgTools.SimulationTools.UnitConversion as U
from scipy.signal import medfilt

TYPE_PROD = 'prod'
TYPE_INJ  = 'inj'
GRAVITY   = 9.8

def unitf(*args):
   '''
   unit function. accepts any number of arguments.
   '''
   return 1.

class WellData(object):
   '''
   Generic class to hold typical well data
   Given individual rates for oil, gas, water, it calculates useful vectors like wct, gvf etc.
   Nomenclature:
      t   = time (number of days since 0001-01-01 00:00:00 UTC (like date2num))
      o   = oil
      g   = gas
      w   = water
      d?  = density
      v?  = viscosity
      q?  = volume flow
      m?  = mass flow
      p   = pressure
      ?l  = liquid (oil + water)
      ?t  = total  (oil + water + gas)
      np  = number of points (length of vectors)
      wc  = water cut
      gvf = gas volume fraction
      T   = temperature
   So, qo is oil rate, ql is liquid rate etc.
   '''
#
   def __init__(self, wellnm, startdate=0, t=[], qo=[], qg=[], qw=[], ql=[], p=[], T=[],
                visc_func=unitf, dens_func=unitf, unit_conv=unitf, welltype=TYPE_PROD,
                z_gauge=0., z_ref=0., rho_ref=1000., dt=None):
      '''
      # input
         wellnm        : well name
         startdate     : number of days since 0001-01-01 00:00:00 UTC (like date2num)
                         if 'today', use today's date
                         if 0, the time-vector should be made by date2num
      '''
      self.wellnm  = wellnm
      self.z_gauge = z_gauge
      self.rho_ref = rho_ref
      self.z_ref   = z_ref
      self.dt      = dt  # in many cases, we assume regular sampling
      if startdate == 'today': startdate = datetime.datetime.now().toordinal()
      #
      self.np = max(len(t), len(qo), len(qg), len(qw), len(ql), len(p), len(T))
      if self.np == 0: raise Exception('No data provdided')
      self.t  = startdate + (t  if len(t)  else pl.arange(self.np))
      self.qo = qo if len(qo) else pl.zeros(self.np)
      self.qg = qg if len(qg) else pl.zeros(self.np)
      self.qw = qw if len(qw) else pl.zeros(self.np)
      self.p  = p  if len(p)  else 200*pl.ones(self.np)
      self.T  = T  if len(T)  else 60*pl.ones(self.np)
      self.typ = welltype
      #
      # adjust pressure to ref-level
      self.p -= rho_ref*GRAVITY*(z_gauge-z_ref) / U.BAR
      #
      # derived data
      self.ql = ql if len(ql) else self.qo + self.qw
      self.qt = self.ql + self.qg
      if max(self.ql) > 0:
         self.wc = self.qw/(self.ql)
      else:
         self.wc = pl.zeros(self.np)
      self.wcp = 100 * self.wc  # wc in percent
      self.gvf = self.qg/self.qt * 100   # in percent
      self.vo = visc_func('oil', self.p, self.T)
      self.vg = visc_func('gas', self.p, self.T)
      self.vw = visc_func('wat', self.p, self.T)
      self.do = dens_func('oil', self.p, self.T)
      self.dg = dens_func('gas', self.p, self.T)
      self.dw = dens_func('wat', self.p, self.T)
      self.mo = self.qo * self.do
      self.mw = self.qw * self.dw
      self.mg = self.qg * self.dg
      self.ml = self.mo + self.mw
      self.mt = self.ml + self.mg
      self.shutins = None
#
   def plot(self, varnm1, varnm2=None, ylim1=None, ylim2=None, sharex=None, sharey=None, fmt='-', ix1=0, ix2=-1):
      '''
      plots one or two variables in one figure.
      # input
         varnm1 : first variable to plot
         varnm2 : second variable to plot (optional)
         ylim1  : min, max for varnm1
         ylim2  : min, max for varnm2
         ix1     : index for start of plotting
         ix2     : index for end of plotting
      # output
         ax1, ax2: axes - one for the first variable and one for the second (if applicable)
                   useful if you need to link zooming in different figures.
      '''
      # plot first variable (always)
      y1 = self.__dict__[varnm1][ix1:ix2]
      ax1 = pl.figure().add_subplot(111, sharex=sharex, sharey=sharey)
      ax1.plot_date(self.t[ix1:ix2], y1, 'k'+fmt, label=varnm1)
      ax1.set_ylabel(varnm1)
      if varnm2: ax1.legend(loc='upper left')
      if not ylim1 == None: ax1.set_ylim(ylim1[0],ylim1[1])
      pl.gcf().autofmt_xdate()                    # beautify
      ax1.set_title(self.wellnm)
      ax1.grid(1)
      #
      # plot second variable (optional)
      if varnm2:
         y2 = self.__dict__[varnm2][ix1:ix2]
         ax2 = ax1.twinx()
         ax2.plot_date(self.t[ix1:ix2], y2, 'r'+fmt, label=varnm2)
         ax2.set_ylabel(varnm2)
         ax2.legend(loc='upper right')
         if not ylim2 == None: ax2.set_ylim(ylim2[0],ylim2[1])
      #
      pl.show()
      if varnm2: return ax1, ax2
      else     : return ax1
#
   def get(self, varnm, minval=-pl.Inf, maxval=pl.Inf):
      '''
      decimate data based on variable values
      '''
      y = self.__dict__[varnm]
      ixs1 = pl.find(minval <= y)
      ixs2 = pl.find(maxval >= y)
      ixs = list(set.intersection(set(ixs1), set(ixs2)))
      return self.t[ixs], y[ixs], ixs
#
   def cum(self, varnm):
      '''
      return cumulative value
      '''
      return UT.cumulative(self.t, self.__dict__[varnm])
#
   def medfilt(self, varnm, dt_filt):
      '''
      median filter
      '''
      winsz = int(dt_filt / self.dt)
      if winsz % 2 == 0: winsz += 1   # make it odd
      y = self.__dict__[varnm]
      return medfilt(y, winsz)
#
   def calc_shutins(self, q_shut, p_min, dt_min, dt_filt, q_preshut, dt_preshut):
      '''
      finds the observed pressure during shut-ins
      # input
         q_shut    : max (total) production (to define shut-ins) [m3/d]
         p_min     : min BHP (for dismissing periods with no BHP-data)
         dt_min    : min time periode to be regarded a shut-in [hours]
         dt_filt   : time periode to filter the pressure [hours]
         q_preshut : min flow-rate in the periode before shut-in (should be high)
         dt_preshut: min length of hi-flowrate periode before shut-in [hours]
      # output
         t_sh : time for shut-ins
         p_sh : shut-in pressures
         ixs  : indices for measured shut-ins
      # note1
         assumes that shut-in starts first time flow-rate is below q_shut
        note2
         shutin-indices, filtered pressure and filtered rates are kept in
         self.shutins, self.pf, self.qlf
        note3
         we will miss last shut-in if it is still ongoing...
      '''
      #
      # some useful parameters
      dt_min     /= 24.                        # convert to days
      dt_preshut /= 24.                        # convert to days
      nps = max(int(dt_preshut / self.dt), 1)  # make sure it's not 0
      #
      # filter data
      pf  = self.medfilt('p', dt_filt/24.)
      qlf = self.medfilt('ql', dt_filt/24.)
      #
      shutins = [] # indice-pairs for shutins
      ix1 = None
      for i, q in enumerate(qlf):
         if q < q_shut:
            if ix1 == None: ix1 = i
         elif ix1 != None:
            if   (self.t[i-1]-self.t[ix1]) > dt_min                \
             and pl.mean(qlf[ix1-nps:ix1]) > q_preshut             \
             and ( (self.typ==TYPE_PROD and pf[i-1] > pf[ix1]) or  \
                   (self.typ==TYPE_INJ  and pf[i-1] < pf[ix1]) )   \
             and pf[i-1] > p_min :
               shutins.append((ix1,i-1))
            ix1 = None
      #
      # keep useful stuff
      self.shutins = shutins
      self.pf      = pf
      self.qlf     = qlf
#
   def get_shutins(self, tp):
      '''
      *must* call calc_shutins before this one
      tp : get pressure at this relative time [hours]
      '''
      tp        /= 24.                         # convert to days
      np  = int(tp / self.dt)
      ixs = [x[0]+np for x in self.shutins]
      return self.t[ixs], self.pf[ixs], ixs

   def resample(self, t):
      '''
      assumes dt has regular sampling
      '''
      for varnm, y in self.__dict__.items():
         if not (type(y) == pl.ndarray and len(y) == self.np) or varnm == 't': continue
         self.__dict__[varnm] = interp(t, self.t, y)
      self.t = t
      self.np = len(t)
      self.dt = (t[1] - t[0]) / 24. # days --> hours
'''
deprecated
   def estimate_PI(self, dp_min, ql_min, rel_smooth_len=0.1):
      tries to estimate PI (flowrate/dP)
      you need to estimate_p_reservoir first
      # input
         dp_min: avoid periods with very small dP
         ql_min: make sure we only consider periods when well is producing
      # output
         t  : time vector
         PI : estimated PI
      if self.p_res == None:
         print 'Need to run estimate_p_reservoir first!'
         return 0
      dp  = abs(self.p - self.p_res)
      sz = int(len(dp)*rel_smooth_len)               # window-size
      dps = pandas.rolling_mean(dp, sz, center=True)  # dp smoothed
      ixs1 = find(dps>dp_min)
      ixs2  = find(self.ql >= ql_min)
      ixs  = list(set.intersection(set(ixs1), set(ixs2)))
      ixs.sort()
      return self.t[ixs], self.qw[ixs]/dps[ixs]
'''



   
if __name__ == '__main__':

   wd = WellData('A11', startdate='today', qo=800+200*rand(500))
   wd.plot('ql')
   pl.ylim(0,1100)
   pl.show()
