'''
Feb 2017
Atle J. Gyllensten
agy@statoil.com

Generic module to hold typical well data. See WellData
Assumes that all variables has the same sampling as the time vector. (This is *not* usually the case)

Latest version:
   2019-09-16
   Latest version:
      2019-09-16 
'''

import sys, os
import pylab as pl
import time, datetime
import AtlejgTools.SimulationTools.UnitConversion as U
import AtlejgTools.Utils                          as UT
from scipy.signal import medfilt
import pandas as PD

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
      gvf = gas volume fraction [in percent]
      gor = gas-oil ratio
      T   = temperature
   So, qo is oil rate, ql is liquid rate etc.
   In addition, qoc is cumulutive oil rate (etc.) if this is asked for
   '''
#
   def __init__(self, wellnm, startdate=0, t=[], qo=[], qg=[], qw=[], ql=[], p=[], T=[], gor=[], wc=[],
                visc_func=unitf, dens_func=unitf, unit_conv=unitf, welltype=TYPE_PROD,
                z_gauge=0., z_ref=0., rho_ref=1000., dt=None, cumul=True, dt_medfilt=None):
      '''
      # input
         wellnm        : well name
         startdate     : number of days since 0001-01-01 00:00:00 UTC (like date2num)
                         if 'today', use today's date
                         if 0, the time-vector should be made by date2num
         dt            : if regular sampling. [days]
         cumul         : calculate cumulative flowrates?
         dt_medfilt    : length of median filter window. no filtering if None
                         note that if filtering is done here, the raw data is not kept
                         if you want both series, use medfilt explicitely later
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
      # initialize
      self.qo = qo if len(qo) else pl.zeros(self.np)
      if len(qg):
         self.qg = qg
      else:
         if len(gor): self.qg = self.qo*gor
         else:   self.qg = pl.zeros(self.np)
      if len(qw):
         self.qw = qw
      else:
         if len(wc): self.qw = self.qo*wc/(1-wc)
         else:  self.qw = pl.zeros(self.np)
      self.p  = p  if len(p)  else 200*pl.ones(self.np)
      self.T  = T  if len(T)  else 60*pl.ones(self.np)
      if dt_medfilt != None:
         self.qo = self.medfilt('qo', dt_medfilt)
         self.qg = self.medfilt('qg', dt_medfilt)
         self.qw = self.medfilt('qw', dt_medfilt)
         self.p  = self.medfilt('p', dt_medfilt)
         self.T  = self.medfilt('T', dt_medfilt)
      self.typ = welltype
      #
      # adjust pressure to ref-level
      self.p -= rho_ref*GRAVITY*(z_gauge-z_ref) / U.BAR
      #
      # derived data
      self.ql = ql if len(ql) else self.qo + self.qw
      self.qt = self.ql + self.qg
      if pl.nanmax(self.ql) > 0:
         self.wc = self.qw/(self.ql)
      else:
         self.wc = pl.zeros(self.np)
      self.wcp = 100 * self.wc  # wc in percent
      self.gvf = self.qg/self.qt * 100   # in percent
      self.gor = self.qg/self.qo
      self.glr = self.qg/self.ql
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
      if cumul: self.calc_cumulatives()
#
   def calc_cumulatives(self):
      self.qoc = UT.cumulative(self.t, self.qo)
      self.qwc = UT.cumulative(self.t, self.qw)
      self.qgc = UT.cumulative(self.t, self.qg)
      self.qlc = UT.cumulative(self.t, self.ql)
      self.qtc = UT.cumulative(self.t, self.qt)

   def plot(self, varnm1, varnm2=None, ylim1=None, ylim2=None, sharex=None,
            sharey=None, fmt='-', ix1=0, ix2=-1, fname=None, showit=True,
            titl='', grid=True, plot_date=True):
      '''
      plots one or two variables in one figure.
      # input
         varnm1 : first variable to plot
         varnm2 : second variable to plot (optional)
         ylim1  : min, max for varnm1
         ylim2  : min, max for varnm2
         ix1    : index for start of plotting
         ix2    : index for end of plotting
         showit : show figure (boolean)
         fname  : save figure to this filename
      # output
         ax1, ax2: axes - one for the first variable and one for the second (if applicable)
                   useful if you need to link zooming in different figures.
      '''
      # plot first variable (always)
      y1 = self.get(varnm1)[ix1:ix2]
      ax1 = pl.figure().add_subplot(111, sharex=sharex, sharey=sharey)
      plotfunc = ax1.plot_date if plot_date else ax1.plot
      plotfunc(self.t[ix1:ix2], y1, 'k'+fmt, label=varnm1)
      ax1.set_ylabel(varnm1)
      if varnm2: ax1.legend(loc='upper left')
      if not ylim1 == None: ax1.set_ylim(ylim1[0],ylim1[1])
      pl.gcf().autofmt_xdate()                    # beautify
      ax1.set_title('Well %s. %s' % (self.wellnm, titl))
      ax1.grid(grid)
      #
      # plot second variable (optional)
      if varnm2:
         y2 = self.get(varnm2)[ix1:ix2]
         ax2 = ax1.twinx()
         plotfunc = ax2.plot_date if plot_date else ax2.plot
         plotfunc(self.t[ix1:ix2], y2, 'r'+fmt, label=varnm2)
         ax2.set_ylabel(varnm2)
         ax2.legend(loc='upper right')
         if not ylim2 == None: ax2.set_ylim(ylim2[0],ylim2[1])
      #
      if fname: pl.savefig(fname)
      if showit: pl.show()
      else:      pl.close()
      if varnm2: return ax1, ax2
      else     : return ax1
#
   def get(self, varnm, minval=-pl.Inf, maxval=pl.Inf):
      '''
      a clean way of getting data.
      will decimate data based on variable values if asked for
      '''
      y = self.__dict__[varnm]
      if minval == -pl.Inf and maxval == pl.Inf: return y
      ixs1 = pl.find(minval <= y)
      ixs2 = pl.find(maxval >= y)
      ixs = list(set.intersection(set(ixs1), set(ixs2)))
      return self.t[ixs], y[ixs], ixs
#
   def cum(self, varnm):
      '''
      return cumulative value
      '''
      return UT.cumulative(self.t, self.get(varnm))
#
   def medfilt(self, varnm, dt_filt, keepit=False):
      '''
      median filter
      # input
         varnm  : variable name ('p', 'qw' ....)
         dt_filt: length of median filter window [days]
         keepit : if True, a filtered variable will be kept with the extension 'f' - f.ex. qo -> qof
      '''
      winsz = int(dt_filt / self.dt)
      if winsz % 2 == 0: winsz += 1   # make it odd
      fy = medfilt(self.get(varnm), winsz)
      if keepit: self.__dict__[varnm+'f'] = fy
      return fy
#
   def calc_steps(self, flowvar, q_shut, q_min, dp_min, dt0, dt1, dt2, grad_lim, q_monoton=True):
      '''
      finds steps in flowrate, typically before / after shutins.
      units: any. just make sure it's consistent
      the steps are stored in two list: step_ups and step_downs
      # input
         flowvar   : typically 'ql' or 'qt'
         q_shut    : shut-level. ideally 0...                               [m3/d]
         q_min     : discard if production-level is below this              [m3/d]
         q_min     : discard if dp is less than this                        [bar]
         dt0       : smoothing window length when looking for hi/lo levels.
                     must be > sampling-length, but dont set it too low.    [d]
         dt1       : min. stable period before jump. must be > dt0.         [d]
         dt2       : min. stable period after jump. must be > dt0.          [d]
         grad_lim  : how steep should the jump be?                          [m3/d / d]
         q_monoton : boolean. make sure flow-rate is also monotonic
      # output
        none. but it sets self.step_downs and self.step_ups
      # note1
        the basic idea is to compare values of flowrate smoothed left (backwards) to
        flowrate smoothed right (forwards). if the difference is big, it means we have a
        major step in flowrate. then we filter these steps to make sure periods before and
        after the step are stable (enough).
      # note 2
        this could be quite time consuming. may use get_steps to filter 'unwanted' data.
      '''
      #
      # prepare data
      qs = UT.interp_nan(self.t, self.get(flowvar))            # interpolate NaN's and make copy
      qs[pl.find(qs < q_shut)] = 0.
      n0 = int(dt0 / self.dt)
      if n0 % 2: n0 += 1                               # must be even
      n = n0/2  # convinient
      n1 = int(dt1 / self.dt)
      n2 = int(dt2 / self.dt)
      #
      # analyze
      #  - first smooth it
      qfs = PD.rolling_mean(qs, n0, center=True)      # flow-rates filtered. smoothing forward (righty)
      #  - study the difference between 'lefty' and 'righty' smoothed curvs
      dqs = pl.zeros(self.np)
      for i in range(n, self.np-n): dqs[i] = qfs[i+n] - qfs[i+n-n0]
      #  - pad with constants
      dqs[:n0] = dqs[n0]
      dqs[-n0:] = dqs[-n0]
      dqs[pl.find(abs(dqs)<q_min)] = 0.                    # avoid all small variations
      #
      #  - find extremas of the difference. these are ramp-ups / ramp-downs
      rds = UT.find_minima(dqs, n)                         # ramp-down
      rus = UT.find_maxima(dqs, n)                         # ramp-ups
      #
      # want index where flow is below q_shut for the last / first time for ramp-ups / downs, respectively
      for r in [rds, rus]:
         i = -1 if r == rus else 0                         # last/first index
         tmp = []
         for ix in r:
            qq = qs[ix-n:ix+n]                             # study this window
            ii = pl.find(qq<=q_shut)
            if len(ii): tmp.append(ii[i]+ix-n)             # last/first time below limit
         if r == rds: rds = tmp
         else:        rus = tmp
      #
      # make sure the periodes before and after jump is long enough,
      # that the jump is steep enough, and monotonic behavior
      #  - ramp-downs
      self.step_downs = []
      for ix in rds:
         if pl.mean(qs[ix:ix+n2]) > q_shut:     continue   # not shut-in after jump
         grad = abs(dqs[ix]) / self.dt
         if grad < grad_lim:                    continue   # shutting too slowly
         qm  = pl.mean(qs[ix+n2:ix])                       # long open window before shutin
         if qm < q_min:                         continue   # rate too small
         pm1 = pl.mean(self.p[ix+n:ix+2*n])                # short periode in the early part
         pm2 = pl.mean(self.p[ix+n2-n:ix+n2])              # short periode at the end
         if pm1 - pm2 < dp_min:                 continue   # too small step (or non-monotonic)
         if q_monoton:
            qm1 = pl.mean(qfs[ix+n:ix+2*n])                # short periode in the early part
            qm2 = pl.mean(qfs[ix+n2-n:ix+n2])              # short periode at the end
            if qm1 < qm2:                       continue   # not monotonic in flowrate
         # if we get here, this index is ok
         self.step_downs.append(ix)
      #  - ramp-ups
      self.step_ups = []
      for ix in rus:
         if pl.mean(qs[ix-n1:ix]) > q_shut:     continue   # not shut-in before jump
         grad = abs(dqs[ix]) / self.dt
         if grad < grad_lim:                    continue   # opening too slowly
         qm  = pl.mean(qs[ix+1:ix+1+n2])                   # long open window after shutin
         if qm < q_min:                         continue   # rate too small
         pm1 = pl.mean(self.p[ix+n:ix+2*n])                # short periode in the early part
         pm2 = pl.mean(self.p[ix+n2-n:ix+n2])              # short periode at the end
         if pm2 - pm1 < dp_min:                 continue   # too small step (or non-monotonic)
         if q_monoton:
            qm1 = pl.mean(qfs[ix+n:ix+2*n])                # short periode in the early part
            qm2 = pl.mean(qfs[ix+n2-n:ix+n2])              # short periode at the end
            if qm1 > qm2:                       continue   # not monotonic in flowrate
         # if we get here, this index is ok
         self.step_ups.append(ix)
      #
      # also keep the filtered data
      self.__dict__[flowvar+'f'] = qfs
#
   def get_steps(self, dt1, dt2, p1, p2, direction='up'):
      '''
      *must* call calc_steps before this one
      # input
         dt1      : min. stable period before jump. must be > sampling-length
         dt2      : min. stable period after jump. must be > sampling-length
         p1       : average pressure must above this level before shutin (using dt1-period)
         p2       : average pressure must above this level after shutin (using dt2-period)
         direction: 'up' or 'down' for step_ups and step_downs, respectively
      # output
         list of step-indices
      # note:
         -if any of dt1, dt2, p1, p2 is None, then it will just give step_ups or step_downs (unfiltered)
         -calc_steps also takes dt1 and dt2 as input, so it does not make sense to have dt1,dt2 wider than this
         -this method is provided to do quick 'filtering' since calc_steps is soo slow
      '''
      #
      if direction == 'up': ixs = self.step_ups
      else:                 ixs = self.step_downs
      if dt1 == None or dt2==None or p1 == None or p2 == None:
         return ixs
      n1 = int(dt1/self.dt)
      n2 = int(dt2/self.dt)
      steps = []
      for ix in ixs:
         if pl.mean(self.p[ix-n1:ix]) >= p1 and pl.mean(self.p[ix:ix+n2]) >= p2:
            steps.append(ix)
      return steps
#
   def get_shutins(self, dt1=None, dt2=None, p1=None, p2=None):
      '''
      see docs for get_steps
      '''
      if self.typ == TYPE_PROD: self.shutins = self.get_steps(dt1, dt2, p1, p2, direction='up')
      else                    : self.shutins = self.get_steps(dt1, dt2, p1, p2, direction='down')
      return self.shutins
#
   def get_startups(self, dt1=None, dt2=None, p1=None, p2=None):
      '''
      see docs for get_steps
      '''
      if self.typ == TYPE_PROD: self.startups = self.get_steps(dt1, dt2, p1, p2, direction='down')
      else                    : self.startups = self.get_steps(dt1, dt2, p1, p2, direction='up')
      return self.startups
#
   def datestr(self, ix):
      '''
      gives a pleasent datestring for the given time index
      '''
      return pl.num2date(self.t[ix]).ctime()
#
   def inspect_steps(self, ixs, flowvar, dt1, dt2, figsdir, ask_to_keep=True):
      '''
      manually inspect given flowrate-steps. keep only the approved.
      # input
       ixs    : indices of steps. typically found using calc_steps
       flowvar : typically 'ql' or 'qt'
       dt1    : period before step
       dt2    : period after step
       figsdir: it will save pics here. must be unique!
      # output
       filtered index-list
       and pics on the format <figsdir>/<wellnm>_<flowvar>_<index>.png
      note1
       i was hoping to use python-figures, but in interactive mode these are
       'blank'. so i have to save figure to file and show these files.
      '''
      n1 = int(dt1/self.dt)
      n2 = int(dt2/self.dt)
      kept = []
      os.mkdir(figsdir)
      for ix in ixs:
         fname = '%s/%s_%s_%06i.png'%(figsdir, self.wellnm, flowvar, ix)
         self.plot(flowvar, 'p', ix1=ix-n1, ix2=ix+n2, showit=0, fname=fname, titl=self.datestr(ix))
         if 'win' in sys.platform: os.system(fname.replace('/','\\'))
         else                    : UT.tcsh('display %s &'%sfname)
         if ask_to_keep:
            ans = raw_input('Keep this [y]/n ? ')
            if ans and ans in ('n', 'N'):
               os.unlink(fname)
               continue
         print 'Keeping %s' % fname
         kept.append(ix)
      return kept
#
   def extend(self, t, value=0):
      '''
      extend variables with the given value so that it has values for all times t.
      assumes the given t has same sampling as self.t
      '''
      for varnm, y in self.__dict__.items():
         if not (type(y) == pl.ndarray and len(y) == self.np) or varnm == 't': continue
         # initate constant vector
         yy = value*pl.ones(t.shape)
         # find interval where will keep existing values
         ixs = pl.find(t<self.t[0])
         i1 = max(ixs) if len(ixs) else 0
         ixs = pl.find(t>self.t[-1])
         i2 = min(ixs) - 1 if len(ixs) else -1
         yy[i1:i2] = y
         self.__dict__[varnm] = yy
      self.t = t
#
   def resample(self, t):
      '''
      assumes it has regular sampling
      this is linear resampling without any filtering or smoothing. use with care!
      resamples all variables.
      '''
      for varnm, y in self.__dict__.items():
         if not (type(y) == pl.ndarray and len(y) == self.np) or varnm == 't': continue
         self.__dict__[varnm] = interp(t, self.t, y)
      self.t = t
      self.np = len(t)
      self.dt = (t[1] - t[0]) / 24. # days --> hours
#
   def step_vector(self, varnm, lvl, y1=0, y2=1):
      '''
      create a vector that has the value y1 when given variable is below given level, and y2 above
      '''
      y = self.get(varnm)
      step = y1*pl.ones(y.shape)
      ixs = pl.find(y > lvl)
      step[ixs] = y2
      return step
#
   def tix(self, t0, dt_max=pl.Inf):
      '''
      finds index closest to t0. (tix = time-index)
      # input
         t0    : could be datetime.datetime or float-representation
         dt_max: use this to make sure you are 'close enough'
      '''
      if type(t0) is datetime.datetime: t0 = pl.date2num(t0)
      ix = pl.searchsorted(self.t, t0)
      if ix == self.np: ix -= 1     # can not go too far
      if abs(t0-self.t[ix]) > dt_max:
         print 'tix: warning: did not find this t0'
         return None
      return ix

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
    def calc_shutins_old(self, q_shut, p_min, dt_min, dt_filt, q_preshut, dt_preshut):
      deprecated: see calc_shutins instead
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
      #
      # some useful parameters
      dt_min     /= 24.                        # convert to days
      dt_preshut /= 24.                        # convert to days
      nps = max(int(dt_preshut / self.dt), 1)  # make sure it's not 0
      #
      # filter data
      pf  = self.medfilt('p', dt_filt/24., keepit=True)
      qlf = self.medfilt('ql', dt_filt/24., keepit=True)
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
#
     return self.t[ixs], self.qw[ixs]/dps[ixs]
'''



   
if __name__ == '__main__':

   wd = WellData('A11', startdate='today', qo=800+200*rand(500))
   wd.plot('ql')
   pl.ylim(0,1100)
   pl.show()
