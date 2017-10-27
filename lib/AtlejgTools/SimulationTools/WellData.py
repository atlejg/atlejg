'''
Feb 2017
Atle J. Gyllensten
agy@statoil.com

Generic module to hold typical well data. See WellData


'''

import pylab as pl
import time, datetime

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
   def __init__(self, wellnm, startdate=0, t=[], qo=[], qg=[], qw=[], ql=[], p=[], T=[], visc_func=unitf, dens_func=unitf, unit_conv=unitf):
      '''
      # input
         wellnm        : well name
         startdate     : number of days since 0001-01-01 00:00:00 UTC (like date2num)
                         if 'today', use today's date
                         if 0, the time-vector should be made by date2num
      '''
      self.wellnm    = wellnm
      if startdate == 'today':
         now = time.localtime()
         startdate = pl.date2num(datetime.datetime(now.tm_year, now.tm_mon, now.tm_mday))
      #
      self.np = max(len(t), len(qo), len(qg), len(qw), len(ql), len(p), len(T))
      if self.np == 0: raise Exception('No data provdided')
      self.t  = startdate + (t  if len(t)  else pl.arange(self.np))
      self.qo = qo if len(qo) else pl.zeros(self.np)
      self.qg = qg if len(qg) else pl.zeros(self.np)
      self.qw = qw if len(qw) else pl.zeros(self.np)
      self.p  = p  if len(p)  else 200*pl.ones(self.np)
      self.T  = T  if len(T)  else 60*pl.ones(self.np)
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
#
   def plot(self, varnm1, varnm2=None, lim1=None, lim2=None, sharex=None, sharey=None, fmt='-'):
      '''
      plots one or two variables in one figure.
      # input
         varnm1 : first variable to plot
         varnm2 : second variable to plot (optional)
      # output
         ax1, ax2: axes - one for the first variable and one for the second (if applicable)
                   useful if you need to link zooming in different figures.
      '''
      # plot first variable (always)
      y1 = self.__dict__[varnm1]
      ax1 = pl.figure().add_subplot(111, sharex=sharex, sharey=sharey)
      ax1.plot_date(self.t, y1, 'k'+fmt, label=varnm1)
      ax1.set_ylabel(varnm1)
      if varnm2: ax1.legend(loc='upper left')
      if not lim1 == None: ax1.set_ylim(0,lim1)
      pl.gcf().autofmt_xdate()                    # beautify
      ax1.set_title(self.wellnm)
      ax1.grid(1)
      #
      # plot second variable (optional)
      if varnm2:
         y2 = self.__dict__[varnm2]
         ax2 = ax1.twinx()
         ax2.plot_date(self.t, y2, 'r'+fmt, label=varnm2)
         ax2.set_ylabel(varnm2)
         ax2.legend(loc='upper right')
         if not lim2 == None: ax2.set_ylim(0,lim2)
      #
      pl.show()
      if varnm2: return ax1, ax2
      else     : return ax1


   
if __name__ == '__main__':

   wd = WellData('A11', startdate='today', qo=800+200*rand(500))
   wd.plot('ql')
   pl.ylim(0,1100)
   pl.show()
