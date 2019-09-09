from scipy.signal import medfilt
import pylab as pl
import h5py
from datetime import datetime
from scipy.interpolate import interp1d
import AtlejgTools.SimulationTools.WellData as WellData

PPM_M        = 300e3 # PPM_M: ppm in mother-solution
MAX_VISC     = 50.
MAX_INJ_RATE = 1250.
MONTH_MAP    = {'JAN':1, 'FEB':2, 'MAR':3, 'APR':4, 'MAY':5, 'JUN':6, 'JUL':7, 'AUG':8, 'SEP':9, 'OCT':10, 'NOV':11, 'DEC':12} # useful when converting dates

def date2num1(date):                                             # date on format mm/dd/yyyy
   '''
   DATE            WWCTH A-18
   TIME            FRACTION
   2/1/2014 0.005
   2/2/2014 0.006
   '''
   m, d, y = [int(x) for x in date.split('/')]
   return pl.date2num(datetime(y, m, d))

def date2num2(date):                                             # date on format mm/dd/yyyy
   '''
   DATE            WWCTH A-18
   TIME            FRACTION
   15-Aug-2014 0.005
   16-Aug-2014 0.006
   '''
   d, m, y = date.split('-')
   m = MONTH_MAP[m.upper()]
   return pl.date2num(datetime(2000+int(y), m, int(d)))

def get_bsw_wct(fnm, winsz=31, wc_func_only=True, date2num_func=date2num1):
   '''
   gets bsw water-cut from yngve's files.
   fnm: /private/agy/MyPolymerStuff/InputData/WCT_mesurements/A-22_BSW.txt
   retuns
    - wc-function
    and optionally (if wc_func_only is False):
    - date
    - wc array (raw)
    - wc array (filtered)
   '''
   #
   d   = pl.loadtxt(fnm, skiprows=3, converters={0:date2num_func})         # skip 3 first rows (header is 2 or 3 lines long..)
   t   = d[:,0]
   ixs = pl.argsort(t)                                                 # make sure time is increasing. we do that by sorting
   t = t[ixs]
   wc  = d[:,1][ixs]
   wcf = medfilt(wc, winsz)
   wc_func = interp1d(t, wcf)
   if wc_func_only: return wc_func
   else           : return wc_func, t, wc, wcf


def read_pilot_area_wells(db_file, include_mothersolution=True):
   '''
   skid: mother solution is often not available! but it is small, so we can ignore it
   '''
   h5  = h5py.File(db_file)
   # A11 stuff
   #    pressure
   tp = h5['A-11']['WBHP'][:,0]         # time - pressure series
   p  = h5['A-11']['WBHP'][:,1]
   #
   #    clamp-on rate
   tqc = h5['A-11']['WWIRHR'][:,0]       # time - inj-rate series
   qc  = h5['A-11']['WWIRHR'][:,1] * 24  # m3/h -> m3/d
   #
   #    skid: inversion water flow
   tqi = h5['A-11']['WWIRI'][:,0]       # time - inj-rate series
   qi  = h5['A-11']['WWIRI'][:,1] * 24  # m3/h -> m3/d
   #
   #    skid: dilution water flow
   tqd = h5['A-11']['WWIRD'][:,0]       # time - inj-rate series
   qd  = h5['A-11']['WWIRD'][:,1] * 24  # m3/h -> m3/d
   #
   #    skid: pump (mother solution). often not available!
   if include_mothersolution:
      tqm = h5['A-11']['WWIRP'][:,0]       # time - inj-rate series
      qm  = h5['A-11']['WWIRP'][:,1] * 24 / 1000. # l/h -> m3/d
   else:
      tqm = tqd
      qm  = pl.zeros(len(tqm))
   #
   #    skid: total liquid rate
   tqs = h5['A-11']['WWIRS'][:,0]       # time - inj-rate series
   qs  = h5['A-11']['WWIRS'][:,1] * 24  # m3/h -> m3/d
   #
   #    resampling and summing. note: we WILL extrapolate if some tag is missing data in either direction
   t1 = min(tqc[0], tqi[0], tqm[0], tqd[0], tqs[0])
   t2 = max(tqc[-1], tqi[-1], tqm[-1], tqd[-1], tqs[-1])
   n  = max(len(tqc), len(tqi), len(tqm), len(tqd), len(tqs))
   dt = (t2-t1) / (n-1)
   t  = pl.linspace(t1,t2, n)
   print 'last date', pl.num2date(t[-1])
   qc = pl.interp(t, tqc, qc)
   qi = pl.interp(t, tqi, qi)
   qm = pl.interp(t, tqm, qm)
   qd = pl.interp(t, tqd, qd)
   qs = pl.interp(t, tqs, qs)
   p  = pl.interp(t, tp, p)
   #
   #    viscosity
   ppm                          = PPM_M * qm / qs
   ppm[pl.find(ppm>PPM_M/60.)]  = 0.
   ppm[pl.find(ppm<0.)]         = 0.
   ppm[pl.find(qs<200.)]        = 0.
   visc                         = 4e-6*ppm**2 - 0.0029*ppm + 2.4097   # excel trendline from kjetil E spread-sheet 
   visc[pl.find(visc<0)]        = pl.NaN
   visc[pl.find(visc>MAX_VISC)] = pl.NaN
   # do some modifications since data is bad..
   #
   # must handle polymer-period 2 separately :-(
   ix = (t < pl.date2num(datetime(2019,2,13))).nonzero()[0][-1]
   qw = pl.concatenate((qc[:ix], qs[ix:]))
   # also: remove highest values (not realistic)
   qw = pl.minimum(qw, MAX_INJ_RATE)
   # for some reason it has inj-rate before well is opened...
   ixs = t > pl.date2num(datetime(2014,6,1))
   # adjust inj-rate for first polymer-phase
   t1, t2 = [pl.date2num(datetime(2018,1,6)), pl.date2num(datetime(2018,4,15))]
   ixs1 = pl.find(t>=t1)
   ixs2 = pl.find(t<=t2)
   ixs  = list(set.intersection(set(ixs1), set(ixs2)))
   qw[ixs] *= 1.6
   #
   #    create WellData object
   a11 = WellData.WellData('A11', welltype='inj', t=t, qw=qw, p=p, dt=dt)
   #    add some properties
   a11.ppm  = ppm
   a11.visc = visc 
   a11.qi   = qi
   a11.qd   = qd
   a11.qm   = qm
   a11.qs   = qs
   a11.qc   = qc
   # A18
   tp = h5['A-18']['WBHP'][:,0]
   p  = h5['A-18']['WBHP'][:,1]
   to = h5['A-18']['WOPR'][:,0]
   qo = h5['A-18']['WOPR'][:,1]
   tw = h5['A-18']['WWPR'][:,0]
   qw = h5['A-18']['WWPR'][:,1]
   # need to resample data to a common timeline covering all tags
   t1 = max(tp[0], to[0], tw[0])
   t2 = min(tp[-1], to[-1], tw[-1])
   dt = max(pl.mean(pl.diff(tp)), pl.mean(pl.diff(to)), pl.mean(pl.diff(tw)))
   t = pl.arange(t1,t2, dt)
   p = pl.interp(t, tp, p)
   qo = pl.interp(t, to, qo)
   qw = pl.interp(t, tw, qw)
   a18 = WellData.WellData('A18', welltype='prod', t=t, qo=qo, qw=qw, p=p, dt=dt)
   # A22
   tp = h5['A-22']['WBHP'][:,0]
   p  = h5['A-22']['WBHP'][:,1]
   to = h5['A-22']['WOPR'][:,0]
   qo = h5['A-22']['WOPR'][:,1]
   tw = h5['A-22']['WWPR'][:,0]
   qw = h5['A-22']['WWPR'][:,1]
   # need to resample data to a common timeline covering all tags
   t1 = max(tp[0], to[0], tw[0])
   t2 = min(tp[-1], to[-1], tw[-1])
   dt = max(pl.mean(pl.diff(tp)), pl.mean(pl.diff(to)), pl.mean(pl.diff(tw)))
   t = pl.arange(t1,t2, dt)
   p = pl.interp(t, tp, p)
   qo = pl.interp(t, to, qo)
   qw = pl.interp(t, tw, qw)
   a22 = WellData.WellData('A22', welltype='prod', t=t, qo=qo, qw=qw, p=p, dt=dt)
   #
   # logistics
   #
   h5.close()
   #
   return {'A-11':a11, 'A-22':a22, 'A-18':a18}   # poor man's shelving

def read_polymerconc_yngve(fnm='/project/peregrino/users/agy/InputData/PolymerConcentrations/polymer.dat'):
   '''
   yngve has a file with polymer concentration. this routine reads it
   '''
   lines = open(fnm).readlines()
   dates, concs = [],[]
   for line in lines[3:]:
      if line.startswith('--'): continue
      r = line.strip().split()
      if len(r) == 1: concs.append(float(r[0]))
      else:           dates.append(datetime(int(r[2]), MONTH_MAP[r[1].replace("'","")], int(r[0])))
   return dates, concs
