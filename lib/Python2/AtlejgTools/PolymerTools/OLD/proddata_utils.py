#import proddata
import pandas   as pd
#import shelve
import h5py
import AtlejgTools.SimulationTools.WellData as WellData
import AtlejgTools.Utils                    as UT
import datetime
from pylab import *
from scipy.interpolate import interp1d
from scipy.signal import medfilt


PPM_M       = 300e3 # PPM_M: ppm in mother-solution
MAX_VISC    = 35.

_wc_funcs = {}      # caching for performance

def get_bsw_wct(fname, sheet, start_ln=9, t_col='A', wc_col='C', winsz=31):
   '''
   assumes excel file has the format defined in the BS&W files
   retuns
    - wc-function
    - wc array (raw)
    - wc array (filtered)
   '''
   t1900 = date2num(datetime.datetime(1900,1,1))  # excel starts at 1900
   t   = UT.read_excel_column(sheet, t_col, start_ln, 99999, file=fname, return_array=1) + t1900
   ixs = find(diff(t)>0) + 1                      # make sure time is increasing. we do that by throwing out outliers
   t = t[ixs]
   wc  = UT.read_excel_column(sheet, wc_col, start_ln, 99999, file=fname)
   wc  = array([x if isreal(x) else NaN for x in wc])[ixs] / 100.  # dont want percent
   wcf = medfilt(wc, winsz)
   ixs = find(~isnan(wcf))
   wc_func = interp1d(t[ixs], wcf[ixs])
   return wc_func, t[ixs], wc[ixs], wcf[ixs]

def get_bsw_wct2(wellnm, winsz=31, wc_func_only=True, path=r"\\statoil.net\unix_st\Project\peregrino\ressim\2016b\sw_2016b_ver002\eclref\polypilot"):
   '''
   gets bsw water-cut from yngve's files.
   assumes ascii-file has name <wellnm>_BSW.txt (like A-22_BSW.txt) with format
   mm/dd/yyyy WCT like this:
      DATE            WWCTH A-18
      TIME            FRACTION
      2/1/2014	0.005
      2/2/2014	0.006
      2/3/2014	0.01
      2/4/2014	0.012
      2/5/2014	0.01
      .
      .
      .
   retuns
    - wc-function
    and optionally (if wc_func_only is False):
    - date
    - wc array (raw)
    - wc array (filtered)
   '''
   if wc_func_only and _wc_funcs.has_key(wellnm): return _wc_funcs[wellnm] # for performance...
   # date-utility function
   def _date2num(date):                                             # date on format mm/dd/yyyy
      m, d, y = [int(x) for x in date.split('/')]
      return date2num(datetime.datetime(y,m,d))
   #
   if '-' not in wellnm: wellnm = '%s-%s' % (wellnm[0], wellnm[1:]) # convert A22 to A-22 if needed
   fnm = r'%s\%s_BSW.txt' % (path, wellnm)
   d   = loadtxt(fnm, skiprows=3, converters={0:_date2num})         # skip 3 first rows (header is 2 or 3 lines long..)
   t   = d[:,0]
   ixs = argsort(t)                                                 # make sure time is increasing. we do that by sorting
   t = t[ixs]
   wc  = d[:,1][ixs]
   wcf = medfilt(wc, winsz)
   wc_func = interp1d(t, wcf)
   _wc_funcs[wellnm] = wc_func
   if wc_func_only: return wc_func
   else           : return wc_func, t, wc, wcf

def adjust_injrate1(t, q):
   # adjust injection rate in-place
   t1, t2 = [date2num(datetime.datetime(2018,1,6)), date2num(datetime.datetime(2018,4,15))]
   ixs1 = find(t>=t1)
   ixs2 = find(t<=t2)
   ixs  = list(set.intersection(set(ixs1), set(ixs2)))
   q[ixs] *= 1.6

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
      qm  = zeros(len(tqm))
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
   t  = linspace(t1,t2, n)
   print 'last date', num2date(t[-1])
   qc = interp(t, tqc, qc)
   qi = interp(t, tqi, qi)
   qm = interp(t, tqm, qm)
   qd = interp(t, tqd, qd)
   qs = interp(t, tqs, qs)
   p  = interp(t, tp, p)
   #
   #qs = qi + qm + qd   # skid flow rate
   qt = qc + qs        # total flow rate. qc is 0 when qs > 0 and vice versa
   adjust_injrate1(t, qt)
   #
   #    viscosity
   ppm = PPM_M * qm / qs
   ppm[find(ppm>PPM_M/60.)] = 0.
   ppm[find(ppm<0.)]        = 0.
   ppm[find(qs<200.)]       = 0.
   visc = 4e-6*ppm**2 - 0.0029*ppm + 2.4097   # excel trendline from kjetil E spread-sheet 
   visc[find(visc<0)]       = NaN
   visc[find(visc>MAX_VISC)]      = NaN
   #    create WellData object
   a11 = WellData.WellData('A11', welltype='inj', t=t, qw=qt, p=p, dt=dt)
   #    add some properties
   a11.ppm  = ppm
   a11.visc = visc 
   a11.qi   = qi
   a11.qd   = qd
   a11.qm   = qm
   '''
   skip A18 for now
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
   dt = max(mean(diff(tp)), mean(diff(to)), mean(diff(tw)))
   t = arange(t1,t2, dt)
   p = interp(t, tp, p)
   qo = interp(t, to, qo)
   qw = interp(t, tw, qw)
   a18 = WellData.WellData('A18', welltype='prod', t=t, qo=qo, qw=qw, p=p, dt=dt)
   '''
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
   dt = max(mean(diff(tp)), mean(diff(to)), mean(diff(tw)))
   t = arange(t1,t2, dt)
   p = interp(t, tp, p)
   qo = interp(t, to, qo)
   qw = interp(t, tw, qw)
   a22 = WellData.WellData('A22', welltype='prod', t=t, qo=qo, qw=qw, p=p, dt=dt)
   #
   # logistics
   #
   h5.close()
   #
   '''
   shelving fails: out of memory
   if sleeptime >= -1:
      # save to file
      db = shelve.open(SHELVE_FILE)
      db['A11'] = a11
      db['A18'] = a18
      db['A22'] = a22
      db.close()
      print 'Updating shelve database %s' % SHELVE_FILE
   '''
   #return {'A-11':a11, 'A-22':a22, 'A-18':a18}   # poor man's shelving
   return {'A-11':a11, 'A-22':a22}   # poor man's shelving


def t2str(t):
   '''
   handle yngve's date-format which is like 01.09.2014 22:42
   this takes the float-representation into yngve-format
   '''
   t = num2date(t)
   return '%02i-%02i-%02i %02i:%02i:%02i'%(t.day, t.month, t.year, t.hour, t.minute, t.second)

def str2t(s):
   '''
   handle yngve's date-format which is like 01.09.2014 22:42
   this takes the yngve-format into float-representation 
   '''
   dd, t = s.split()
   d, m, y = dd.split('.')
   h, mm  = t.split(':')
   d, m, y, h, mm = int(d), int(m), int(y), int(h), int(mm)
   return date2num(datetime.datetime(y, m, d, h, mm))


def get_indices_from_filenames(pattern, sep='_', pos=-1):
   return [int(UT.basename(x).split(sep)[pos]) for x in UT.glob(pattern, sortit=1)]


def create_shutinpressure_csv(w, pattern, t0):
   '''
   creates csv-file with shutinpressures. assumes that plotshutins has been run
   and then a number of these plots have been removed manually by visual inspection.
   w: WellData instance
   pattern: typically 'Figs/a22_*.png'
   '''
   if not w.shutins: w.calc_shutins(20., 60, 24., 6., 500., 48.)
   #ixs = [int(x.split('_')[1]) for x in UT.glob(pattern, sortit=1)]
   ixs = get_indices_from_filenames(pattern)
   fname = '%s_shutins.csv' % w.wellnm
   f = open(fname, 'w')
   f.write(
'''# COLUMN INFO:
# time;presssure
#
''')
   ts, ps = w.get_shutins(t0)[:2]
   for ix in ixs:
       f.write('%s; %.1f\n' % (t2str(ts[ix]), ps[ix]))
   f.close()
   print fname, 'was created'

def plot_csv_file(csvfile):
   figure()
   to_datetime = lambda x: pd.to_datetime(x, dayfirst=True)
   dd = pd.read_csv(csvfile, sep=';', skiprows=3, header=None, converters={0:to_datetime})
   plot_date(dd[0], dd[1], 'o--', ms=8, lw=2)
   title(csvfile)
   grid(True)
   show()

# main logic 

if __name__ == '__main__':

   if len(sys.argv) > 1:

      t_start = date2num(datetime.datetime(2017, 10, 5, 0, 0))   # use vegard's values before this date?
      t_start = date2num(datetime.datetime(1899, 1, 1, 0, 0))
      wellnm      = sys.argv[1]
      mode        = sys.argv[2]

      if not 'wells' in dir(): wells = read_pilot_area_wells()

      w = wells[wellnm]


      if mode == 'writeshutins':
         if not w.shutins: w.calc_shutins(20., 60, 24., 6., 500., 48.)
         ts, ps0 = w.get_shutins(0.5)[:2]
         ps1     = w.get_shutins(1.0)[1]
         ps2     = w.get_shutins(2.0)[1]
         ps4     = w.get_shutins(4.0)[1]
         fname = '%s_shutins.csv' % wellnm
         f = open(fname, 'w')
         for i in range(len(ts)):
            if ts[i] < t_start: continue
            f.write('%s; %.1f; %.1f; %.1f; %.1f\n' % (t2str(ts[i]), ps0[i], ps1[i], ps2[i], ps4[i]))
         f.close()
         print fname, 'was created'

      if mode == 'plotshutins':
         if not w.shutins: w.calc_shutins(20., 60, 24., 6., 500., 48.)
         dt_before = float(sys.argv[3])  # days
         dt_after  = float(sys.argv[4])  # days
         na = int(dt_after / w.dt)
         nb = int(dt_before / w.dt)
         ts, ps, ixs = w.get_shutins(0)
         for n, ix in enumerate(ixs):
            if ts[n] < t_start: continue
            ax1 = w.plot('qlf', 'pf', ix1=ix-nb, ix2=ix+na)[0]
            ax1.plot_date((w.t[ix],w.t[ix]), (0,max(w.qlf[ix-nb:ix+na])), 'b--') # vertical line
            ylabel('Pressure [bar]')
            fname = 'Figs/%s_%03i_shutin_%s.png'%(wellnm, n, num2date(w.t[ix]).isoformat()[:10])
            savefig(fname)
            print fname, 'created'
            close()
