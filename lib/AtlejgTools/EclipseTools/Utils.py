'''
collection of useful routines when working with eclipse - typically file reader(s) etc.
atle j. gyllensten
agy@statoil.com
'''

import sys, tempfile, os, pdb, cPickle, re, time
import pylab as pl
import numpy
import AtlejgTools.LogTools
import AtlejgTools.SimulationTools.UnitConversion as U
import AtlejgTools.Utils                          as UT
import AtlejgTools.Plotting.Utils                 as PU
import ert
import ert.ecl.ecl                                as ecl
import inspect
import subprocess, signal, psutil
import ert.ecl.ecl as ecl

_log = AtlejgTools.LogTools.Log(AtlejgTools.LogTools.LOG_LEVEL_DEBUG)

SEPARATOR    = ':' # separator for variable names (like WBHP:OP1, where SEPARATOR is ':'. (Used to be '-')
JOINT_LENGTH = 12. # meter
DELIM        = '/' # end-of-record for Eclipse input
SECTIONS = ('RUNSPEC', 'GRID', 'EDIT', 'PROPS', 'REGIONS', 'SOLUTION', 'SUMMARY', 'SCHEDULE', 'OPTIMIZE')

def _read_lines(f, line):
   nvals_pr_column = 4 # hardcoded, will probably not change
   nvals = int(line[20:23])
   nlines = int((nvals-1) / nvals_pr_column) + 1 # clever...
   vals = []
   for lineno in range(nlines):
      line = f.readline()
      vals.extend(map(float, line.split()))
   return vals

def read_frft_file(fname, varnm, t_start=-pl.Inf, t_end=pl.Inf, dotranspose=True):
   '''
   Reads  an frft file, returning time, well-coordinates (md), and values for the requested variable pr time and md.
   See keyword WRFTPLT.
   '''
   f = open(fname)
   z = []; t = [] # init
   while True:
      line = f.readline()
      if not line: break
      if 'TIME' in line:
         t_ = _read_lines(f, line)[0]
         if   t_ > t_end  : break
         elif t_ < t_start: continue
      if varnm in line:
         # this way we make sure len(t) == len(z)
         z.append(_read_lines(f, line))
         t.append(t_)
      elif len(t) == 1 and 'CONLENEN' in line: # do it only once
         x = _read_lines(f, line)
   f.close()
   t = pl.array(t); x = pl.array(x); z = pl.array(z)
   if dotranspose: z = z.T
   return t, x, z

def read_rft_file_old(varnm, ef, rftfile=None, nwells=1, segments=None):
   '''
   get the requested variable. result will be put into a matrix
   with one row per timestep, and one column per segment.
   varnm is typically SEGORAT, CONPRES etc.
   a list of matrices is returned, one matrix for each well
   #
   really wanted to use ecl.EclRFTFile(rftfile), but it does not support SEGMENTS (yet).
   so this is sort of a hack. EclFile expects a restart file, but accepts RFT as well.
   problem is that multiple wells is not properly handled.
   so, here we split the data to get one matrix pr. well. this is potentially dangerous,
   since the way data is organized in the file could be subject to changes, i guess.
   #
   ef is obtained using ecl.EclFile. if not given, rftfile must be given.
   #
   if not segments is given, it will try to find which segments to use.
   #
   use keyword WRFTPLT to get the requested file from eclipse.
   '''
   if not ef: ef = ecl.EclFile(rftfile)
   n_tsteps   = ef.num_named_kw(varnm)/nwells  # ef.num_report_steps() fails for RFT-file
   if segments is None:
      n_well_segments = len(ef.iget_named_kw('CONKH', 0))
      n_segments = len(ef.iget_named_kw(varnm, 0))   # is n_well_segments or 2*n_well_segments+1
      if n_segments > n_well_segments:
         # usually we only want values from "icd-segments"
         j1 = n_well_segments + 1
      else:
         j1 = 1
      if varnm in ['TIME', 'DATE']:
         # some variables not defined for each segment
         j2 = 2
      else:
         j2 = j1 + n_well_segments
      segments = pl.arange(j1, j2)
   m = pl.zeros((n_tsteps*nwells, len(segments)))
   for tstep in range(n_tsteps*nwells):
      kw = ef.iget_named_kw(varnm, tstep)
      for j in range(len(segments)):
         m[tstep, j] = kw[segments[j]-1]
   # make a list - one matrix for each well
   l = []
   for i in range(nwells):
      indx = range(i,n_tsteps*nwells, nwells)
      l.append(m[indx,:])
   return l

def read_rft_file(varnm, ef, rftfile=None, nwells=1):
   '''
   get the requested variable. result will be put into a matrix
   with one row per timestep, and one column per segment.
   varnm is typically SEGORAT, CONPRES etc.
   a list of matrices is returned, one matrix for each well
   #
   really wanted to use ecl.EclRFTFile(rftfile), but it does not support SEGMENTS (yet).
   so this is sort of a hack. EclFile expects a restart file, but accepts RFT as well.
   problem is that multiple wells is not properly handled.
   so, here we split the data to get one matrix pr. well. this is potentially dangerous,
   since the way data is organized in the file could be subject to changes, i guess.
   #
   ef is obtained using ecl.EclFile. if not given, rftfile must be given.
   #
   use keyword WRFTPLT to get the requested file from eclipse.
   #
   it will find which segments to use.
   '''
   if not ef: ef = ecl.EclFile(rftfile)
   n_tsteps   = ef.num_named_kw(varnm)/nwells  # ef.num_report_steps() fails for RFT-file
   # make a list - one matrix for each well
   l = []
   for wellid in range(nwells):
      n_well_segments = len(ef.iget_named_kw('CONKH', wellid))
      n_segments = len(ef.iget_named_kw(varnm, wellid))   # is n_well_segments or 2*n_well_segments+1
      if n_segments > n_well_segments:
         # usually we only want values from "icd-segments"
         j1 = n_well_segments + 1
      else:
         j1 = 1
      if varnm in ['TIME', 'DATE']:
         # some variables not defined for each segment
         j2 = 2
      else:
         j2 = j1 + n_well_segments
      segments = pl.arange(j1, j2)
      m = pl.zeros((n_tsteps, len(segments)))
      for tstep in range(n_tsteps):
            kw = ef.iget_named_kw(varnm, wellid+nwells*tstep)
            for j in range(len(segments)):
               m[tstep, j] = kw[segments[j]-1]
      l.append(m)
   return l

class RFT_data():
   '''
   see documentation on read_rft_file.
   use get(varnm) to obtain data.
   use varnames() to get a list of available variables.
   note special methods time() and measured_depth() to get time and measured depth, respectively.
   we can decide if we will use values per segment, per valve or per meter. data is stored per segment.
   note: this ability to plot per segment/valve/meter is the main motiviation for having this piece of code.
         if the segments lengths do no vary, one could use the summary data (SOFR, SWCT etc.). one advantage
         with the summary file is that it is written continuously, whereas the RFT file is dumped at the end
         of the simulation.
   '''
   PR_SEGM  = 'pr_segm'
   PR_METER = 'pr_meter'
   PR_JOINT = 'pr_joint'
   PR_VALVE = 'pr_valve'
   def __init__(self, fname, wells, vpjs=None, sel_well=None, mode=PR_SEGM):
      '''
      vpjs is valves pr. joint pr. well. could be constant pr. well, or a vector.
      sel_well is the "selected well" (useful if only 1 well is of interest)
      mode: data is PR_SEGM, PR_METER, PR_JOINT, or PR_VALVE. Default is PR_SEGM.
      we operate with a list of well names. this list *must* have the order as defined in WELSPECS
      
      example:
      r = ECL.RFT_data('A1.RFT', ['P1', 'P2'], vpjs=[1., 1.], sel_well='P1')
      r.contours('SEGORAT')            # this is per segment
      r.mode = ECL.RFT_data.PR_METER
      r.contours('SEGORAT')            # this is per meter
      '''
      self._ef      = ecl.EclFile(fname)
      self.wells    = wells              # list of well names. *must* have the order given in WELSPEC
      self.nwells   = len(wells)
      self.vpjs     = vpjs               # valves pr. joint (for each well). could be None.
      self.sel_well = sel_well           # selected well. could be None.
      self.mode     = mode               # could be PR_SEGM, PR_METER, PR_JOINT, or PR_VALVE. 
      self._data    = {}                 # container for caching results
#
   def _wellid(self, well):
      i = 0
      for w in self.wells:
         if w == well: break
         i += 1
      return i
#
   def get(self, varnm, well=None):
      '''
      we obtain data pr. segment and convert afterwards - if requested
      well: self.sel_well is used if not given
      '''
      if not varnm in self.varnames():
         print "no such variable: ", varnm
         return None
      if well is None: well = self.sel_well
      wellid = self._wellid(well)  # convinent
      nm = '%s:%s' % (well, varnm)
      if not self._data.has_key(nm):
         self._data[nm] = read_rft_file(varnm, self._ef, nwells=self.nwells)[wellid]
      m = self._data[nm]           # matrix with data
      l = ['CONORAT', 'SEGGRAT', 'CONWRAT', 'SEGORAT', 'CONGRAT', 'SEGWRAT', 'CONPRES',
           'SEGPRES', 'PRESSURE', 'SEGGVEL', 'SEGOVEL', 'SEGWVEL'] # to be continued...
      if self.mode == RFT_data.PR_SEGM or varnm not in l: return m # dont scale these
      if self.mode == RFT_data.PR_METER                 : return m / self.segm_length(well)
      if self.mode == RFT_data.PR_JOINT                 : return m / self.segm_length(well) * JOINT_LENGTH
      if self.mode == RFT_data.PR_VALVE                 : return m / self.segm_length(well) * (JOINT_LENGTH/pl.array(self.vpjs))
#
   def time(self):
      '''
      get time for result
      '''
      return self.get('TIME')[:,0]
#
   def well_segments(self, well=None):
      if well is None: well = self.sel_well
      n_well_segments = len(self._ef.iget_named_kw('CONKH', 0))
      return pl.arange(n_well_segments) + 1
#
   def measured_depth(self, well=None):
      '''
      get measured depth (starting at 0) for given well
      well: self.sel_well is used if not given
      note: will fail if CONLENEN is not defined
      '''
      if not 'CONLENEN' in self.varnames(): return None
      if well is None: well = self.sel_well
      md = self.get('CONLENEN', well)[0]
      return md - md[0]
#
   def md(self, well=None):
      '''
      just an alias for measured_depth
      '''
      return self.measured_depth(well)
#
   def segm_length(self, well=None):
      '''
      well: self.sel_well is used if not given
      note: will fail if CONLENEN and CONLENST is not defined
      '''
      if not ('CONLENEN' in self.varnames() and 'CONLENST' in self.varnames()): return None
      if well is None: well = self.sel_well
      return self.get('CONLENEN', well)[0] - self.get('CONLENST', well)[0]
#
   def varnames(self):
      return self._ef.keys()
#
   def xplot(self, varnm1, varnm2, tsteps=None, well=None, linestyle='-', marker=None):
      '''
      cross-plot 2 variables different time steps.
      note! if tsteps is None, all time steps will be used.
      well: self.sel_well is used if not given
      '''
      if well is None: well = self.sel_well
      m1 = self.get(varnm1, well)
      m2 = self.get(varnm2, well)
      pl.figure()
      t  = self.time()
      if tsteps is None: tsteps = pl.arange(len(t))
      for tstep in tsteps:
         pl.plot(m1[tstep,:], m2[tstep,:], label='%i days'%t[tstep], marker=marker, linestyle=linestyle)
      pl.legend(loc='best')
      pl.grid(True)
      pl.xlabel('%s (%s)' % (varnm1, self.mode))
      pl.ylabel('%s (%s)' % (varnm2, self.mode))
      pl.title('%s - %s' % (UT.basename(self._ef.name), well))
      pl.show()
#
   def plot(self, varnm, tsteps=None, well=None, linestyle='-', marker=None):
      '''
      plot profiles along well for different time steps.
      note! if tsteps is None, all time steps will be used.
      well: self.sel_well is used if not given
      '''
      if well is None: well = self.sel_well
      m = self.get(varnm, well)
      pl.figure()
      md = self.md(well)
      if md is None: md = self.well_segments(well)
      t  = self.time()
      if tsteps is None: tsteps = pl.arange(len(t))
      for tstep in tsteps:
         pl.plot(md, m[tstep,:], label='%i days'%t[tstep], marker=marker, linestyle=linestyle)
      pl.legend(loc='best')
      pl.grid(True)
      pl.ylabel('%s (%s)' % (varnm, self.mode))
      pl.xlabel(self.xlabel())
      pl.title('%s - %s' % (UT.basename(self._ef.name), well))
      pl.show()
#
   def xlabel(self):
      if not 'CONLENEN' in self.varnames(): return 'segment #'
      else                                : return 'md [m]'
   def contours(self, varnm, well=None, zmin=None, zmax=None):
      '''
      well: self.sel_well is used if not given
      '''
      if well is None: well = self.sel_well
      m = self.get(varnm, well)
      md = self.md(well)
      if md is None: md = self.well_segments(well)
      t = self.time()
      # sometimes the data for t=0 is not included
      delta = len(t) - m.shape[0]
      if delta > 0: t = t[delta:]
      ax = PU.contourf(md, t, m, zmin=zmin, zmax=zmax)
      titl = '%s (%s). Case %s Well %s' % (varnm, self.mode, UT.basename(self._ef.name), well)
      pl.matplotlib.artist.setp(ax, title=titl, ylabel='TIME [days]', xlabel=self.xlabel())
      pl.show()
      return ax

class TimeSeries() :
   '''
   simple container for keeping and accessing time series, typically found in summary file
   '''
   def __init__(self, data, varnames, units, nm):
      '''
      data    : matrix of data (time series). each column is a time serie. the name of the variable
                is found in varnames
      varnames: list of variable names. order *must* be the same as in data
      units   : list of variable units. order *must* be the same as in data
      nm      : name of this time serie
      '''
      self._data = data
      self._map       = {} # keeps track of which unit goes with each variable
      self.nm         = nm
      self.wells      = [] # list of wells
      self._segms      = {} # list of all segments for each well
      self._icd_segms  = {} # list of icd-segments for each well
      self._well_segms = {} # list of well-segments for each well
      i = -1
      for varnm in varnames:
         i += 1
         self._map[varnm.rstrip()] = (i, units[i].rstrip())
         if varnm.startswith('W'):
            wname = varnm.split(SEPARATOR)[1]
            if not wname in self.wells: self.wells.append(wname)
         if varnm.startswith('S'):
            if varnm == 'STEPTYPE': continue
            if '%s-'%SEPARATOR in varnm: continue   # avoid negative segments
            wname, segm = varnm.split(SEPARATOR)[1:]
            segm = int(segm)
            if not self._segms.has_key(wname): self._segms[wname] = []
            if segm <= 1: continue  # should not be included. starting at 2.
            if not segm in self._segms[wname]: self._segms[wname].append(segm)
      # for convinience
      if not 'TIME' in varnames:
         self.time = self._data[:,self._map['YEARS'][0]] * 365.25
      else:
         self.time = self._data[:,self._map['TIME'][0]]
      self.wells.sort()
      for wname in self._segms.keys():
         self._segms[wname].sort()
         # assume well-segments are the first half - then icd-segments
         nsegms = len(self._segms[wname])
         self._well_segms[wname] = self._segms[wname][:nsegms/2]
         self._icd_segms[wname]  = self._segms[wname][nsegms/2:]
#
   def segments(self, wellnm):
      return pl.array(self._segms[wellnm])
#
   def icd_segments(self, wellnm):
      return pl.array(self._icd_segms[wellnm])
#
   def well_segments(self, wellnm):
      return pl.array(self._well_segms[wellnm])
#
   def unit(self, varnm) :
      if   varnm in self._map.keys(): return self._map[varnm][1]
      elif varnm in ('TIME', 'DAYS'): return 'days'
      else                          : return '-'
#
   def get(self, varnm, cumul=False, scaler=1., Bo=1, Bg=1, Bw=1, Rs=1) :
      '''
      Get the timeserie
      '''
      if 'SLFR' in varnm:           # total liquid production for well segments (not available in Eclipse)
         postfix = varnm.split('SLFR')[1] # f.ex. SLFR-Q21-17 -> -Q21-17
         y = self._data[:,self._map['SOFR'+postfix][0]] \
           + self._data[:,self._map['SWFR'+postfix][0]]
         self._map['SLFR'+postfix] = self._map['SWFR'+postfix]
      elif 'SFFR' in varnm:           # total flow production for well segments (not available in Eclipse)
         postfix = varnm.split('SFFR')[1] # f.ex. SLFR-Q21-17 -> -Q21-17
         y = self._data[:,self._map['SOFR'+postfix][0]].copy()*Bo # gonna do +=, so need a copy()
         if self._map.has_key('SGFR'+postfix):
            free_gas = self._data[:,self._map['SGFR'+postfix][0]] -y*Rs
            y += free_gas*Bg
         if self._map.has_key('SWFR'+postfix):
            y += self._data[:,self._map['SWFR'+postfix][0]]*Bw
         self._map['SFFR'+postfix] = self._map['SOFR'+postfix]
      elif varnm.startswith('FSN'): # sigurds number: cumulative oil / cumulative water or gas
         if varnm.endswith('W'):
            y = self._data[:,self._map['FOPT'][0]] / self._data[:,self._map['FWPT'][0]] # FSNW: water
         else:
            y = self._data[:,self._map['FOPT'][0]] / self._data[:,self._map['FGPT'][0]] # FSNG: gas
         self._map[varnm] = (-1, '-')
      elif varnm == 'TIME':
         y = self.time    # could be 'YEARS'...
      else:
         y = self._data[:,self._map[varnm][0]]
      if cumul: return scaler*UT.cumulative(self.time, y)
      else    : return scaler*y
#
   def varnames(self):
      '''
      Get list of variable names
      '''
      return self._map.keys()
#
   def has_varname(self, varnm):
      '''
      Get list of variable names
      '''
      return self._map.has_key(varnm)
#
   def cross_plot(self, varnm1, varnm2, newfig=False, do_clf=False, scaler1=1., scaler2=1., marker=None):
      '''
      cross-plot two variables
      '''
      y1 = self.get(varnm1, scaler=scaler1)
      y2 = self.get(varnm2, scaler=scaler2)
      if newfig : pl.figure()
      elif do_clf: pl.clf()
      pl.plot(y1, y2, label=self.nm, marker=marker)
      pl.xlabel('%s [%s]' % (varnm1, self.unit(varnm1)))
      pl.ylabel('%s [%s]' % (varnm2, self.unit(varnm2)))
      pl.legend(loc='best')
      pl.grid(True)
#
   def plot(self, varnm, newfig=False, do_clf=False, scaler=1., marker=None):
      '''
      Plot the timeserie
      '''
      self.cross_plot('TIME', varnm, newfig=newfig, do_clf=do_clf, scaler2=scaler, marker=marker)
#
   def savetxt(self, var_list=None, scaler=1.):
      '''always writes time first. saves a txt-file'''
      m = self.time
      if var_list is None: var_list = self.varnames()
      for var in var_list:
         m = pl.vstack((m, self.get(var, scaler=scaler)))
      fname = self.nm+'.txt'
      pl.savetxt(fname, m.T)
      print 'saved data to file %s' % fname
#
   def get_segm_data(self, varnm, wellnm, segments, scaler=1.,Bo=1., Bg=1., Bw=1., Rs=1., tsteps=None):
      d = None
      for segm in segments:
         varname = '%s%s%s%s%i' % (varnm, SEPARATOR, wellnm, SEPARATOR, segm)
         if d is None: d = self.get(varname, scaler=scaler, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs)
         else        : d = pl.vstack((d, self.get(varname, scaler=scaler, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs)))
      if not tsteps is None: return d[:, tsteps]
      else                 : return d
#
   def plot_segm_data(self, varnm, wellnm, segments, tsteps=None, newfig=False):
      if newfig: pl.figure()
      pl.plot(self.get_segm_data(varnm, wellnm, segments, tsteps=tsteps))
      if newfig:
         pl.title(self.nm)
         pl.ylabel(varnm)
         pl.xlabel('segments [-]')
#
   def contour_plot(self, varnm, wellnm, segments,
                    zmin=None, zmax=None, scaler=1., Bo=1., Bg=1., Bw=1., Rs=1.,
                    accum=False, relative=False):
      z = self.get_segm_data(varnm, wellnm, segments, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs, scaler=scaler)
      if relative:
         z /= max(z.flatten())
      if accum:
         for i in range(len(segments)):
            z[i,:] = UT.cumulative(self.time, z[i,:])
      ax =  PU.contourf(segments, self.time, z.T, zmin=zmin, zmax=zmax)
      titl = '%s - %s - %s' % (self.nm, wellnm, varnm)
      if accum: titl += ' (accum.)'
      pl.matplotlib.artist.setp(ax, title=titl, ylabel='TIME [days]', xlabel='segment # [-]')
      return ax
#
   def get_interp_value(self, varnm, t0):
      '''
      interpolate to given time.
      '''
      t = self.time
      y = self.get(varnm)
      return pl.interp([t0], t, y)[0]
#
   def aver_valve_dp(self, wellnm, segments):
       return pl.mean(self.get_segm_data('SPRD', wellnm, segments).ravel())
#
   def npv(self, oilprice, watprice, gasprice, opex, int_rate):
      '''
      npv = net present value
      based on NPV calculations by Heriot-Watt. see mail from eltazy 1/12-2015.
      prices are in $/bbl
      opex uses only liquid rates (not gas rate).
      watprice and opex are typically negative.
      prices could be scalars or vectors.
      int_rate is interest rate in % (tyically 8).
      '''
      dt = pl.diff(self.time)
      t  = self.time[1:]
      qo = self.get('FOPR')[1:]*dt
      qw = self.get('FWPR')[1:]*dt
      qg = self.get('FGPR')[1:]*dt
      revenue = qo*oilprice + qw*watprice + qg*gasprice + (qo+qw)*opex
      c = 1. if not 'm3' in self.unit('FOPR') else 1./U.BARREL         # unit conversion
      npv = c*pl.cumsum(UT.npv_scale(revenue, t, int_rate))
      return pl.concatenate(([0], npv))                                # pad a zero for t=0

class SummaryVectors() :
   '''
   simple container for accessing eclipse summary vectors.
   should at some point replace TimeSeries. this one is reading vectors as they are used -
   do not need to load everything initially
   TODO: get block names right etc (now using standard eclipse notation)
   '''
   def __init__(self, nm, separator=':'):
      '''
      nm      : name of this case
      '''
      self.nm      = nm.split('.')[0]  # remove extension
      self.shortnm = UT.basename(nm)
      self.wells       = [] # list of wells
      self._segms      = {} # list of all segments for each well
      self._icd_segms  = {} # list of icd-segments for each well
      self._well_segms = {} # list of well-segments for each well
      self.separator   = separator
      self.sum = ecl.EclSum(nm)
      if 'TIME' in self.sum.keys():
         self.time = self.sum.get_vector('TIME').values
      else:
         self.time = self.sum.get_vector('YEARS').values*365.25
      for varnm in self.sum.keys():
         if varnm.startswith('W'):
            wname = varnm.split(separator)[1]
            if not wname in self.wells: self.wells.append(wname)
         if varnm.startswith('S'):
            if varnm == 'STEPTYPE': continue
            if '%s-'%separator in varnm: continue   # avoid negative segments
            wname, segm = varnm.split(separator)[1:]
            segm = int(segm)
            if not self._segms.has_key(wname): self._segms[wname] = []
            if segm <= 1: continue  # should not be included. starting at 2.
            if not segm in self._segms[wname]: self._segms[wname].append(segm)
      # for convinience
      if not 'TIME' in self.sum.keys():
         self.time = self.sum.get_vector('YEARS').values * 365.25
      else:
         self.time = self.sum.get_vector('TIME').values
      self.wells.sort()
      for wname in self._segms.keys():
         self._segms[wname].sort()
         # assume well-segments are the first half - then icd-segments
         nsegms = len(self._segms[wname])
         self._well_segms[wname] = self._segms[wname][:nsegms/2]
         self._icd_segms[wname]  = self._segms[wname][nsegms/2:]
#
   def segments(self, wellnm):
      return pl.array(self._segms[wellnm])
#
   def icd_segments(self, wellnm):
      return pl.array(self._icd_segms[wellnm])
#
   def well_segments(self, wellnm):
      return pl.array(self._well_segms[wellnm])
#
   def unit(self, varnm) :
      if varnm in self.sum.keys(): return self.sum.unit(varnm)
      else                        : return '-'
#
   def get(self, varnm, cumul=False, scaler=1., Bo=1, Bg=1, Bw=1, Rs=1) :
      '''
      Get the vector
      '''
      if not self.sum.has_key(varnm):
         print 'WARNING: no such varnm:', varnm
         return 0
      if 'SLFR' in varnm:           # total liquid production for well segments (not available in Eclipse)
         postfix = varnm.split('SLFR')[1] # f.ex. SLFR:Q21:17 -> :Q21:17
         y = self.sum.get_vector('SOFR'+postfix).values \
           + self.sum.get_vector('SWFR'+postfix).values
      elif 'SFFR' in varnm:           # total flow production for well segments (not available in Eclipse)
         postfix = varnm.split('SFFR')[1] # f.ex. SLFR:Q21:17 -> -Q21:17
         y = self.sum.get_vector('SOFR'+postfix).values.copy()*Bo # gonna do +=, so need a copy()
         if self.sum.has_key('SGFR'+postfix):
            free_gas = self.sum.get_vector('SGFR'+postfix).values -y*Rs
            y += free_gas*Bg
         if self.sum.has_key('SWFR'+postfix):
            y += self.sum.get_vector('SWFR'+postfix).values*Bw
      elif varnm.startswith('FSN'): # sigurds number: cumulative oil / cumulative water or gas
         if varnm.endswith('W'):
            y = self.sum.get_vector('FOPT').values / self.sum.get_vector('FWPT').values # FSNW: water
         else:
            y = self.sum.get_vector('FOPT').values / self.sum.get_vector('FGPT').values # FSNG: gas
      elif varnm == 'TIME':
         y = self.time    # could be 'YEARS'...
      else:
         y = self.sum.get_vector(varnm).values
      if cumul: return scaler*UT.cumulative(self.time, y)
      else    : return scaler*y
#
   def varnames(self):
      '''
      Get list of vector names
      '''
      return list(self.sum.keys())
#
   def has_varname(self, varnm):
      '''
      Check if vector is available
      '''
      return self.sum.has_key(varnm)
#
   def cross_plot(self, varnm1, varnm2, newfig=True, do_clf=False, scaler1=1., scaler2=1., marker=None):
      '''
      cross-plot two variables
      '''
      y1 = self.get(varnm1, scaler=scaler1)
      y2 = self.get(varnm2, scaler=scaler2)
      if newfig : pl.figure()
      elif do_clf: pl.clf()
      pl.plot(y1, y2, label=self.shortnm, marker=marker)
      pl.xlabel('%s [%s]' % (varnm1, self.unit(varnm1)))
      pl.ylabel('%s [%s]' % (varnm2, self.unit(varnm2)))
      pl.legend(loc='best')
      pl.grid(True)
#
   def plot(self, varnm, newfig=True, do_clf=False, scaler=1., marker=None):
      '''
      Plot the vector
      '''
      self.cross_plot('TIME', varnm, newfig=newfig, do_clf=do_clf, scaler2=scaler, marker=marker)
#
   def savetxt(self, var_list=None, scaler=1.):
      '''always writes time first. saves a txt-file'''
      m = self.time
      if var_list is None: var_list = self.varnames()
      for var in var_list:
         m = pl.vstack((m, self.get(var, scaler=scaler)))
      fname = self.shortnm+'.txt'
      pl.savetxt(fname, m.T)
      print 'saved data to file %s' % fname
#
   def get_segm_data(self, varnm, wellnm, segments, scaler=1.,Bo=1., Bg=1., Bw=1., Rs=1., tsteps=None):
      d = None
      for segm in segments:
         varname = '%s%s%s%s%i' % (varnm, self.separator, wellnm, self.separator, segm)
         if d is None: d = self.get(varname, scaler=scaler, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs)
         else        : d = pl.vstack((d, self.get(varname, scaler=scaler, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs)))
      if not tsteps is None: return d[:, tsteps]
      else                 : return d
#
   def plot_segm_data(self, varnm, wellnm, segments, tsteps=None, newfig=True):
      if newfig: pl.figure()
      pl.plot(self.get_segm_data(varnm, wellnm, segments, tsteps=tsteps))
      if newfig:
         pl.title(self.shortnm)
         pl.ylabel(varnm)
         pl.xlabel('segments [-]')
#
   def get_block_data(self, varnm, blocks, scaler=1.,Bo=1., Bg=1., Bw=1., Rs=1., tsteps=None):
      d = None
      for block in blocks:
         varname = '%s%s%i,%i,%i' % (varnm, self.separator, block[0], block[1], block[2])
         if d is None: d = self.get(varname, scaler=scaler, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs)
         else        : d = pl.vstack((d, self.get(varname, scaler=scaler, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs)))
      if not tsteps is None: return d[:, tsteps]
      else                 : return d
#
   def plot_block_data(self, varnm, blocks, tsteps=None, newfig=True):
      if newfig: pl.figure()
      d = self.get_block_data(varnm, blocks, tsteps=tsteps)
      ntsteps = d.shape[1]
      x = pl.arange(len(blocks)) + 1  # just referring to blocks by indices (1,2,...)
      for i in range(ntsteps):
         y = d[:,i]
         if ntsteps > len(UT.COLOURS):
            red = float(i) / ntsteps
            col = (red,0,1-red)
         else:
            col = UT.COLOURS[i]
         pl.plot(x, y, color=col)
      if newfig:
         pl.title(self.shortnm)
         pl.ylabel(varnm)
         pl.xlabel('blocks [-]')
         pl.grid(True)
#
   def contour_plot_blocks(self, varnm, blocks,
                    zmin=None, zmax=None, scaler=1., Bo=1., Bg=1., Bw=1., Rs=1.,
                    accum=False, relative=False):
      z = self.get_block_data(varnm, blocks, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs, scaler=scaler)
      if relative:
         z /= max(z.flatten())
      if accum:
         for i in range(len(blocks)):
            z[i,:] = UT.cumulative(self.time, z[i,:])
      ax =  PU.contourf(range(len(blocks)), self.time, z.T, zmin=zmin, zmax=zmax)
      titl = '%s - %s' % (self.shortnm, varnm)
      if accum: titl += ' (accum.)'
      pl.matplotlib.artist.setp(ax, title=titl, ylabel='TIME [days]', xlabel='block # [-]')
      return ax
#
   def contour_plot(self, varnm, wellnm, segments,
                    zmin=None, zmax=None, scaler=1., Bo=1., Bg=1., Bw=1., Rs=1.,
                    accum=False, relative=False):
      z = self.get_segm_data(varnm, wellnm, segments, Bo=Bo, Bg=Bg, Bw=Bw, Rs=Rs, scaler=scaler)
      if relative:
         z /= max(z.flatten())
      if accum:
         for i in range(len(segments)):
            z[i,:] = UT.cumulative(self.time, z[i,:])
      ax =  PU.contourf(segments, self.time, z.T, zmin=zmin, zmax=zmax)
      titl = '%s - %s - %s' % (self.shortnm, wellnm, varnm)
      if accum: titl += ' (accum.)'
      pl.matplotlib.artist.setp(ax, title=titl, ylabel='TIME [days]', xlabel='segment # [-]')
      return ax
#
   def get_interp_value(self, varnm, t0):
      '''
      interpolate to given time.
      '''
      t = self.time
      y = self.get(varnm)
      return pl.interp([t0], t, y)[0]
#
   def aver_valve_dp(self, wellnm, segments):
       return pl.mean(self.get_segm_data('SPRD', wellnm, segments).ravel())
#
   def npv(self, oilprice, watprice, gasprice, opex, int_rate):
      '''
      npv = net present value
      based on NPV calculations by Heriot-Watt. see mail from eltazy 1/12-2015.
      prices are in $/bbl
      opex uses only liquid rates (not gas rate).
      watprice and opex are typically negative.
      prices could be scalars or vectors.
      int_rate is interest rate in % (tyically 8).
      '''
      dt = pl.diff(self.time)
      t  = self.time[1:]
      qo = self.get('FOPR')[1:]*dt
      qw = self.get('FWPR')[1:]*dt
      qg = self.get('FGPR')[1:]*dt
      revenue = qo*oilprice + qw*watprice + qg*gasprice + (qo+qw)*opex
      c = 1. if not 'm3' in self.unit('FOPR') else 1./U.BARREL         # unit conversion
      npv = c*pl.cumsum(UT.npv_scale(revenue, t, int_rate))
      return pl.concatenate(([0], npv))                                # pad a zero for t=0

def _build_varnames(vnames, wellnms, rec1, rec2):
   varnames = []
   for (varnm, wnm, r1, r2) in zip(vnames, wellnms, rec1, rec2):
      varnm = varnm.strip()
      wnm   = wnm.strip()
      r1 = r1.strip()
      r2 = r2.strip()
      if len(wnm)   > 0: varnm += SEPARATOR + wnm
      if len(r1) > 0:
         if ' ' in r1:
            rec = r1.split()
            r1 = ','.join(rec)
         varnm += SEPARATOR + r1
      if len(r2) > 0:
         if ' ' in r2:
            rec = r2.split()
            r2 = ','.join(rec)
         varnm += SEPARATOR + r2
      varnames.append(varnm)
   return varnames

def read_summary(sumryfile):
   '''
   reads summary file, returning a TimeSeries container with all data, just like read_rsm (below).
   uses ert.
   consider using read_summary2 in stead...
   '''
   casenm = UT.basename(sumryfile)
   varnames = []
   units    = []
   data = None
   sum = ecl.EclSum(casenm)
   for varnm in sum.keys():
      y = sum.get_vector(varnm).values
      if data is None: data = y
      else           : data = pl.vstack((data, y))
      # dont know where to find units. hardcoded!!!
      if   varnm[-2:] in ('PR'): unit = 'Sm3/d'
      elif varnm[-2:] in ('PT'): unit = 'Sm3'
      elif 'BHP' in varnm      : unit = 'bar'
      elif varnm == 'TIME'     : unit = 'days'
      else                     : unit = '-'
      units.append(unit)
      varnames.append(varnm.replace(':',SEPARATOR)) # same format as read_rsm uses (for segment data)
   return TimeSeries(data.T, varnames, units, casenm)

def read_summary2(sumryfile):
   '''
   reads summary file, returning a SummaryVectors container
   NB! using ECL.SummaryVectors. Not thoroughly debugged!
   '''
   return SummaryVectors(sumryfile)

def read_rsm(fname, t_start=-pl.Inf, t_end=pl.Inf, dt=-1.) :
   '''
   reads rsm file, returning a TimeSeries container with all data
   rsm-file is divided into sections of 10 (currently) variables that are written together as columns.
   this function writes these columns to a tempfile which is then read by loadtxt (a bit cumbersome...)
   segment-data has variable names given by <varnm>-<wellnm>-<segmno>
   '''
   rsmfile = open(fname)
   tmpnm   = tempfile.mkstemp()[1]
   _log.debug('tmpfile= %s' % tmpnm)
   #
   n = 0
   section  = 0
   # time is always first column of each section
   # and must be treated specially
   varnames = ['TIME']
   units    = ['DAYS']
   t        = -pl.Inf
   t_prev   = 0
   data = None
   # parse file
   while True:
      n += 1
      line = rsmfile.readline()
      if not line: break
      if len(line) <= 1: continue
      # new section ?
      if line[1:15] == 'SUMMARY OF RUN' :
         if section > 0 :
            tmpfile.close()
            d = pl.loadtxt(tmpnm)
            if section == 1 : data = d
            else            : data = pl.hstack((data, d[:,1:]))
            t_prev = 0
         tmpfile = open(tmpnm, 'w')
         section += 1
         #
         # read header
         line = rsmfile.readline()
         vnames = line.split('\t')[2:]
         line = rsmfile.readline()
         units.extend(line.split('\t')[2:])
         line = rsmfile.readline()
         if '**' in line: line = rsmfile.readline() # check if multiplier
         wellnms = line.split('\t')[2:]
         line = rsmfile.readline()
         rec1 = line.split('\t')[2:]
         if pl.any([vname.startswith('L') for vname in vnames]):
            line = rsmfile.readline()
            rec2 = line.split('\t')[2:]
         else:
            rec2 = ['' for vname in vnames] # list of empty strings
         varnames.extend(_build_varnames(vnames, wellnms, rec1, rec2))
         continue
      # write data to tempfile (unless empty line)
      if len(line.rstrip()) <= 1: continue
      t = float(line.split()[0])
      if t_start <= t <= t_end and t - t_prev > dt:
         tmpfile.write(line)
         t_prev = t
   #
   # remember saving last data
   tmpfile.close()
   d = pl.loadtxt(tmpnm)
   if not data is None: data = pl.hstack((data, d[:,1:]))
   else               : data = d
   #
   os.unlink(tmpnm)
   rsmfile.close()
   #
   return TimeSeries(data, varnames, units, fname.split('.')[0])

def read_rates_from_prtfile(fname, t_start=-pl.Inf, t_end=pl.Inf):
   f = open(fname)
   t = []
   oilr = []
   gasr = []
   watr = []
   cnt = pl.Inf
   for line in f:
      cnt += 1
      if ' TIME= ' in line:
         rec = line.split()
         t_ = float(rec[3])
         if t_ < t_start: continue
         if t_ > t_end  : break
         t.append(t_)
         cnt = 0
         continue
      if cnt == 5:
         rec = line.split()
         oilr.append(float(rec[1]))
         continue
      if cnt == 6:
         rec = line.split()
         watr.append(float(rec[1]))
         continue
      if cnt == 7:
         rec = line.split()
         gasr.append(float(rec[1]))
         continue
   f.close()
   data = pl.vstack((pl.array(t), pl.array(oilr), pl.array(watr), pl.array(gasr))).T
   varnames = ('TIME', 'OILR', 'WATR', 'GASR')
   units = ('days', 'm3/day', 'm3/day', 'm3/day')
   return TimeSeries(data, varnames, units, fname.split('.')[0])

def _read_3D_data(casenm, nx, ny, nz, file_ext):
   '''
   uses Ert
   '''
   fname = UT.glob('%s.*%s' % (casenm, file_ext))[0] # handle both .FUNRTS and .UNRST
   rst = ecl.EclFile(fname)
   ncells = nx*ny*nz
   varnames = set()
   for head in rst.headers:
      if head[1] == ncells: varnames.add(head[0])
   d = {} # gonna collect all data here
   for varname in varnames:
      d[varname] = []
      for tstep in range(rst.num_named_kw(varname)):
         data = pl.zeros((nx,ny,nz))
         rstd = rst.iget_named_kw(varname, tstep)
         n = -1
         for k in range(nz):
            for j in range(ny):
               for i in range(nx):
                  n += 1
                  data[i,j,k] = rstd[n]
         d[varname].append(data)
   d['dates'] = rst.dates[:]
   return d

def _read_3D_data_old(casenm, nx, ny, nz, file_ext):
   '''
   reads ascii file. pure python. a quick test using iPython's %timeit, showed
   that this routine 3 times as fast as the new routine for formatted result-file
   (and twice as fast as the new for a unformatted file...)
   '''
   fname = casenm+'.'+file_ext
   print 'reading', fname
   f = open(fname)
   alldata  = {}
   ncells = nx*ny*nz
   searchstr = "%i 'REAL'" % ncells
   readthis = False
   n = j = k = 0
   i = -1
   for line in f:
      if searchstr in line:
         varname = line.split()[0].replace("'","")
         data = pl.zeros((nx,ny,nz))
         if not alldata.has_key(varname): alldata[varname] = []
         alldata[varname].append(data)
         readthis = True
         continue
      if not readthis: continue
      rec = line.split()
      for num in rec:
         i += 1
         n += 1
         if i == nx:
            i = 0
            j += 1
            if j == ny:
               j = 0
               k += 1
         data[i,j,k] = float(num)
      if n == ncells:
         n = j = k = 0
         i = -1
         readthis = False
   return alldata

class Three_D_Data():
   def __init__(self, casenm, nx, ny, nz, file_ext='UNRST'):
      self.casenm = casenm
      d = _read_3D_data(casenm, nx, ny, nz, file_ext=file_ext)
      self.varnames   = d.keys()
      self.data       = d
      self.dates      = d['dates']
      self.ntimesteps = len(self.dates)
      self.time       = pl.zeros(self.ntimesteps)
      self.nx         = nx
      self.ny         = ny
      self.nz         = nz
      for i in range(1, self.ntimesteps):
         self.time[i] = (self.dates[i] - self.dates[0]).days
#
   def get(self, varnm, tstep):
      return self.data[varnm][tstep]
#
   def dump(self):
      f = open(self.casenm+'.pck', 'w')
      cPickle.dump(self, f)
      f.close()
#
   @staticmethod
   def load(casenm):
      f = open(casenm+'.pck')
      d = cPickle.load(f)
      f.close()
      return d
#
   def contour_plot(self, varnm, direction, imin, imax, jmin, jmax, kmin, kmax, dgrid=None):
      '''
      direction should be 'X', 'Y', or 'Z'
      if 'X' is given, it plots the time-development along the line
      given by imin to imax at jmin and kmin. similarly for 'Y' and 'Z'
      dgrid is used if the grid is not regular. then this should be the grid size
      for each grid cell for the relevant direction.
      '''
      if   direction.upper() == 'X':
         ii = pl.arange(imin, imax)
         ij = jmin
         ik = kmin
         n = len(ii)
         igrid = ii
      elif direction.upper() == 'Y':
         ii = imin
         ij = pl.arange(jmin, jmax)
         ik = kmin
         n = len(ij)
         igrid = ij
      else:
         ii = imin
         ij = jmin
         # want positive k downwards
         ik = pl.arange(kmax-1, kmin-1, -1)
         n = len(ik)
         igrid = ik
      d = self.data[varnm]
      z = pl.zeros((n, len(d))); tstep = 0
      for m in d:
         z[:,tstep] = m[ii, ij, ik]
         tstep += 1
      if dgrid is None:
         grid_coord = pl.arange(z.shape[0])
      else:
         grid_coord = pl.cumsum(dgrid)[igrid]
      ax =  PU.contourf(grid_coord, pl.arange(z.shape[1]), z.T)
      titl = '%s - %s' % (self.casenm, varnm)
      pl.matplotlib.artist.setp(ax, title=titl, ylabel='Time [tstep]', xlabel=direction)

def _read_inflow_data(casenm, segm_start, segm_end):
   nvars     = 6
   n_segm = segm_end - segm_start + 1
   data = pl.zeros([n_segm, nvars])
   potential_solving = False
   f = open(casenm+'.DBG')
   t = t_prev = -1.
   datalist = []; times = []
   #
   # read file
   for line in f :
      if ' 0--TIME STEP STARTING' in line:
         rec = line.split()
         t_prev = t
         t = float(rec[7])
         potential_solving = False
         if t_prev > 0:
            datalist.append(data.copy())
            times.append(t_prev)
      if 'CALCULATING WELL POTENTIAL' in line:
         potential_solving = True
      if potential_solving : continue
      # array filled redundantly (for each iteration), only last one is interesting
      if 'WCFRSG' in line:
         rec = line.split()
         for j in range(nvars) :
            i = int(rec[1])-segm_start 
            if i < 0 or i > n_segm-1: continue
            data[i, j] = float(rec[j+2])
   # keep data from last timestep
   datalist.append(data.copy())
   times.append(t)
   # we now have all the data. for ease of access, we'll process it
   n_tsteps = len(times)
   rho = pl.zeros([n_segm, n_tsteps])
   q   = pl.zeros([n_segm, n_tsteps])
   dp  = pl.zeros([n_segm, n_tsteps])
   mu  = pl.zeros([n_segm, n_tsteps]) # viscosity
   for n in range(n_tsteps):
      data = datalist[n]
      i = -1
      for row in data :
         # each row is data for one segment with the following format:
         # row[0] is total mass flow rate
         # row[1] is viscosity
         # row[2] is density
         # row[3] not used
         # row[4] not used
         # row[5] is dp
         i += 1
         rho[i,n] = row[2]
         q[i,n]   = row[0] / rho[i,n]
         mu[i,n]  = row[1]
         dp[i,n]  = row[5]
   return pl.array(times), rho, q, dp, mu

class InflowData() :
   '''
   simple container for accessing inflow data obtained from dbg-file.
   could be useful if we want to *know* what data eclipse uses for calculating
   dp for icd/aicd valves etc.
   '''
#
   def __init__(self, casenm, segm_start, segm_end):
      '''
      segm_start: first icd segment
      segm_end  : last icd segment
      units : given by eclipse (typially: t=days, dp=bar, q=m3/day/segment, mu=cP, rho=kg/m3)
      must have:
      -- debug info ala Fredric Spiegel
      WELDEBUG
       P1 3000 /
      /
      #
      2016-11-01: simplified it so that it only picks up the provided data. if you need conversion
                  to flowrate per valve etc, you need to do it yourself...
      '''
      self.casenm = casenm
      self.t, self.rho, self.q, self.dp, self.mu = _read_inflow_data(casenm, segm_start, segm_end)
#
   def dump(self):
      f = open(self.casenm+'.pck', 'w')
      cPickle.dump(self, f)
      f.close()
#
   def write_asciifile(self, all_vars=False):
      fname = self.casenm + '_dbg.asc'
      f = open(fname, 'w')
      f.write('#  time density viscosity flowrate dp\n')
      for n in len(self.t):
         f.write('%.1f %.2f %.2f %.1f %.2f\n' % (self.t[n], self.rho[n], self.mu[n], self.q[n], self.dp[n]))
      print '\n'
      f.close()
#
   @staticmethod
   def load(casenm):
      f = open(casenm+'.pck')
      infl = cPickle.load(f)
      f.close()
      return infl
#
   def contour_plot(self, varnm, zmin=None, zmax=None, nlvls=20):
      z = self.__dict__[varnm]
      ax = PU.contourf(self.t, self.md, z,
            zmin=zmin, zmax=zmax, nlvls=nlvls)
      titl = '%s. %s' % (self.casenm, varnm)
      pl.matplotlib.artist.setp(ax, title=titl, xlabel='Time [days]', ylabel='MD [m]')
      pl.show()
#
   def cross_plot(self, varnm1, varnm2, tstep):
      y1 = self.__dict__[varnm1]
      y2 = self.__dict__[varnm2]
      pl.figure()
      pl.plot(y1[:,tstep], y2[:,tstep], '*')
      pl.xlabel('%s' % varnm1)
      pl.ylabel('%s' % varnm2)
      pl.title('%s (time = %.1f)' % (self.casenm, tstep))
      pl.show()
#
   def varnames(self):
      return filter(lambda nm: not nm.startswith('_'), self.__dict__.keys())

def peaceman_radius(dy, dz, Ky, Kz):
   return 0.28 * pl.sqrt(dz**2*pl.sqrt(Ky/Kz) + dy**2*pl.sqrt(Kz/Ky)) / ((Ky/Kz)**.25 + (Kz/Ky)**.25)

def PI_eclipse(mu, r_drainage, r_well, h, K, theta=2*pl.pi):
   '''
   PI for one phase well.
   formula found in Chap 76 (p.1175ff) in Eclipse Technical Description 2009.2
   K is in md, all length in meters, viscosity in cP
   '''
   return (0.008527 * theta * K * h) / pl.log(r_drainage/r_well) / mu

def PI_joshi(mu, Bo, r_drainage, r_well, h, L, K):
   '''
   PI for one phase well.
   formula found in Tarek, Reservoir Engineer Handbook 3rd ed, p.535, eq 7-58
   K is in md, all length in meters, viscosity in cP
   '''
   # unit conversions m->ft
   L    = L/U.FEET
   r_drainage = r_drainage/U.FEET
   r_well  = r_well/U.FEET
   h    = h/U.FEET
   #
   a = (L/2)*pl.sqrt((0.5 + pl.sqrt(0.25+(2*r_drainage/L)**4)))
   R = (a+pl.sqrt(a**2 - (L/2.)**2)) / (L/2.)
   PI = 0.0078*h*K / (mu*Bo*(pl.log(R) + (Bo**2*h/L)*pl.log(h/(2*r_well))))
   return PI*0.1589873 # bbl/d -> m3/d

def runtime_analyzis(PRT_files):
   if type(PRT_files) is str: PRT_files = (PRT_files, ) # accept string also
   for fname in PRT_files:
      cmd = 'echo %s' % fname
      os.system(cmd)
      cmd = 'grep PROBLEM %s | wc' % fname
      os.system(cmd)
      cmd = 'grep -A1 Problems %s' % fname
      os.system(cmd)
      cmd = 'grep "TIME= " %s | tail -1' % fname
      os.system(cmd)

def lsf_run(datafiles, nCPUs=1, version='2012.2', queue='normal', interconn_type='scali', sleeptime=None):
   if type(datafiles) is str: datafiles = [datafiles]
   for datafile in datafiles:
      base = os.path.splitext(datafile)[0]
      if nCPUs > 1: opt = '-lsfhtype any -%s -procs %i' % (interconn_type, nCPUs)
      else        : opt = ''
      cmd =  '/prog/ecl/grid/macros/@eclipse -lsf -ver %s -lsfqueue %s %s -local -data `pwd` -file %s' % (version, queue, opt, base)
      print cmd
      if sleeptime: time.sleep(sleeptime)
      os.system(cmd)

def lsf_run2(datafiles, version='2012.2', max_running=0, lsf=True):
   '''
   based on tips from Joachim Hove.
   doesnt work on RH3 (python 2.5)
   '''
   if type(datafiles) is str: datafiles = [datafiles]
   if lsf: drv = ert.job_queue.driver.LSF_DRIVER
   else  : drv = ert.job_queue.driver.LOCAL_DRIVER
   que = ert.ecl.ecl.EclQueue(driver_type=drv, max_running=max_running , size=len(datafiles), ecl_version=version)
   for datafile in datafiles:
      if not datafile.startswith(DELIM): datafile = '%s/%s' % (os.getcwd(), datafile) # must have full path
      que.submit(datafile)
      print datafile
   que.block_waiting()    # que is running in separate thread. must make sure all jobs are started (i think ...)

def lsf_run3(datafiles):
   if type(datafiles) is str: datafiles = [datafiles]
   cmd = 'run_eclipse ' + ' '.join(datafiles)
   UT.tcsh(cmd)

def lcl_run(datafile, version='2016.1'):
   '''
   running single process in foreground
   one annyoing detail is the ECL.CFG file. sometimes the simulation doesnt start because it asks what
   to do with this file. i try hard to avoid this...
   '''
   os.system('rm -f ECL.CFG')
   cmd = '@eclipse -vers %s -override -file %s' % (version, UT.basename(datafile)) 
   os.system(cmd)
   os.system('rm -f ECL.CFG')

def relperm_oil(swi, sor, kr_max, expon, npoints):
   '''uses a corey exponential function. see http://jgmaas.com/scores/facts.html
   swi    : irred. water sat
   sor    : remaining oil sat
   kr_max : max relperm
   expon  : corey exponent
   npoints: number of points in output arrays
   '''
   if npoints > 4:
      sw = numpy.linspace(swi, 1.-sor, npoints-2)
      kr = kr_max * ((1.-sw - sor)/(1.-sor - swi))**expon 
      kr = numpy.concatenate((numpy.r_[1.], kr, numpy.r_[0.]))
      sw = numpy.concatenate((numpy.r_[0.], sw, numpy.r_[1.]))
   else:
      # degenerate case
      sw = numpy.r_[0., swi, (1.-sor), 1.]
      kr = numpy.r_[1., 1., kr_max, 0.]
   return (sw, kr)

def relperm_water(swi, sor, kr_max, expon, npoints):
   '''uses a corey exponential function. see http://jgmaas.com/scores/facts.html
   swi    : irred. water sat
   sor    : remaining oil sat
   kr_max : max relperm
   expon  : corey exponent
   npoints: number of points in output arrays
   '''
   if npoints > 4:
      sw = numpy.linspace(swi, 1.-sor, npoints-2)
      kr = kr_max * ((sw -swi)/(1.-sor - swi))**expon 
      kr = numpy.concatenate((numpy.r_[0.], kr, numpy.r_[1.]))
      sw = numpy.concatenate((numpy.r_[swi], sw, numpy.r_[1.]))
   else:
      # degenerate case
      sw = numpy.r_[0., swi, (1.-sor), 1.]
      kr = numpy.r_[0., 0., kr_max, 1.]
   return (sw, kr)

def create_relpermfile_oil_n_water(swi, sor, krw_max, kro_max, expon_w, expon_o, npoints, fname):
   '''uses a corey exponential function. see http://jgmaas.com/scores/facts.html
   swi    : irred. water sat
   sor    : remaining oil sat
   krw_max: max relperm for water
   kro_max: max relperm for oil
   expon_w: corey exponent for water
   expon_o: corey exponent for oil
   npoints: number of points in output arrays
   fname  : filename of output file (to be included in eclipse data file)
   '''
   f = open(fname, 'w')
   f.write('''
-- Corey relperm curve.
-- Swi     = %.3f
-- Sorw    = %.3f
-- Krw_max = %.3f
-- Kro_max = %.3f
-- w_exp   = %.3f
-- o_exp   = %.3f

SWOF
--Sw   Krw   Kro
''' % (swi, sor, krw_max, kro_max, expon_w, expon_o))
   Sw, Krwn = relperm_water(swi, sor, krw_max, expon_w, npoints)
   Sw, Kron = relperm_oil(swi, sor, kro_max, expon_o, npoints)
   for (s,kw,ko) in zip(Sw,Krwn,Kron):
      f.write(' %.3f %.3f %.3f 1*\n' % (s,kw,ko))
   f.write('/\n')
   f.close()

def analyze_run(logfiles, messagetype='PROBLEM', newfig=True, use_basename=True, **kwargs):
   import collections
   if newfig: pl.figure()
   for logfile in logfiles:
      f = open(logfile)
      t = []; cnt = collections.deque([0])
      for line in f:
         match = re.search('%s  AT'%messagetype, line)
         if match: # new problem. add it up
            cnt[-1] += 1
         match = re.search('^ STEP.*TIME=', line)
         if match: # new timestep
            rec = line.split('TIME=')[1].split()
            t.append(float(rec[0]))
            cnt.append(0)
      cnt.pop() # remove last element (it shouldnt be there)
      cumcnt = pl.cumsum(cnt)
      # plot linear or log?
      if cumcnt[-1] > 1000 : plotfnc = pl.semilogy
      else                 : plotfnc = pl.plot
      if use_basename: lbl = UT.basename(logfile)
      else           : lbl = logfile
      plotfnc(t, cumcnt, label=lbl, **kwargs)
      f.close()
      print logfile, cumcnt[-1]
   pl.xlabel('time [days]')
   pl.ylabel('# %s accumulated' % messagetype)
   pl.legend(loc='best')
   pl.grid()
   if   messagetype == 'PROBLEM':
      pl.title('Problems are usually caused by non-converging reservoir calculations')
   elif messagetype == 'WARNING':
      pl.title('Warnings are usually caused by non-converging well calculations')

def dp_along_well(segments, rho, mu, diam, q_segm, roughn=1.5e-05, dy=1.):
   '''
   uses eclipse equations to calculate dp along well.
   based on Ali's matlab script.
   segments: array of individual segment lentghts (in m)
   [rho]  = kg/m3
   [mu]   = cp
   q_segm = array of indivdual segment flows [m3/d/segm]
   [diam] = m
   '''
   C_f = 4.343e-15  # friction factor constant
   C_r = 0.01474    # conversion factor
   #
   well_length = pl.sum(segments)                 # length of the horizontal well
   ny          = int(well_length / float(dy))
   segments_   = dy*pl.ones(ny)                   # small segments for numerical integration
   q_segm_     = pl.interp(segments_, segments, q_segm/(segments*dy))
   q_acc       = pl.zeros(ny)                     # accumul flow from toe to heel
   f           = pl.zeros(ny)                     # Friction coeff. (combinations)
   dp_         = pl.zeros(ny)                     # pressure inside tubing (well)
   #
   # downscaling and calculating the friction coeff. 
   for i in range(ny):
       q_acc[i] = q_segm_[i] * (i+1)
       Re = (C_r * rho) / (mu * diam) * q_acc[i]
       if Re < 4000:        # laminar
           f[i] = 16. / Re
       else:                # turbulent
           f[i] = 1 / (3.6 * pl.log10(6.9 / Re + (roughn / (3.7 * diam)) ** (10./9.))) **2
   # pressure  from toe to heel 
   dp_[0] = ((C_f * rho/diam**5) * dy * (f[1]*q_acc[1]**2 + f[0]*q_acc[0]**2)) / 2.
   for i in pl.arange(ny-1):
       dp_[i+1] = dp_[i] + ((C_f * rho/diam**5) * dy * (f[i+1]*q_acc[i+1]**2 + f[i]*q_acc[i]**2)) / 2.
   # upscale from every meter to segment lengths
   dp = pl.zeros(segments.shape)
   ind = 0
   for i in range(len(segments)):
      ind += int(segments[i] / float(dy))
      dp[i] = dp_[ind]
   #
   return dp[::-1] # reverse it to get heel to toe

def read_reporttimes(log_fname):
   f = open(log_fname)
   reporttimes = []
   while(True):
      line1 = f.readline()
      if not line1: break
      if not 'MESSAGE  AT TIME' in line1: continue
      line2 = f.readline()
      if not 'RESTART FILE WRITTEN' in line2: continue
      reporttimes.append(float(line1.split()[3]))
   f.close()
   return pl.array(reporttimes)

def corey_relperm(Swi, Sorw, Krw_max, Kro_max, w_exp, o_exp, n_entries=20):
   '''
   might as well use relperm_water and/or relperm_oil above
   use with care: normalizes Sw - dont remember why...
   Equations and notation taken from
   http://sp-st02.statoil.com/sites/2ed91913-85c4-4f83-a610-c8178fdc1773/Pregrino/Document%20library/Peregrino%20%20LET%20relperm%20and%20Pc%20for%20Eclipse%20input.xls
   Sw        : Irred. water
   Sorw      : Residual oil
   Kro_max   : max relperm for oil
   Krw_max   : max relperm for water
   w_exp     : corey exponent for water
   o_exp     : corey exponent for oil
   n_entries : length of vectors
   '''
   Sw  = numpy.linspace(0., 1., n_entries)
   Krw = Sw**w_exp
   Kro = (1.-Sw)**o_exp
   Swn = Sw / (1.-Sorw-Swi)
   return (Swn, Krw*Krw_max, Kro*Kro_max)

def create_corey_relperm(fname, Swi, Sorw, Krw_max, Kro_max, w_exp, o_exp, n_entries=20):
   '''
   might as well use relperm_water and/or relperm_oil above
   use with care: normalizes Sw - dont remember why...
   Equations and notation taken from
   http://sp-st02.statoil.com/sites/2ed91913-85c4-4f83-a610-c8178fdc1773/Pregrino/Document%20library/Peregrino%20%20LET%20relperm%20and%20Pc%20for%20Eclipse%20input.xls
   Sw        : Irred. water
   Sorw      : Residual oil
   Kro_max   : max relperm for oil
   Krw_max   : max relperm for water
   w_exp     : corey exponent for water
   o_exp     : corey exponent for oil
   n_entries : length of vectors
   '''
   f = open(fname, 'w')
   f.write('''
-- Corey relperm curve.
-- Swi     = %.3f
-- Sorw    = %.3f
-- Krw_max = %.3f
-- Kro_max = %.3f
-- w_exp   = %.3f
-- o_exp   = %.3f

SWOF
--Sw   Krw   Kro
''' % (Swi, Sorw, Krw_max, Kro_max, w_exp, o_exp))
   Sw, Krwn, Kron = corey_relperm(Swi, Sorw, Krw_max, Kro_max, w_exp, o_exp, n_entries=n_entries)
   for (s,kw,ko) in zip(Sw,Krwn,Kron):
      f.write(' %.3f %.3f %.3f *\n' % (s,kw,ko))
   f.write('/\n')
   f.close()
   print 'file %s created' % fname
   return Sw, Krwn, Kron

def read_property(fname, prop, ncols=4):
   '''
   assumes ascii-file has structure
   BOX
      i1 i2 j1 j2 k1 k2
   /
   PROPERTY_X
      p1 p2 p3 p4
      p5 ...
      .
      .
      pNNN
   /
   note: floviz exports properties on this format.
   '''
   f = open(fname)
   while True:
      line = f.readline()
      if line.startswith('BOX'):
         # read dimensions and initialize
         line = f.readline()
         rec = line.split()
         i1 = int(rec[0]); i2 = int(rec[1])
         j1 = int(rec[2]); j2 = int(rec[3])
         k1 = int(rec[4]); k2 = int(rec[5])
         m = pl.zeros((i2-i1+1,j2-j1+1,k2-k1+1))
         continue
      if line.startswith(prop):
         # read data
         n = ncols - 1
         for k in pl.arange(k2-k1+1):
            for j in pl.arange(j2-j1+1):
               for i in pl.arange(i2-i1+1):
                  n += 1
                  if n == ncols:
                     line = f.readline()
                     rec = line.split()
                     n = 0
                  m[i,j,k] = rec[n]
         break   # done
   f.close()
   return m

def _readline(f):
   line = f.readline()
   line = line.strip()
   line = line.replace(DELIM, '')
   return line

def read_input_property(fname, prop, ni, nj, nk):
   '''
   reads input property, typically from include file.
   '''
   f = open(fname)
   m = pl.zeros((ni,nj,nk))
   while True:
      line = _readline(f)
      if not line.startswith(prop):
         continue
      # if we get here, were gonna read data
      rec = []
      k = 0
      while k < nk:
         j = 0
         while j < nj:
            i = 0
            while i < ni:
               if len(rec) == 0:
                  line = _readline(f)
                  rec = line.split()
               m[i,j,k] = rec.pop(0)
               i += 1
            j += 1
         k += 1
      break   # done
   f.close()
   return m

def read_input_property2(fname, prop, nrows=pl.Inf, stop_sign='', identif='', skip=0):
   '''
   #
   NB! consider using DataDeck instead (see below), which is much more robust.
   #
   reads input property, typically from include file.
   this one reads data as written - as a matrix
   use 'identif' to select for example only data assiosiated with a specific well
   prop : name of Eclipse keyword (like DENSITY or COMPDAT)
   note1: one of nrows or stop_sign should be given
   note2: stop_sign is used with line.startswith()
   note3: a '' matches anything
   note4: it may fail totally since it doesnt really read the records like Eclipse does
   '''
   f = open(fname)
   m = []
   found_prop, found_identif = False, False
   while True:
      line = _readline(f)
      if line.startswith(prop):
         found_prop   = True
         found_identif  = False
      if identif in line        :
         found_identif  = True
      if line.startswith(stop_sign):
         found_prop   = False
         found_identif  = False
      if not (found_prop and found_identif): continue
      # if we get here, were gonna read data
      i = -1
      while True:
         i += 1
         line = _readline(f)
         if i < skip: continue
         if i > nrows + skip: break
         if line.startswith(stop_sign): break
         m.append([float(x) for x in line.split()])
      break
   f.close()
   return pl.array(m)

def write_property(m, prop, fname, fmt):
   '''
   m is matrix with shape (nx,ny,nz)
   fmt is typically '%i' or '%.2f' etc
   uses EQUALS (so not very effective)
   '''
   fmtstr = ' %s %s %%i %%i %%i %%i %%i %%i /\n' % (prop, fmt)
   file = open(fname, 'w')
   file.write("EQUALS\n")
   for k in range(m.shape[2]):
      for j in range(m.shape[1]):
         for i in range(m.shape[0]):
            file.write(fmtstr % (m[i,j,k],i+1,i+1,j+1,j+1,k+1,k+1))
   file.write("/\n")
   file.close()
   print fname, 'was written'

def write_property2(m, prop, fname, fmt, ncols=6):
   '''
   m is matrix with shape (nx,ny,nz)
   fmt is typically '%i' or '%.2f' etc
   does not use EQUALS
   '''
   file = open(fname, 'w')
   file.write("%s\n"%prop)
   colno = 0
   for k in range(m.shape[2]):
      for j in range(m.shape[1]):
         for i in range(m.shape[0]):
            file.write(' '+fmt%m[i,j,k])
            colno += 1
            if colno == ncols:
               file.write('\n')
               colno = 0
   file.write("/\n")
   file.close()
   print fname, 'was written'

def _read_relperm(swof_file):
   sfile = open(swof_file)
   header = True
   sw = []; kr_w = []; kr_o = []
   for line in sfile:
      if 'SWOF' in line:
         header = False
         continue
      if header: continue
      if line.startswith(DELIM) : break     # done
      if line.startswith('--'): continue  # skip comments
      rec = line.split()
      if len(rec) == 4:                   # allow blank lines
         sw.append(float(rec[0]))
         kr_w.append(float(rec[1]))
         kr_o.append(float(rec[2]))
   sfile.close()
   return (pl.array(sw), pl.array(kr_w), pl.array(kr_o))

def plot_relperms(swof_files):
   pl.figure()
   i = -1
   nfiles = len(swof_files)
   for swof_file in swof_files:
      i += 1
      sw, kr_w, kr_o = _read_relperm(swof_file)
      pl.plot(sw[:-1], kr_w[:-1], color=(i/float(nfiles), 0, 1-i/float(nfiles)), label=UT.basename(swof_file))
      pl.plot(sw[:-1], kr_o[:-1], color=(i/float(nfiles), 0, 1-i/float(nfiles)), label='_nolegend_')
   pl.xlabel('water saturation [-]')
   pl.ylabel('relperm [-]')
   pl.legend()
   pl.grid(True)
   pl.show()

def create_vfp_table1(q, dp_exponent, dp_max, ref_depth, table_number=1, fname=None):
   '''
   vfp-tables are for lift-curves.
   creates an oil-water vfp table where there is a linear wc-dependency.
   dp is given by c*q**dp_expon s.t. max(dp) = dp_max for q = max(q)
   q should be an array.
   '''
   c = dp_max/q[-1]**dp_exponent
   q_     = ' '.join(['%.1f'%x  for x in q])
   dp_wc1 = ' '.join(['%.1f'%x for x in c*q**dp_exponent])
   dp_wc0 = ' '.join(['%.1f'%x for x in pl.zeros(len(q))])
   txt = '''
VFPPROD
-- dp_exponent = %.1f dp_max = %.1f
 %i %i LIQ WCT GOR THP 3* /
 %s  /                   flow vals
 1 /                        THP
 0 1 /                      linear wct
 0 /                        gfr
 0 /                        alq
 1 1 1 1   %s   /  wc = 0
 1 2 1 1   %s   /  wc = 1
''' % (dp_exponent, dp_max, table_number, ref_depth, q_, dp_wc0, dp_wc1)
   if not fname is None:
      f = open(fname, 'w')
      f.write(txt)
      f.close()
   return txt

def create_vfp_table2(q, q_func, wc, wc_func, dp_max, ref_depth, shift=0.,table_number=1, fname=None):
   '''
   creates a oil-water vfp table where wc & flow rate dependency is user specified.
   dp is given by c*q_func(q)*wc_func(wc) s.t. max(dp) = dp_max for q = max(q) & wc = max(wc)
   q should be an array.
   the shift term is added.
   note: flow is pr. segment. so length of segments and valves pr. joint must be handled elsewhere.
   '''
   q_     = ' '.join(['%.1f'%x  for x in q])
   wc_    = ' '.join(['%.2f'%x  for x in wc])
   txt = '''
VFPPROD
-- q-function : %s
-- wc-function: %s
 %i %i LIQ WCT GOR THP 3* /
 %s / flow vals
 1. /                        THP
 %s / wct
 0. /                        gfr
 0. /                        alq
''' % (q_func.func_name, wc_func.func_name, table_number, ref_depth, q_, wc_)
   c = (dp_max-shift)/(q_func(q[-1])*wc_func(wc[-1]))
   dp = pl.zeros((len(wc), len(q)))
   for i in range(len(wc)):
      txt += '1 %i 1 1 ' % (i+1)
      for j in range(len(q)):
         dp[i,j] = c*q_func(q[j])*wc_func(wc[i]) + shift
         txt += '%.1f ' % dp[i,j]
      txt += ' /\n'
   if not fname is None:
      f = open(fname, 'w')
      f.write(txt)
      f.close()
   return dp, txt

def create_vfp_table3(q, q_func, wc, wc_func, gfr, gfr_func, dp_max, ref_depth, shift=0.,table_number=1, fname=None, thp=0.):
   '''
   creates a oil-water-gas vfp table where wc & flow rate dependency is user specified.
   dp is given by c*q_func(q)*wc_func(wc)*gfr_func(gfr) s.t. max(dp) = dp_max for q = max(q)
   q should be an array.
   the shift term is added.
   gfr: note! must be @ standard conditions
   note: flow is pr. segment. so length of segments and valves pr. joint must be handled elsewhere.
   '''
   q_   = format_long_line(' '.join(['%.1f'%x  for x in q]))
   wc_  = format_long_line(' '.join(['%.2f'%x  for x in wc]))
   gfr_ = format_long_line(' '.join(['%.2f'%x  for x in gfr]))
   txt = '''
VFPPROD
-- q-function : %s
-- wc-function: %s
 %i %i LIQ WCT GOR THP 2* BHP /
 %s / flow vals
 %.2f / THP
 %s / wct
 %s / gfr
 0. / alq
''' % (q_func.func_name, wc_func.func_name, table_number, ref_depth, q_, thp, wc_, gfr_)
   dp = pl.zeros((len(wc), len(gfr), len(q)))
   for i in range(len(wc)):
      for j in range(len(gfr)):
         txt += '1 %i %i 1 ' % ((i+1), (j+1))
         line = ''
         for k in range(len(q)):
            dp[i,j,k] = min(q_func(q[k])*wc_func(wc[i])*gfr_func(gfr[j]) + shift, dp_max)
            line += '%.1f ' % dp[i,j,k]
         txt += format_long_line(line) + ' /\n'
   if not fname is None:
      f = open(fname, 'w')
      f.write(txt)
      f.close()
   return dp, txt

def create_vfp_table4(s):
   '''
   this one is meant to be very flexible.
   s is a struct with following properties:
   #
   p_res : reservoir pressures (array). increasing
   rho_f : calculates density based on fluid-properties and phase saturations
   mu_f  : calculates viscosity based on fluid-properties and phase saturations
   dp_f  : calculates dp based on rho, mu, and flowrates
   Rs_f  : calculates dissolved gas as a function of pressure
   Bi_f  : calculates formation volume factor (i = o,w,g) as function of pressure
   eos_f : calculates gas density as a function of pressure
   q_l   : liquid rates (array)
   wct   : water cuts (array)
   gor   : surfacew gor (array)
   rho_i : density for fluid i (i = o,w,g)
   mu_i  : viscosity for fluid i (i = o,w,g)
   p_ref
   ref_depth
   table_number=1
   fname=None
   gor: note! must be @ standard conditions
   note: flow is pr. valve, and valves pr. segments (vps) must be given (so if length of segments varies, this must be handled elsewhere)
   here is some code to test it:
   #
dp_max         = 150.
valve_type     = 'ar3_peregrino'
def dp1(rho, mu, q):
   q *= 1000/24. # m3/d -> l/h
   return min(dp_max, VC.rcp_dp2(valve_type, rho, mu, q))
s              = UT.Struct()
s.p_res        = linspace(150,250, 11)                                    # reservoir pressures (array). increasing
s.rho_f        = lambda ro,rw,rg, so,sw,sg: ro*so + rw*sw + rg*sg         # calculates density based on fluid-properties and phase saturations
s.mu_f         = lambda mo,mw,mg, so,sw,sg: mo*so + mw*sw + mg*sg         # calculates viscosity based on fluid-properties and phase saturations
s.dp_f         = dp1                                                      # calculates dp based on rho, mu, and flowrates
s.Rs_f         = lambda p: 13.                                            # calculates dissolved gas as a function of pressure
s.Bo_f         = lambda p: 1.06678                                        # calculates oil formation volume factor as function of pressure
s.Bw_f         = lambda p: 1.012                                          # calculates water formation volume factor as function of pressure
s.Bg_f         = lambda p: 0.01                                           # calculates gas formation volume factor as function of pressure
s.q_l          = linspace(10, 2000, 15) * 24 / 1000.                      # typical liquid rates (array) for one valve [m3/d]
s.wct          = [0]                                                      # water cut (array)
s.gor          = linspace(13, 2000, 30)                                   # surface gor (array)
s.rho_o        = 800.                                                     # oil density
s.rho_w        = 1000.                                                    # water density
s.rho_g        = 1.2                                                      # gas density
s.mu_o         = 180.0                                                    # oil viscosity
s.mu_w         = 0.5                                                      # water viscosity
s.mu_g         = 0.015                                                    # gas viscosity
s.ref_depth    = 2000
s.table_number = 2
s.vps          = 2.                                                       # valves pr. segment
s.fname        = 't'                                                      # write VFPPROD to this file
s.append       = True                                                     # append to the given file
reload(ECL)
txt = ECL.create_vfp_table4(s)
   '''
   # convinent to create some strings
   q_   = format_long_line(' '.join(['%.1f'%x  for x in s.q_l * s.vps]))  # scaling flow
   wct_ = format_long_line(' '.join(['%.2f'%x  for x in s.wct]))
   gor_ = format_long_line(' '.join(['%.2f'%x  for x in s.gor]))
   thp_ = format_long_line(' '.join(['%.2f'%x  for x in s.p_res]))
   #
   txt = '''
VFPPROD
-- VPS = %f
-- dp-function : %s
 %i %i LIQ WCT GOR THP 2* BHP /
 %s / flow vals
 %s / THP
 %s / wct
 %s / gor
 0. / alq
''' % (s.vps, s.dp_f.func_name, s.table_number, s.ref_depth, q_, thp_, wct_, gor_)
   for i in range(len(s.p_res)):
      p = s.p_res[i]
      for j in range(len(s.wct)):
         for k in range(len(s.gor)):
            line = ' %i %i %i 1   ' % ((i+1), (j+1), (k+1))
            qw = s.q_l[0]*s.wct[j]
            qo = s.q_l[0] - qw
            qg = qo*s.gor[k] - qo*s.Rs_f(p) # want free gas
            # reservoir conditions
            qo *= s.Bo_f(p)
            qw *= s.Bw_f(p)
            qg *= s.Bg_f(p)
            rho_o = s.rho_o / s.Bo_f(p)
            rho_w = s.rho_w / s.Bw_f(p)
            rho_g = s.rho_g / s.Bg_f(p)
            # phase fractions
            qt = qo + qw + qg
            so = qo/qt
            sw = qw/qt
            if s.gor[k] < s.Rs_f(p): sg = 0.     # gor cannot be below Rs
            else                   : sg = qg/qt
            # calculate rho & mu (mixed)
            rho = s.rho_f(rho_o, rho_w, rho_g, so, sw, sg)
            mu = s.mu_f(s.mu_o(p), s.mu_w(p), s.mu_g(p), so, sw, sg)
            for q_l in s.q_l:
               q = qt*q_l/s.q_l[0]               # just need to scale it...
               dp = s.dp_f(rho, mu, q)
               line += '%.1e ' % dp
            line += ' /'
            txt += format_long_line(line) + '\n'
            #if k == 5: stop
            #print 'p=%i: %i %i %i %.1f' % (p,s.gor[k],sg*100,(q*1000/24.),dp)
   if not s.fname is None:
      mode = 'a' if s.append else 'w'
      f = open(s.fname, mode)
      f.write(txt)
      f.close()
   return txt

def create_vfp_table5(q, q_func, wc, wc_func, dpo_max, dp_cut, ref_depth, table_number=1, fname=None, thp=0., dp_xtra=0.):
   '''
   creates a oil-water vfp table where wc & flow rate dependency is user specified.
   dpo_max is the dp at wc=0, q = q_max
   dp_cut is just to avoid rediculously high values for high watercuts
   dp is given by c*q_func(q) where c = c(wc) s.t. dp = dpo_max for q = q - dq
   where dq is the difference between the oil flow rate and the flow rate for the given wc.
   q is an array.
   q_func is a scalar function, typically q**2
   wc is an array from 0 to 1
   wc_func gives values from 1 (for wc=1 => water) to some value > 1 (for wc=0 => oil)
   note1: flow is pr. segment. so length of segments and valves pr. joint must be handled elsewhere.
   note2: this is made to create liftcurves which are easy to make for the transition from pure
          water to pure oil. a linear wc_func gives a linear transition (as function of wc)
   note3: the default thp-value is 0, but it used to be hardcoded to 1. (so watch out if you have some old code)
   [q]  = m3/d
   [dp] = bar
   '''
   q_     = format_long_line(' '.join(['%.2f'%x  for x in q]))
   wc_    = format_long_line(' '.join(['%.2f'%x  for x in wc]))
   txt = '''
VFPPROD
-- q-function : %s
-- wc-function: %s
 %i %i LIQ WCT GOR THP 3* /
 %s / flow vals
 %.2f /                        THP
 %s / wct
 0. /                        gfr
 0. /                        alq
''' % (q_func.func_name, wc_func.func_name, table_number, ref_depth, q_, thp, wc_)
   dp = pl.zeros((len(wc), len(q)))
   # for clarity...
   dq_max = q[-1] - q[-1]/wc_func(0)  # difference between oil and water @ dpo_max
   y = wc_func(0) - wc_func(1)        # length of range of wc_func values
   for i in range(len(wc)):
      if y > 0: dq = dq_max * (wc_func(0)-wc_func(wc[i]))/y
      else:     dq = 0.
      c = dpo_max/q_func(q[-1]-dq)
      txt += ' 1 %i 1 1\n ' % (i+1)
      line = ''
      for j in range(len(q)):
         dp[i,j] = min(dp_cut, c*q_func(q[j])) + dp_xtra
         line += '%.2e ' % dp[i,j]
      txt += format_long_line(line) + ' /\n'
   if not fname is None:
      f = open(fname, 'w')
      f.write(txt)
      f.close()
   return dp, txt

def create_vfp_table6(q, q_func, gvf, gvf_func, dpo_max, dp_cut, ref_depth, table_number=1, fname=None, thp=0., dp_xtra=0.):
   '''
   creates a oil-gas vfp table where gvf & flow rate dependency is user specified.
   dp is given by c*q_func(q) where c = c(gvf) s.t. dp = dpo_max for q = q - dq
   where dq is the difference between the oil flow rate and the flow rate for the given gvf.
   gvf: gas volume fraction = Qg/(Qg+Qo+Qw). note: gvf should not be 1 (gas only)
   dpo_max is the dp at gvf=0, q = q_max
   dp_cut is just to avoid rediculously high values for high rates
   q is a vector.
   q_func is a scalar function, typically q**2
   gvf is an array from 0 to 1
   gvf_func gives values from 1 (for gvf=1 => gas) to some value > 1 (for gvf=0 => oil)
   note1: flow is pr. segment. so length of segments and valves pr. joint must be handled elsewhere.
   note2: this is made to create liftcurves which are easy to make for the transition from pure
          water to pure oil. a linear gvf_func gives a linear transition (as function of gvf)
   note3: the default thp-value is 0, but it used to be hardcoded to 1. (so watch out if you have some old code)
   [q]  = m3/d
   [dp] = bar
   '''
   gvf = pl.array(gvf)  # for convinience
   q_   = format_long_line(' '.join(['%.2f'%x  for x in q]))
   glr_ = format_long_line(' '.join(['%.2f'%x  for x in (gvf/(1.-gvf))])) # gvf to glr
   txt = '''
VFPPROD
-- q-function : %s
-- wc-function: %s
 %i %i LIQ WCT GLR THP 3* /
 %s / flow vals
 %.2f /                    thp
 0  /                        wct
 %s /                        gor
 0. /                        alq
''' % (q_func.func_name, gvf_func.func_name, table_number, ref_depth, q_, thp, glr_)
   dp = pl.zeros((len(gvf), len(q)))
   # for clarity...
   dq_max = abs(q[-1] - q[-1]/gvf_func(0))  # difference between oil and gas @ dpo_max
   y = gvf_func(0) - gvf_func(1)            # length of range of gvf_func values
   dbg = open('t', 'w')
   for i in range(len(gvf)):
      if y != 0: dq = dq_max * (gvf_func(0)-gvf_func(gvf[i]))/abs(y)
      else:      dq = 0.
      c = dpo_max/q_func(q[-1]-dq)
      txt += ' 1 1 %i 1\n ' % (i+1)
      dbg.write('gvf %.2f dq_max %.2f dq %.2f q %.2f y %.2f c %.4e\n' % (gvf[i], dq_max, dq, q[-1]-dq, y, c))
      line = ''
      sc = 1-gvf[i] if gvf[i] < 1 else 1.   # eclipse uses liquid rate - we want liq + gas
      for j in range(len(q)):
         dp[i,j] = min(dp_cut, c*q_func(q[j]/sc)) + dp_xtra
         line += '%.2e ' % dp[i,j]
      txt += format_long_line(line) + ' /\n'
   if not fname is None:
      f = open(fname, 'w')
      f.write(txt)
      f.close()
   dbg.close()
   return dp, txt

def icd_welsegs(nsegments, fname=None, diam=0.2, wellnm='wellnm'):
   '''
   just a lazy version of icd_helper.
   nsegments is number of segments for the well. (it is 11 if last entry
   of WELSEGS is segment number 12).
   put this section into WELSEGS.
   make sure COMPSEGS is updated so that the connected segments are the
   'icd-segments'. in practice, it means that the branch-number is increasing
   down the column.
   '''
   txt = '-- here comes the (A)ICD branches\n'
   for n in range(nsegments):
      segm = nsegments + 2 + n
      branch = n + 2
      txt += "  %i %i %i %i 0.1 0. %.2f 5E-05 /\n" % \
             (segm,segm, branch,branch, diam)
   txt += '''
WSEGAICD
  %s  %i %i (TO BE COMPLETED) /
/
''' % (wellnm, nsegments+2, 2*nsegments+1)
   if not fname is None:
      f = open(fname, 'w')
      f.write(txt)
      f.close()
      print 'created file', fname
   return txt

def read_trajectory(casenm, fname, is_ijk):
   from ert.ecl.ecl_grid import EclGrid
   g = EclGrid(casenm)
   f = open(fname)
   traj = []
   for line in f:
      rec = line.split()
      if not len(rec): continue
      p = (float(rec[0]), float(rec[1]), float(rec[2]))
      if is_ijk: x,y,z = g.get_xyz(ijk=(int(p[0])-1,int(p[1])-1,int(p[2])-1)) # the usual off-by-one
      else     : x,y,z = p
      traj.append([x,y,z])
   f.close()
   return traj

def _writeProd(traj, wname):
   well_definition = "%s.tmp.eclpost" % (wname)
   f = open(well_definition, 'w')
   print >>f, "----------------------------------------"
   print >>f, "WELL %s /" % (wname)
   print >>f, "----------------------------------------"    
   print >>f, "WELSPECS  COMPDAT WELSEGS"
   print >>f, "OIL"
   print >>f, "RWE 0.108 /"
   print >>f, "WELLXYZ"
   for p in traj:
      print >>f, "  %.1f   %.1f   %.1f" % (p[0],p[1],p[2])
   print >>f, "/"
   f.close()
   return well_definition

def create_well_schedule(casenm, traj_fname, wname, is_ijk=True, diam=0.1, roughness=1e-5, transf_by_mapaxes="Yes"):
   '''
   based on scrips createWell.py and extractTraj.py by Vegard Kippe.
   usage: create_well_schedule('0001_GI_OP_SHALLOW', 'OP-1.ijk', 'OP-1', True, 0.2, 5e-5)
   nb! must probably make sure no MAPAXES is active.
   the format of the traj_fname is simply
   i1 j1 k1
   i2 j2 k2
   ..
   iN jN kN
   '''
   #pdb.set_trace()
   traj = read_trajectory(casenm, traj_fname, is_ijk)
   well_definition = _writeProd(traj, wname)
   output_fn = '%s.SCH' % wname
   tmp_fn = "tmpXX.eclpost"
   f = open(tmp_fn, 'w')
   print >>f, "PROCESS INIT DATA"
   print >>f, "DEFINE WELL CONNECTION DATA"
   print >>f, "Global"
   print >>f, well_definition
   print >>f, casenm
   print >>f, output_fn
   print >>f, transf_by_mapaxes
   f.close()
   cmd = "eclpost %s" % (tmp_fn)
   os.system(cmd)
   cmd = 'sed -i s/DIAM/%f/ %s' % (diam, output_fn)
   os.system(cmd)
   cmd = 'sed -i s/ROUGH/%f/ %s' % (roughness, output_fn)
   os.system(cmd)
   print 'schedule section for well', wname, 'written to file', output_fn

def get_dimensions(datafile):
   f = open(datafile)
   rightplace = False
   for line in f:
      line = line.rstrip()
      if len(line) == 0       : continue
      if line.startswith('--'): continue
      if line.startswith('DIMENS'):
         rightplace = True
         continue
      if not rightplace: continue
      nx, ny, nz = line.split()[:3]
      break
   f.close()
   return (int(nx), int(ny), int(nz))

def get_icd_segments(datafile, wellnm):
   '''
   useful to get icd-segment numbering from data file.
   assumes that segments are sequential.
   assumes that section is ended with a DELIM on a separate line.
   assumes only one WSEGAICD or WSEGVALV section in file.
   '''
   smin = pl.Inf; smax = -pl.Inf
   f = open(datafile)
   rightplace = False
   for line in f:
      line = line.rstrip()
      if len(line) == 0       : continue
      if line.startswith('--'): continue
      if line.startswith('WSEGAICD') or line.startswith('WSEGVALV'):
         rightplace = True
         continue
      if not rightplace: continue
      if line.startswith(DELIM): break
      if not wellnm in line  : continue
      s1, s2 = line.split()[1:3]
      smin = min(int(s1),smin)
      smax = max(int(s2),smax)
   f.close()
   return pl.arange(smin, smax+1)

def read_variable(casenm, varnm, report_only=False, force_creation=True, store_at='.', skiprows=2, delete_afterwards=False):
   '''
   reads a variable from an Eclipse summary file using
   symmary.x (Unix-tool). creates a file called <store_at/casenm>.<varnm>
   if force_creation is False, it will create file only if it does not exist already.
   returns time and variable-vector as a matrix
   use store_at if you wanna collect data-files somewhere else (more clean)
   '''
   varfile = '%s/%s.%s' % (store_at, casenm, varnm)
   if force_creation or not os.path.exists(varfile):
      opt = '--report-only' if report_only else ''
      cmd = 'summary.x %s %s %s > %s' % (opt, casenm, varnm, varfile)
      os.system(cmd)
   d = pl.loadtxt(varfile, skiprows=skiprows, usecols=(0,2))
   if delete_afterwards: os.unlink(varfile)
   return d

def format_long_line(line, max_length=80, prepend=' '):
   items = line.split()
   lines = []
   newline = ''
   for item in items:
      if len(newline) >= max_length - len(prepend):
         lines.append(newline)
         newline = prepend
      newline += (item + ' ')
   lines.append(newline)      # dont forget this one
   return '\n'.join(lines)

def read_keyword_section(fname, keyword):
   '''
   reads eclipse data-file (or include-file) searching for given keyword.
   returns the data connected to this keyword.
   note1: will get in trouble if DELIM is part of the data (like a filename)
   note2: consider using DataDeck in stead
   '''
   f = open(fname)
   lines = []
   right_place = False
   while True:
      line = f.readline()
      if not line: break
      if right_place or line.startswith(keyword):
         if not right_place:
            right_place = True
            continue
         line = line.strip()
         while line.startswith('--') or len(line) == 0: # skip comments and empty lines
            line = f.readline()
            line = line.strip()
         while not DELIM in line:     # append lines until DELIM is found
            line2 = f.readline()
            line2 = line2.strip()
            while line2.startswith('--') or len(line2) == 0:
               line2 = f.readline()
            line += ' ' + line2
         line = line.strip()
         if line == DELIM: break
         line, xtra = line.split(DELIM,1)  # remove comments after DELIM
         lines.append(line)
         if xtra.strip().startswith(DELIM): break
   return lines


def wsegvalv(l_segms, r, vpj=1, Cv=0.6):
   '''
   creates the WSEGVALV keyword for a nozzle.
   returns a formatted string and a list of the areas used.
   l_segms : length of each segment along the well
   r       : nozzle radius [m]
   vpj     : number of valves per joint
   '''
   A0 = pl.pi*r**2
   s = 'WSEGVALV\n'
   segno = len(l_segms) + 3 
   As = []   # for convinience
   for l in l_segms:
      A = vpj*l/JOINT_LENGTH*A0
      s += '  %i  %.2f  %.3e /\n' % (segno, Cv, A)
      As.append(A)
      segno += 1
   return s + '/\n', As

def convert_egrid(casenm):
   inpfile = '/tmp/inp.txt'
   print inpfile, 'created'
   f = open(inpfile, 'w')
   s = '''
U
%s
1
8
EGRID
N
''' % casenm
   f.write(s)
   f.close()
   cmd = '@convert < %s' % inpfile
   print cmd
   os.system(cmd)

def utm_coordinates(fegrid):
   '''
   uses eclpost to write COORD in utm coordinates.
   assumes that the fegrid file has mapaxes included.
   eclpost creates a fegrid-like file with extension FILLED
   NOTE! Deprecated. for some reason eclpost scales the grid
         by a factor 0.2, so it becomes useless.
   '''
   bname = UT.basename(fegrid)
   filled = bname+'.FILLED'
   if os.path.exists(filled): return filled
   s = '''
Process Grid Data
Grid data extraction
Coordinate data
Z
Eclipse input format
Z
A
%s
%s
Yes
''' % (bname, bname)
   f = open(UT.basename(UT.tmpfile()), 'w')
   f.write(s)
   f.close()
   os.system('eclpost %s'%f.name)
   os.unlink(f.name)
   return filled



class SimpleSectorModel(object):
   def __init__(self, casenm, tstep=0,
                grkws=['PORO', 'PERMZ', 'PERMY', 'PERMX', 'NTG', 'ActiveCell'],
                edkws=['TRANX', 'TRANY', 'TRANZ'],
                regkws=['SATNUM'],
                solkws=['SWAT', 'SGAS', 'PRESSURE', 'RS', 'RV']):
      '''
      This is a class that creates a sector model cut from a given full field model.
      The advantage is that if the full field model is itself a sector model or gets the initial state
      from some other simulation (using RESTART etc.), it is problematic/impossible to create a 'standard'
      sector model. In addition, the files produced by this sector model are small.
      What it does, is that it gets all grid and reservoir properties from EGRID and UNRTST files, and
      the dynamic properties (pressure and saturations) from UNRST.
      It exports the needed properties using FloViz.
      You will get a stand-alone model that is put in a separate directory (see create()).
      #
      Note1: Will not work for LGRs
      #
      casenm : full field model to create sector from
      tstep  : time step to read initial solution from. default is 0 - which is the initial solution
               (simulate full field model for just 1 day if you dont have the initial solution you want)
      *kws   : if not-standard properties are to be used. for example, solkws should be set to
               ['SWAT', 'PRESSURE', ...] for two-phase oil/water cases (SGAS is not available)
      '''
      self.casenm   = casenm
      self.cmdfile  = None
      self.files    = None
      self.tstep    = tstep
      self.propsdir = 'Props_%s_%i' % (casenm, tstep)
      self.grkws    = grkws
      self.edkws    = edkws
      self.regkws   = regkws
      self.solkws   = solkws
      if not os.path.exists(self.propsdir): os.mkdir(self.propsdir)
      # dimensions of full field...
      dd = DataDeck(casenm+'.DATA')
      self.nx,self.ny,self.nz = [int(x) for x in dd.get_raw_data('DIMENS')[0]]
      #
      # run floviz to export properties
      self._create_export_file()
      logfile = 'FloViz.log'
      if os.path.exists(logfile): os.unlink(logfile)
      cmd = ('@floviz', '-local', '-version', '2009.2', '-play', self.cmdfile)
      p = subprocess.Popen(cmd)
      #
      while(True):
         if os.path.exists(logfile): break                                         # floviz has started...
         time.sleep(1)
      while(True):
         if UT.grep_column(logfile, 'Program finished ok',1,is_float=False): break # floviz is done
         print 'waiting for FloViz to finsh...'
         time.sleep(1)
      # kill floviz. p.terminate() does not work since @floviz forks another process
      ps = [p for p in psutil.process_iter() if len(p.cmdline())>1 and 'floviz' in p.cmdline()[0]]
      for proc in ps:
         if len(proc.cmdline()) == 4 and proc.cmdline()[3] == self.cmdfile:
            print 'killing floviz'
            os.kill(proc.pid, signal.SIGKILL)
            break
#
   def _create_grid(self, box, mapaxes):
      box = pl.array(box)   # for convinience
      # make sure we have a FEGRID-file
      fegrid = '%s.FEGRID'%self.casenm
      if not os.path.exists(fegrid): convert_egrid(self.casenm)
      gridfile = '%s/GRID.grdecl' % self._incname(box)
      c,z  = read_fegrid(fegrid, self.nx,self.ny,self.nz)
      write_decimated_grid(c,z, box, gridfile, mapaxes=mapaxes, ncols=4)
      return gridfile
#
   def create(self, sectornm, box, datafile=None, cleanup=True, mapaxes=None):
      '''
      creates a sector model based on box. it will create a stand-alone model put into a
      directory called sectornm. here, all include files are copied, so that they can be
      changed if necesarry (we are using pack_sim for this).
      the sector datafile will be based on the full field model (or by the provided datafile),
      where a lot of the properties are now just read (included) from the full field model.
      in the SCHEDULE section, all wells are shifted according to the box indices.
      mapaxes: the six values needed for the MAPAXES keyword of the full field model (or None)
      Note: to get the UTM coordinates right, mapaxes must be provided here, or it must be put
            manually in the GRID-section of the data file.
      '''
      if not datafile: datafile = '%s.DATA' % self.casenm
      # create new include files for our box
      sec = {}
      #
      # GRID section
      if self.files.gr:
         s = self._create_includefiles(self.files.gr, box)
         gridfile = self._create_grid(box, mapaxes)
         sec['GRID'] = s + "INCLUDE\n '%s' /\n\n"%gridfile
         sec['GRID'] += 'GRIDFILE\n 0 1 /\n\nINIT\nRPTGRID\n /\n'    # always useful
      #
      # EDIT section
      if self.files.ed:
         s = self._create_includefiles(self.files.ed, box)
         sec['EDIT'] = s
      #
      # REGIONS section
      if self.files.reg:
         s = self._create_includefiles(self.files.reg, box)
         sec['REGIONS'] = s
      #
      # SOLUTION section
      if self.files.sol:
         s = self._create_includefiles(self.files.sol, box)
         sec['SOLUTION'] = s
      #
      # SOLUTION section
      sec['SCHEDULE'] = shift_schedule_section(datafile, box)
      #
      # put it together
      datafile2 = '%s_SM.DATA'%UT.basename(datafile)
      rewrite_datafile(sec, box, datafile, datafile2)
      #
      # finally, use pack_sim to put it in a separate directory - and cleanup
      if os.path.exists(sectornm): os.system('rm -rf %s' % sectornm)
      os.system('pack_sim %s %s' % (datafile2, sectornm))
      os.rename('%s/%s'%(sectornm,datafile2), '%s/%s.DATA'%(sectornm,UT.basename(datafile)))
      #
      # for QC
      print sectornm, 'has been created'
      cmd = 'find %s' % sectornm
      os.system(cmd)
      if cleanup:
         cmd = 'rm -rf %s %s' % (self._incname(box), datafile2)
         os.system(cmd)
#
   def _create_export_file(self):
      s1 = '''
      ImportGrid( Vendor = "Eclipse",   &
         Location = &
         "%s.EGRID" &
         ,  Grid = "Grid1",  PropManager = "PropManager1" )
      ReadProperties( Vendor = "Eclipse",   &
         Location = &
         "%s.INIT, &
   InitialProperty",  PropManager = "PropManager1",  DefaultToPropUnits = "TRUE",   &
         ReadNow = "TRUE" )
      ReadProperties( Vendor = "Eclipse",   &
         Location = &
         "%s.UNRST, &
   MinReportStep, %i, MaxReportSTep, %i, RestartProperty",   &
         PropManager = "PropManager1",  DefaultToPropUnits = "TRUE",   &
         ReadNow = "FALSE" )
      ImportWells( Vendor = "Eclipse",   &
         Location = &
         "%s.UNRST, &
   MinReportStep, %i, MaxReportSTep, %i",  Grid = "Grid1",  Wells = "Wells1",   &
         ReadNow = "FALSE" )
      ReadRestartFiles(  )
      UpdateGridUnits( Grid = "Grid1" )
      UpdateWellUnits( Wells = "Wells1" )
      ReportGridDimensions( Grid = "Grid1" )
      AddSimulationToView( Viewer = "3DViewer1",  Grid = "Grid1",   &
         GridComposite = "GridComposite1",  Wells = "Wells1",   &
         WellComposite = "WellComposite1" )
      ''' % (self.casenm, self.casenm, self.casenm, self.tstep, self.tstep, self.casenm, self.tstep, self.tstep)
      cmdfile = 'export_%s.CMDLOG' % (self.casenm)
      f = open(cmdfile, 'w')
      f.write(s1)
      templ = '''
      ExportProps( Vendor = "Eclipse",   &
         Location = &
         "%s" &
         ,  Grid = "Grid1",  PropManager = "PropManager1",  PropNames = " %s &
            ",  Components = "COARSE_AND_LGRS" %s)'''
      files = UT.Struct()
      files.gr, files.ed, files.reg, files.sol = (None, None, None, None)
      for kws in [self.grkws, self.edkws, self.regkws, self.solkws]:
         if not kws: continue
         if kws == self.solkws:
            fname = '%s/%s_%i_solution.txt' % (self.propsdir, self.casenm, self.tstep)
            rstep = ' ,  ReportStep = %i ' % self.tstep
            files.sol = fname
         else:
            rstep = ''
         if kws == self.grkws:
            fname = '%s/%s_%i_grid.txt' % (self.propsdir, self.casenm, self.tstep)
            files.gr = fname
         if kws == self.edkws:
            fname = '%s/%s_%i_edit.txt' % (self.propsdir, self.casenm, self.tstep)
            files.ed = fname
         if kws == self.regkws:
            fname = '%s/%s_%i_regions.txt' % (self.propsdir, self.casenm, self.tstep)
            files.reg = fname
         kwstr = ''.join(["&\n          '%s'    " % x for x in kws])
         f.write(templ % (fname, kwstr, rstep))
         print fname, 'will be created'
      f.close()
      print cmdfile, 'has been created'
      self.cmdfile = cmdfile
      self.files   = files
#
   def _incname(self, box):
      incname = 'Include_%s_%i__%i_%i_%i_%i_%i_%i' % (self.casenm, self.tstep, box[0],box[1], box[2],box[3], box[4],box[5])
      if not os.path.exists(incname): os.mkdir(incname)
      return incname
#
   def _create_includefiles(self, propfile, box):
      '''
      creates a sector model based on exported properties and solution variables from floviz.
      '''
      #
      # PART 1: read property-file into data matrices
      kws = []
      f = open(propfile)
      data = []
      active = False
      while True:
         line = f.readline()
         if not line: break
         if line.startswith('BOX'):
            # read dimensions and initialize
            line = f.readline()
            rec = line.split()
            i1 = int(rec[0]); i2 = int(rec[1])
            j1 = int(rec[2]); j2 = int(rec[3])
            k1 = int(rec[4]); k2 = int(rec[5])
            m = pl.zeros((i2-i1+1,j2-j1+1,k2-k1+1))
            f.readline()
            kw = _readline(f)
            kws.append(kw)
         line = f.readline()
         rec = line.split()
         ncols = len(rec)
         n = -1
         for k in pl.arange(k2-k1+1):
            for j in pl.arange(j2-j1+1):
               for i in pl.arange(i2-i1+1):
                  n += 1
                  if n == ncols:
                     line = f.readline()
                     rec = line.split()
                     ncols = len(rec)
                     n = 0
                  m[i,j,k] = rec[n]
                  if kw in ('SGAS' ,'SWAT') and m[i,j,k] < 0: m[i,j,k] = 0.0
         f.readline()
         f.readline()
         data.append(m)
      f.close()
      #
      # PART 2: write decimated properties
      i1,i2, j1,j2, k1,k2 = box
      includestr = ''
      for i in range(len(kws)):
         kw = kws[i]
         fmt = '%i' if kw in ('SATNUM', 'ACTCELL', 'PVTNUM', 'DOMAINS') else '%.3e'
         if kw == 'ACTCELL': kw = 'ACTNUM'
         fname = '%s/%s.grdecl' % (self._incname(box), kw)
         write_property2(data[i][i1-1:i2,j1-1:j2,k1-1:k2], kw, fname, fmt)
         includestr += "INCLUDE\n '%s' /\n\n" % fname
      self._kws = kws    # for debug
      return includestr
#
   def cleanup(self):
      '''
      after this, this instance cannot be used...
      '''
      [os.unlink(x) for x in UT.glob(['%s.FEGRID'%self.casenm, '%s.FILLED'%self.casenm, 'FloViz*', self.cmdfile])]
      os.system('rm -rf %s' % self.propsdir)

def read_fegrid(fname, nx, ny, nz):
   '''
   reads FEGRID file into two matrices:
    -coord with shape = [nx+1,ny+1,6]
     coord[0,0,:] is x1,y1,z1 + x2,y2,z2 of two point p1 & p2 defining a 'vertical' grid-line.
    -zcorn with shape [2nx, 2ny, 2nz]
   '''
   f = open(fname)
   ni = nx+1
   nj = ny+1
   #
   # reading COORD
   foundit = False
   coord = pl.zeros((ni,nj,6))
   n = 0
   i,j,k = (0,0,0)
   while n < ni*nj*6:
      line = f.readline().strip()
      if 'COORD' in line:
         foundit = True
         line = f.readline().strip()
      if not foundit: continue
      for v in line.split():
         coord[i,j,k] = float(v)
         if k == 5:
            k = -1
            i += 1
         if i == ni:
            i = 0
            j += 1
         k += 1
         n += 1
   #
   # reading ZCORN
   foundit = False
   zcorn = pl.zeros((2*nx,2*ny,2*nz))
   n = 0
   i,j,k,= (-1,0,0)
   #while n < 8*(nx-1)*(ny-1)*(nz-1):
   while n < 8*nx*ny*nz:
      line = f.readline().strip()
      if 'ZCORN' in line:
         foundit = True
         line = f.readline().strip()
      if not foundit: continue
      for v in line.split():
         i += 1
         n += 1
         if i == 2*nx:
            i = 0
            j += 1
         if j == 2*ny:
            j = 0
            k += 1
         zcorn[i,j,k] = float(v)
   f.close()
   return coord, zcorn

def write_decimated_grid(coord, zcorn, box, fname, mapaxes=None, ncols=6):
   '''
   writes a sector (given by box) of the provided grid.
   usually, coord and zcorn is read from an FEGRID file using read_fegrid.
   useful when creating sector models
   mapaxes: the six values needed for the MAPAXES keyword of the full field model (or None)
   '''
   i1,i2, j1,j2, k1,k2 = pl.array(box)-1
   f = open(fname, 'w')
   if mapaxes is not None:
      f.write('MAPAXES\n %.3f %.3f\n %.3f %.3f\n %.3f %.3f \n/\n\n' %
            (mapaxes[0],mapaxes[1],mapaxes[2],mapaxes[3],mapaxes[4],mapaxes[5]))
   f.write("COORD\n")
   colno = 0
   for j in range(j1, j2+2):
      for i in range(i1, i2+2):
         for n in range(6):
            f.write(' %.8e'%coord[i,j,n])
            colno += 1
            if colno == ncols:
               f.write('\n')
               colno = 0
   f.write("\n/\n\n")
   f.write("ZCORN\n")
   i1,i2, j1,j2, k1,k2 = 2*(pl.array(box)-1)
   colno = 0
   for k in range(k1, k2+2):
      for j in range(j1, j2+2):
         for i in range(i1, i2+2):
            f.write(' %.8e'%zcorn[i,j,k])
            colno += 1
            if colno == ncols:
               f.write('\n')
               colno = 0
   f.write("\n/\n")
   f.close()
   print fname, 'was written'


def shift_schedule_section(datafile, box, fname=None, sched_only=True, remove_after_delim=True):
   '''
   Note1: schedule part must be in the datafile - not in an include
   Note2: will not work for LGRs
   '''
   dd = DataDeck(datafile, remove_after_delim=remove_after_delim)
   f = open(datafile)
   s1 = ''
   #
   # before SCHEDULE => keep all
   while True:
      line = f.readline()
      s1 += line
      kw = line.split()[0].rstrip() if line[0].isupper() else ' '
      if kw == 'SCHEDULE': break
      else               : continue
   #
   # in SCHEDULE => keep everything except WELSPECS, COMPDAT, and COMPSEGS
   all_kws = [x.split('-')[0] for x in dd._chunks.keys()]
   well_kws = ('WELSPECS', 'COMPDAT', 'COMPSEGS')
   s3 = ''
   while True:
      line = f.readline()
      if not line: break
      kw = line.split()[0].rstrip() if line[0].isupper() else ' '
      if kw in well_kws:  # skip these...
         while True:
            line = f.readline()
            if not line: break
            kw = line.split()[0].rstrip() if line[0].isupper() else ' '
            if not kw or kw in well_kws: continue
            if kw in all_kws:
               s3 += line
               break
      else: s3 += line
   f.close()
   #
   # define well - shifting WELSPECS, COMPDAT, and COMPSEGS
   #
   #  -- WELSPECS
   data = dd.get_raw_data('WELSPECS')
   s2 = '\nWELSPECS\n'
   wells = []
   for x in data:
      wells.append(x[0])
      s2 += ' %s %s' % (x[0], x[1])
      for i in range(2,5):
         if x[i] == '1*': s2 += ' %s' % x[i]
         else           : s2 += ' %i' % (int(x[i])-box[2*(i-2)]+1)
      s2 += ' ' + ' '.join(x[5:]) + ' /\n'
   s2 += '/\n'
   #
   #  -- COMPDAT
   data = dd.get_raw_data('COMPDAT')
   s2 += '\nCOMPDAT\n'
   for x in data:
      s2 += ' %s' % x[0]
      for i in range(1,4):
         if x[i] == '1*': s2 += ' %s' % x[i]
         else           : s2 += ' %i' % (int(x[i])-box[2*(i-1)]+1)
      if x[4] == '1*': s2 += ' %s' % x[4]
      else           : s2 += ' %i' % (int(x[4])-box[4]+1)
      s2 += ' ' + ' '.join(x[5:]) + ' /\n'
   s2 += '/\n'
   #
   #  -- COMPSEGS
   for well in wells:
      data = dd.get_raw_data('COMPSEGS', well)
      s2 += '\nCOMPSEGS\n %s /\n' % data[0][0]
      for x in data[1:]:
         for i in range(0,3):
            if x[i] == '1*': s2 += ' %s' % x[i]
            else           : s2 += ' %i' % (int(x[i])-box[2*i]+1)
         s2 += ' ' + ' '.join(x[3:]) + ' /\n'
      s2 += '/\n'
   #
   # compile
   if sched_only: s = s2 + s3
   else:          s = s1 + s2 + s3
   if fname:
      f = open(fname, 'w')
      f.write(s)
      f.close()
      print fname, 'created'
   #stop
   return s


def rewrite_datafile(sections, box, datafile, datafile2=None):
      '''
      replaces sections of the data-file with the given content. the other
      sections are left untouched.
      sections: {'SUMMARY' : "INCLUDE\n 'include/sum.inc' /\n"}
      '''
      f = open(datafile)
      skiprest = False
      s = '''
-- Based on %s
-- Box = [%i %i  %i %i  %i %i]\n\n''' % (datafile, box[0], box[1], box[2], box[3], box[4], box[5])
      while True:
         line = f.readline()
         if not line: break
         kw = line.split()[0].rstrip() if line[0].isupper() else ' '
         if box and kw == 'DIMENS':        # replace DIMENS
            i1,i2, j1,j2, k1,k2 = box
            s += '\nDIMENS\n %i  %i  %i /\n' % (i2-i1+1, j2-j1+1, k2-k1+1)
            f.readline()
            continue
         if kw in SECTIONS:
            skiprest = False
            if kw in sections.keys():
               s += '%s\n\n%s' % (kw, sections[kw])
               skiprest = True
               continue
         if skiprest: continue
         s += line
      f.close()
      if datafile2:
         print 'creating', datafile2
         f = open(datafile2, 'w')
         f.write(s)
         f.close()
         return
      return s

def _create_qc_datafile(repl, fname, ntsteps):
   templstr = '''
NOECHO
-- ====================================================================
RUNSPEC
-- ====================================================================
 
TITLE
 Model for QC of inflow characterstics
 
--Model dimensions 
DIMENS
-- NX  NY  NZ
   2   3   2 /
 
-- Phases present
OIL
GAS
WATER
 
-- Units
METRIC
 
-- unified and ascii-formatted files in&out
UNIFIN 
UNIFOUT 
--FMTOUT

-- Run start date
START
--  DAY   MONTH  YEAR
    31 'Aug' 2002   /

-- Linear solver stack size (had problems with convergence)
NSTACK
 50 /
 
-----------------------------------------------------
-- Dimensions
-----------------------------------------------------
 
TABDIMS
-- NTS  NTPVT  NSS  NPPVT  NRPVT
    1    1     23    24     11  /
 
REGDIMS
-- NTFIP  NMFIPR  NRFREG  NTFREG
     1     1       0       0  /
 
-- Well dimensions
WELLDIMS
-- NWMAXZ  NCWMAX  NGMAXZ  NWGMAX
     4      156       2      3  /

WSEGDIMS
-- NSWLMX   NSEGMX   NLBRMX
    1        240      135   /

-----------------------------------
-- Run control settings
-----------------------------------
--NOSIM
NOINSPEC        -- Prevent output of unnecessary files
NORSSPEC

VFPPDIMS
 _VFPPDIMS_ /
 
-- ========================================================================
--                                            GRID SECTION
-- ========================================================================
GRID 

DYV
 3*12 /
DXV
 2*12 /
DZ
 12*1 /

-- set params for the complete domain
EQUALS
  PORO   0.3 /      -- tried lower poro on well cell (less transient effects) but no good
  PERMX  1000 /
  TOPS   1499 1 2 1 3 1 1  /
/

COPY
 'PERMX'  'PERMY'   /
 'PERMX'  'PERMZ'   /
/

INIT
 
GRIDFILE 
  2 / 

RPTGRID
 22*0 /

ACTNUM
 0 0   1 0   0 0
 1 0   1 1   1 0
/

-- ========================================================================
--                                            EDIT SECTION
-- ========================================================================
EDIT
-- =======================================================================
--                                             PROPS SECTION
-- =======================================================================
PROPS

 
-- Pbub      Bo@Pbub    Muo@Pbub
PVDO 
     1      1.0001      _MUO_1_
 10000      1.0000      _MUO_2_
/

--  Pg      Bg     Mug
PVDG 
     1      1.0001     _MUG_1_
 10000      1.0000     _MUG_2_
/
  
-- Pref       Bw  Cw         Muw      Cv
PVTW
   1          1.0  4.30E-5  _MUW_     0.0 / 


 
--    Do      Dw       Dg    (Kg/m3)
DENSITY
     _RHO_O_ _RHO_W_ _RHO_G_  /


-- Saturation tables - Linear rel.perm tables

SWFN
--   Sw         krw        Pcwo
   0.0          0.0        0.0
   1.0          1.0        0.0
/

SGFN
--   Sg         krg        Pcgo
   0.0          0.0        0.0
   1.0          1.0        0.0
/

SOF3
--   So        krow        krog
   0.0          0.0        0.0
   1.0          1.0        1.0
/

----------------------------------------------------------------------------
-- Formation compressibility data
----------------------------------------------------------------------------
-- Pref     Cr
ROCK
   158.2    2.85E-4 /
 
-- ======================================================================
--                                           REGIONS SECTION
-- ======================================================================
REGIONS

SOLUTION

MESSAGES 
 12*99999999 /
 
--EQUIL
-- DATUM      P@DATUM   OWC    PC@OWC  GOC   PC@GOC   Rs   Rv   N
--  1544.50       150.0    9999.0   0.0   9    0.0    0    0   20 /
PRESSURE
 12*100 /

SWAT
 0 0 0 0 0 0
 0 0 0 1 0 0
/

SGAS
 0 0 0 0 0 0
 0 0 0 0 1 0
/

-- INITIAL RESTART FILE
RPTSOL
  RESTART=2            -- Create initial restart file
  FIP=2                -- Volumes for all FIPNUMs
/

-- ======================================================================
--                                            SUMMARY SECTION
-- ======================================================================
SUMMARY

ELAPSED
TELAPDAY
DATE

RUNSUM
EXCEL

FPR
FGOR
FOPR
FGPR
FGPRF
FWPR
FOPT
FGPT
FOE
FOIP
FGIP
FVPR
FVPT

WBHP
/
WBP
/
WTHP
/
WPI
/
WMCTL
/
WOPR
/
WOIR
/
WWIR
/
WGIR
/

TCPU
TCPUTS
TIMESTEP
TCPUDAY
NEWTON
MSUMNEWT
NLINEARS
MSUMLINS

SPRD
/
SOFR
/
SWFR
/
SGFR
/
SGFRF
/
SWCT
/
SGOR
/

-- ======================================================================
--                                           SCHEDULE SECTION
-- ======================================================================
SCHEDULE
RPTRST         Controls on output to the RESTART file
 BASIC=3
 FREQ=%i
/ 

RPTSCHED
 'NEWTON=1'      -- linear eq summary for each time step to prt-file (=1)
 'SUMMARY=1'     -- solution summary for each time step to prt-file (=1)
 'WELSPECS'
/

WELSPECS
 P1    PROD   1 2 1* OIL   2* SHUT YES 1* SEG 0 /
 I_OIL INJ    1 1 1* OIL   2* SHUT YES /
 I_GAS INJ    1 3 1* GAS   2* SHUT YES /
 I_WAT INJ    2 2 1* WATER 2* SHUT YES /
/

COMPDAT   Well completion specification data
--                                             inflow-
--Name	 I  J  K  K   status satnum   CF	     diam     Kh  Skin Dfac Direction	pressureRadius	
 P1       1* 1* 1 1  OPEN     1*      1*     0.114  1*      2*      Y              1*    /
 I_OIL    1* 1* 2 2  OPEN     1*      1*     0.114  1*      2*      Z              1*    /
 I_GAS    1* 1* 2 2  OPEN     1*      1*     0.114  1*      2*      Z              1*    /
 I_WAT    1* 1* 2 2  OPEN     1*      1*     0.114  1*      2*      Z              1*    /
/

WELSEGS     Defines the segment structure of a multisegment well
-- name  BHP      BHP     Volume  Len/Dpt
--       ref-dpt  ref-MD          type
P1      1500.0   1500.0	   1*   	 INC	 /				
--                                       Hydraulic-
-- First Last Branch Outlet Length Depth Diam     Rough
-- segm. segm num    segm          chng     
      2    2    1     1      12    0.00  0.114 1.50E-05   /
-- here comes the (a)icd-branches
   3    3    2      2        0.1   0.0   0.1    1.50E-05   /
/

COMPSEGS    Defines location of completions in a multisegment well
P1 /
--I  J  K branch start end dir ind depth
 1  2  1  2    1500.0    1*   Y   1  1500.0 /
/

-- Basic tuning 
------------------------------------------------------------------------
-- Initial time step
--         Max time step
TUNING
-- Record 1: Time stepping controls
-- TSINIT  TSMAXZ  TSMINZ  TSMCHP  TSFMAX  TSFMIN  TSFCNV  TFDIFF
/
-- Record 2: Time truncation and convergence controls
-- TRGTTE  TRGCNV  TRGLCV  XXXTTE  XXXCNV  XXXMBE  XXXLCV  XXXWFL  TRGFIP   
/
-- Record 3: Control of Newton and linear iterations
-- NEWTMX  NEWTMN  LITMAX  LITMIN  MXWSIT  MXWPIT  DDPLIM  DDSLIM  TRGDPR  XXXDPR  
      1*     1*      50
/

-- Iteration parameters for multi-segment wells
WSEGITER
 40 1* 0.2 /

------------------------------------------------------------------------
-- running the simulation
------------------------------------------------------------------------

-- Name     Open?   Control       Orat     Wrat     Grat  Lrat RFV  BHP
WCONPROD
 P1      OPEN      BHP            1*       1*       1*    1*   1*   100 /
/

-- create RFT file which contains well profiles (per timestep)
--WRFTPLT
--Well     RFT   PLT   Segment
--Name     Data  Data  Data
-- P1       REPT  REPT  REPT   /   
--/

-- debug info ala Fredric Spiegel
--WELDEBUG
-- P1 3000 /  values > 0 gives debug output (pr. position)
--/

_SCHEDULE_

''' % ntsteps
   newstr = templstr
   for r in repl:
      newstr = newstr.replace(r[0], r[1])
   f = open(fname, 'w')
   f.write(newstr)
   f.close()
   print f.name, 'is created'

def create_qc_datafile(rho, mu, qs, gvfs, wcts, vfppdims, vfpprod, wsegaicd=None, fname='TT.DATA', ntsteps=12, dt=7.):
   '''
   creates a datafile for veryfying inflow characteristics, especially when using VFP-tables.
   the model is a small 'pad' with 5 active cells, 3 for injecting oil, gas, and water, respectively, and one is producing.
   it will run the model len(qs)*len(gvfs)*len(wcts) times - each TSTEP ntsteps*dt times.
   either wsegaicd or vfpprod must be given.
   wsegaicd-coefficient are given in the order required by the WSEGAICD keyword: a_aicd, rho_cal, mu_cal, x, y.
   vfpprod is a string giving the full Eclipse VFPPROD input for *1* valve.
   vfpprod could be made by ECL.create_vfp_table5(linspace(1,79,40), lambda x: x**2, (0,1), lambda x: 2-x, 40, 150, 1500, 1)[1].
   vfppdims is (#flowrates, #wcts, # gvfs)
   rho = (rho_o, rho_g, rho_w).
   mu  = (mu_o, mu_g, mu_w).
   units: m3/d, bar.
   note: does not handle gvf = 1 (pure gas)
   '''
   rho_o, rho_g, rho_w = rho
   mu_o, mu_g, mu_w = mu
   if type(gvfs) in (float, int): gvfs = [gvfs]
   if type(wcts) in (float, int): wcts = [wcts]
   k,j = -1,-1
   if wsegaicd:
      schedule = 'WSEGAICD\n P1   3   3   %.3e   12  %i %.2f 0.7 1* 1*  1* 1* %.2f %.2f  /\n/' % tuple(wsegaicd)
   else:
      schedule = vfpprod
      schedule += '\nWSEGTABL\n P1 3 3 1 F- EXT NO /\n/'
   for gvf in gvfs:
      k += 1
      for wct in wcts:
         j += 1
         i = -1
         dpv = []
         qv  = []
         for q in qs:
            i += 1
            # translate wct/gvf to rates
            q_gas = q * gvf
            q_liq = q - q_gas
            q_wat = wct * q_liq
            q_oil = q_liq - q_wat
            sgas = q_gas/q
            swat = q_wat/q
            schedule += '''
WCONINJE
 I_OIL OIL   OPEN RESV 1*  %.1f 25000 /
 I_GAS GAS   OPEN RESV 1*  %.1f 25000 /
 I_WAT WATER OPEN RESV 1*  %.1f 25000 /
/
''' % (q_oil, q_gas, q_wat)
            schedule += 'TSTEP\n %i*%.1f /\n' % (ntsteps, dt)
   repl = [('_MUO_1_', str(mu_o)),
           ('_MUO_2_', str(mu_o+0.001)),    # need monotonic behaviour
           ('_MUG_1_', str(mu_g)),
           ('_MUG_2_', str(mu_g+0.00001)),  # need monotonic behaviour
           ('_MUW_', str(mu_w)),
           ('_RHO_O_', str(rho_o)),
           ('_RHO_G_', str(rho_g)),
           ('_RHO_W_', str(rho_w)),
           ('_VFPPDIMS_', '%i 1 %i %i 1 1' % tuple(vfppdims)),
           ('_SCHEDULE_', schedule)]
   _create_qc_datafile(repl, fname, ntsteps)
   return fname

def actions_for_onoff_valves(wellnm, segm_var, lim, sched_file, fname=None):
   '''
   compiles ACTIONS keyword for simulating the effect of on-off inflow valves.
   will shut when some value (typially SWCT) is above some limit, and will
   open up when it goes below.
   assumes 1-1 connection between segments and grid cells.
   uses the raw compdat-data
   #
   # input
   wellnm    : well name
   segm_var  : segment variable to check. typically 'SWCT' or 'SGOR'
   lim       : when to shut/open. typially in [0,1]
   sched_file: where the COMPDAT for the well is defined
   #
   # returns
   string to be put into DATA-file with all ACTIONS keywords
   '''
   dd = DataDeck.DataDeck(sched_file)
   compdats = [x for x in dd.get_raw_data('COMPDAT') if x[0]==wellnm]
   n_segm = len(compdats)
   segments = range(n_segm+2, 2*n_segm+2)
   s = '-- MUST HAVE ACTDIMS  %i /\n/' % (2*n_segm)
   for n,segment in enumerate(segments):
      compdat2 = ' '.join(compdats[n])
      compdat1 = compdat2.replace('OPEN', 'SHUT')
      s1 = '''
ACTIONS
 SHUT_%03i %s %i %s > %.3f /
COMPDAT
 %s /
/
ENDACTIO\n''' % (segment, wellnm, segment, segm_var, lim, compdat1)
      s2 = '''
ACTIONS
 OPEN_%03i %s %i %s < %.3f /
COMPDAT
 %s /
/
ENDACTIO\n''' % (segment, wellnm, segment, segm_var, lim, compdat2)
      s += s1 + s2 
   if fname:
      f = open(fname, 'w')
      f.write(s)
      f.close()
      print fname, 'created'
      return
   return s

def plot_openvalves(pattern, indx1, indx2):
   '''
   for a case with on-off valves (see actions_for_onoff_valves), this
   will analyze the printfile to find how many valves are open at any given
   time
   pattern = file name pattern for printfiles
   indx1   = index for first valve. make sure it is correct!
   indx2   = index for last valve. make sure it is correct!
   '''
   printfiles = UT.glob(pattern, sortit=True)
   pl.figure()
   for i,printfile in enumerate(printfiles):
      if not printfile.endswith('.PRT'):
         print printfile, 'is not a print-file. skipped'
         return
      f = open(printfile)         # print-file
      status = {}
      for n in range(indx1, indx2+1):
         status['%03i'%n] = 1       # start with all open
      openvalves = [len(status)]    # keeps track of number of open valves
      t          = [0]              # time for events (shutting / opening of valves)
      for line in f:
         if 'WILL BE ACTIVATED AT TIME' in line:
            rec = line.split()
            mode, valveno = rec[3].split('_')
            if mode == 'SHUT': status[valveno] = 0 
            else             : status[valveno] = 1 
            t.append(float(rec[-1].replace('*','')))
            openvalves.append(sum(status.values()))
      f.close()
      pl.plot(t, openvalves, 'o--', color=UT.COLOURS[i], label=UT.basename(printfile))
   pl.xlabel('Time [days]')
   pl.ylabel('# open valves [-]')
   pl.grid(True)
   pl.legend(loc='best')
   pl.show()
   return t, openvalves

def unify_results(casenm, mode, first_tstep, last_tstep, cleanup=True):
   '''
   create unified resultfiles
   mode: 'results' or 'summary'
   '''
   mode_id = 2 if mode == 'results' else 4
   f = open('conv.txt', 'w')
   s = '''U
%s
2
%i
%i
%i
N
''' % (casenm, mode_id, first_tstep, last_tstep)
   f.write(s)
   f.close()
   os.system('@convert < %s' % f.name)
   if cleanup: os.unlink(f.name)

class EclipseCoupling(object):
   '''
   A simple class for coupling eclipse to python. it uses READDATA.
   Make sure that:
      - You dont use UNIFOUT
      - You do at least one time-step before entering READDATA
      - Use RPTONLY to force eclipse to write summary-files (eg A1.S000x) only at report-times
        (not at every time-step)
      - There is at least one statement (other than TSTEP or END) after READDATA
      - You provide this class with a simulation state instance. This class *must* have
         casenm                    : just the basename of the simulation
         next_tstep(tstep_no, trn) : returns next tstep. return None when you want to stop.
                                     trn: an open transcript file (or None)
         new_schedule(tstep, trn)  : returns a formatted string that is the schedule-part. tstep is from next_tstep.
                                     typically, this function will read output from eclipse (should use read_summary2 for this)
                                     to adjust some input in the schedule section
                                     trn: an open transcript file (or None)
   #
   Example usage (this naive example will only run a simulation with fixed time-steps)
   #
   class SimState(object):
      def __init__(self, casenm, tstep=1., maxtime=10):
         self.casenm       = casenm
         self.maxtime      = maxtime
         self.tstep        = tstep
         self.simtime      = 0.
      #
      def next_tstep(self, tstep_no):
         self.simtime += self.tstep
         return self.tstep if self.simtime < self.maxtime else None
      #
      def new_schedule(self):
         return '\nTSTEP\n %.2e \n/' % self.tstep
   #
   ss = SimState(casenm)
   ec = ECL.EclipseCoupling(ss, trn=True)
   ec.run()
   ec.unify_results()
   ec.cleanup()
   '''
#
   def __init__(self, simstate, tstep_no=1, trn=False, sleeptm1=0.1, sleeptm2=2.0):
      '''
      simstate: see above
      trn     : use transcript file (boolean)
      tstep_no: must be >= 1, which implies that simulation has been run at least one timestep
                before coupling is activated (i.e. need TSTEP before READDATA). that's how it is:
                we are lagging one step behind the simulator
      sleeptm1: how long to wait before checking again if eclipse has done its part (creating eg. A1.S000x-file)
                default value should be appropriate.
      sleeptm2: wait a bit after detecting the file to make sure eclipse has finished writing to the file
                default value should be appropriate.
      '''
      self.casenm   = simstate.casenm   # for convinience
      self.simstate = simstate
      self.tstep_no = tstep_no
      self.sleeptm1 = sleeptm1
      self.sleeptm2 = sleeptm2
      self.trn = open('%s.trn'%self.casenm, 'w') if trn else None
#
   def _write_schedule_for_readdata(self, txt):
      '''
      writes a schedule file for READDATA
      txt should be a valid eclipse schedule section - or empty
      '''
      fnm =  '%s.I%04i' % (self.casenm, self.tstep_no)     # ala A01.I0002
      f = open(fnm, 'w')
      f.write(txt)
      f.close
      if self.trn:
         msg = 'created file %s' % fnm
         if txt: msg += ' (with text)'
         print >>self.trn, msg
#
   def _wait_for_reportstep(self, tstep_no):
      '''
      make sure eclipse has simulated the given time-step
      '''
      sumfile = '%s.S%04i' % (self.casenm, tstep_no)
      if self.trn: print >>self.trn, 'waiting for', sumfile
      while not os.path.exists(sumfile):
         time.sleep(self.sleeptm1)
      time.sleep(self.sleeptm2)  # make sure eclipse is done...
      if self.trn: print >>self.trn, 'found', sumfile
#
   def run(self, run_eclipse=True, ecl_version='2016.2'):
      '''
      this is the main engine of the EclipseCoupling. does the ping-pong.
      if run_eclipse is True, it will fork out a process to run the eclipse case.
      this is sometimes unstable...
      (TODO??: To make sure eclipse is done with its part, we could check simulation time from the report-file)
      '''
      print '## main loop'
      if not run_eclipse: print 'Time to start eclipse-case'
      while True:
         self._write_schedule_for_readdata('')                                  # empty schedule-file so eclipse will wait for new instructions
         if run_eclipse:
            run_eclipse = False                                                 # must only be done once
            # run eclipse in parallel
            pid = os.fork()
            if pid == 0:                                                        # this is the child process
               lcl_run(self.casenm+'.DATA', version=ecl_version)
               os._exit(0)                                                      # hard stop of the child process
         self._wait_for_reportstep(self.tstep_no)                               # make sure eclipse has finished its time-step
         tstep = self.simstate.next_tstep(self.tstep_no, trn=self.trn)
         if self.trn: print >>self.trn, 'tstep=', tstep
         if tstep:
            self._write_schedule_for_readdata(self.simstate.new_schedule(tstep, trn=self.trn)) # new schedule part to eclipse
         else:
            self._write_schedule_for_readdata('RPTONLY')                        # dummy schedule part to eclipse
         os.system('touch %s.OK'%self.casenm)                                   # a signal to eclipse to pick up schedule file
         print ' # main loop', self.tstep_no, tstep
         self.tstep_no += 1
         if not tstep: break
      print 'DONE'
#
   def unify_results(self):
      unify_results(self.casenm, 'results', 1, self.tstep_no-1)  # -1 since last command in loop is a dummy
      unify_results(self.casenm, 'summary', 1, self.tstep_no-1)
#
   def cleanup(self):
      if self.trn and not self.trn.closed: self.trn.close()
      os.system('cat %s.I0* >! %s.sched' % (self.casenm, self.casenm))     # cat all sched-files created into one file

# aliases
load_rsm = read_rsm
from AtlejgTools.EclipseTools.DataDeck import DataDeck
