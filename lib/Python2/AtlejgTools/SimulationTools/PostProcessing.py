#!/usr/bin/env python

import sys
import time
from pylab import *
import os
import ConfigParser
import glob, re, fnmatch
from datetime import datetime
import AtlejgTools.Utils                             as UT
import AtlejgTools.EclipseTools.Utils                as ECL
import AtlejgTools.Plotting.Utils                    as PU
import AtlejgTools.RCP_project.Valve_characteristics as VC

def _funrst_reader(fname, nx=0,ny=0,nz=0):
   print 'reading', fname
   return ECL.Three_D_Data(UT.basename(fname), int(nx), int(ny), int(nz))

def _init_reader(fname, nx=0,ny=0,nz=0):
   print 'reading', fname
   return ECL.Three_D_Data(UT.basename(fname), int(nx), int(ny), int(nz), file_ext='INIT')


def _date(datestr):
   ''' assumes datestr is like 1972-9-24 '''
   ds = datestr.split('-')
   return datetime(int(ds[0]), int(ds[1]), int(ds[2]))

def _get_parameters(case, strings=[]):
   '''
   accepts 
      - A1 (where A1 is directory, and A1/A1.ini is the inifile)
      - A1 (A1.ini is inifile)
      - A1.ini
   strings is a list of variable names that will be read as strings.
   the rest will be read as floats
   '''
   if os.path.isdir(case)    : inifile = '%s/%s.ini' % (case, case)
   elif case.endswith('.ini'): inifile = case   # for historic reasons
   else                      : inifile = case + '.ini'
   if not os.path.exists(inifile):
      raise Exception("cannot find inifile '%s'" % inifile)
   cf = ConfigParser.ConfigParser()
   cf.read(inifile)
   p = UT.Struct()
   # useful to have numbers as floats
   for opt in cf.options('constants'):
      if   opt in strings:
         # these are strings
         txt = "p.%s = cf.get('constants','%s')" % (opt, opt)
      else:
         # the rest are numbers
         txt = "p.%s = float(cf.get('constants','%s'))" % (opt, opt)
      exec(txt)
   # useful...
   p.name = UT.basename(inifile)
   return p

def _startdate(startdate):
   months = {'JAN':1, 'FEB':2, 'MAR':3, 'APR':4, 'MAY':5, 'JUN':6, 'JUL':7, 'JLY':7, 'AUG':8, 'SEP':9, 'OCT':10, 'NOV':11, 'DEC':12}
   y = int(startdate[2])
   m = months[startdate[1].upper()] # 'JAN' etc.
   d = int(startdate[0])
   return y, m, d


def get_summary(case):
   '''
   case could be a directory, a file-name, or a casename. if it is a filename with
   an extension, this extension is replaced with UNSRMY so that summary may be read.
   if it is a directory, it assumes there is only one summary file in this directory.
   for performance reasons, it will hold startdate in a separate array since the 
   DataDeck operation could be very slow.
   '''
   if os.path.isdir(case):
      localdir = case
      case = '*' # assumes only one summary file in directory
   else:
      localdir = '.'
      case = UT.basename(case)
   cwd = os.getcwd()
   os.chdir(localdir)
   smryfile = UT.glob('%s.*UNSMRY' % case)[0] # accepts UNSRMY or FUNSMRY extension
   if not os.path.exists(smryfile):
      raise Exception("No summary file for case %s" % case)
   s = CM1.add(smryfile, reread=UI.reread)
   os.chdir(cwd)
   if UI.startdate: ymd = UI.startdate       # year,month,date
   else:
      if STARTDATES.has_key(smryfile):
         ymd = STARTDATES[smryfile]
      else:
         d = ECL.DataDeck(case+'.DATA').get_raw_data('START')[0]
         ymd = _startdate(d)
         STARTDATES[smryfile] = ymd
   s.startdate = date2num(datetime(*ymd))
   s.date = s.time + s.startdate             # very useful when using plot_date etc.
   return s

def _fix(s, varnm):
   '''
   automatically give WBHP-P1 if WBHP is asked for
   '''
   if varnm.startswith('W') and not s.separator in varnm: varnm += s.separator+WELLNM
   return varnm

def get_val(varnm, cases, tstep=-1):
   UI.scaler       = 1.    # if you want to scale any variable. use with care!!
   y = []
   if type(cases) is str: cases = (cases,) # make sure we have a list
   for case in cases:
      s = get_summary(case)
      y.append(s.get(_fix(s,varnm))[tstep])
   return y

def _get_params(parnm, cases):
   y = []
   for case in cases:
      if parnm == 'licd':
         val = JOINT_LENGTH/_get_param_value(UT.basename(case), 'valves_pr_joint')
      else:
         val = _get_param_value(UT.basename(case), parnm)
      y.append(val)
   return y

def _icd_segments(case, wellnm):
   s = get_summary(case)
   return s.icd_segments(wellnm)

def _well_segments(case, wellnm):
   s = get_summary(case)
   return s.well_segments(wellnm)

# constants
COLOURS   = UT.COLOURS
GRAVITY   = 9.8
JOINT_LENGTH = 12.

# globals
# useful to have globals as objects
UI               = UT.Struct()  # user input
UI.sortit        = True  # sort list of files alphabetically?
UI.accum         = False # accumulate some data as a function of time?
UI.md_accum      = False # accumulate some data as a function of md?
UI.legend_type   = 1     # 1 :show completion-type and casenm when showing legend
                         # 2 :show casenm only when showing legend
                         # 3 :show completion-type only when showing legend
UI.plot_md       = False # if true, it will plot measured depth in stead of segment-# for segments
UI.plot_kwargs   = None  # typically PP.UI.plot_kwargs = lambda case, cases: {'marker':'o', 'ls':'--'}
UI.legend        = None  # function. should be set externally.
UI.use_in_situ   = False # plot rates at reservoir conditions. to be implemented!!
UI.Bo           = 1.    # formation factor oil. to be done: should not be constant
UI.Bg           = 1.    # formation factor gas. to be done: should not be constant
UI.Bw           = 1.    # formation factor wat. to be done: should not be constant
UI.Rs           = 1.    # resolved gas.         to be done: should not be constant
UI.COMP_VARS    = []      # variables to be plotted with the 'compiled' functionality (like ['FOPR', 'WBHP'])
UI.yscaler      = None    # if you want to scale any variable in 'xplot' and in 'segments'. use with care!!
UI.startdate    = None    # use (yyyy, mm, dd)  if you want to overwrite START in DATA-file
UI.reread       = True    # if True: will reread summary if file has changed since last access
UI.plot_dates   = True    # if True: will use plot_date in xplot when x-variable is TIME
UI.show_fig     = True    # if True: will display figure. set to False if you want to save hardcopy without showing figure
UI.silent       = False   # if True: will be as quiet as possible
UI.segm_length  = None    # if plot_md is True, you need to set this. segm_length = segm_length(casenm)
UI.plot_hist    = False   # if plot_hist is True, it will plot the historical data (e.g. WOPRH) for the last case

# useful for debugging and accessing data outside the program
DBG = UT.Struct()

#CM1 = UT.CacheManager(ECL.read_summary)
CM1 = UT.CacheManager(ECL.read_summary2) # experimental ...
CM2 = UT.CacheManager(_funrst_reader)
CM3 = UT.CacheManager(_init_reader)
y_ = '' # useful for debugging or interactive use

# for performance reasons - avoid calling DataDeck too often. setting the UI.startdate is the fastest
STARTDATES = {}

#icd_segments     = None  # function. should be set externally. usually _icd_segments below is what you want
#well_segments    = None  # function. should be set externally. usually _well_segments below is what you want
icd_segments     = _icd_segments
well_segments    = _well_segments
well_indices     = None   # should be set externally.
WELLNM           = 'none'
get_parameters   = None   # should be set externally. typically create.get_parameters

rcParams['lines.linewidth'] = 4 # 1 is default
rcParams['lines.markersize'] = 12 # 6 is default

def _get_param_value(case, parnm):
   p = get_parameters(case)
   return p.__dict__[parnm]

def _segm2j(case, segm):
   '''
   maps well segment as defined in eclipse to j-index
   '''
   p = get_parameters(case)
   well_j = well_indices(p)[1]
   return segm + well_j - 2

def _j2segm(case, j, use_icd_segm=False):
   '''
   maps j-index to well segment as defined in eclipse
   if use_icd_segm is True, it will make sure you get that one.
   '''
   p = get_parameters(case)
   well_j = well_indices(p)[1]
   xtra = p._nsegm if use_icd_segm else 0
   return j - well_j + 2 + xtra

def _compl_nm(case):
   casenum = re.split('[A-Z]*', case)[1]  # re.split: A132 -> 132
   if casenum[0] == '1': return 'none'
   if casenum[0] == '2': return 'RCP'
   if casenum[0] == '3': return 'ICD'

def _groups(cases, ngroups, index=1):
   '''
   divides the list of cases into ngroups based on the casename index.
   example: cses = (A1, B1, C1, A2, B2, C2) -> (A1, A2), (B1, B2), (C1, C2)
   '''
   groups = []
   for i in range(ngroups): groups.append([])
   for case in cases:
      i = int(re.split('[A-Z]*', case)[1][index-1])   # re.split: A132 -> 132
      groups[i-1].append(case)
   return groups

def _welldp(cases, plot_kwargs=None):
   global y_
   figure()
   red = linspace(0, 1, len(cases))
   n = -1
   for case in cases:
      n += 1
      s = get_summary(case)
      wsegm = well_segments(case, WELLNM)
      dp = s.get('SPR-%s-%i' % (WELLNM, wsegm[-1])) - s.get('SPR-%s-%i' % (WELLNM, wsegm[0]))
      if plot_kwargs: kwargs = plot_kwargs(case, cases)
      else          : kwargs = {'color':(red[n], 0, 1-red[n])}
      plot(s.time, dp, label=os.path.splitext(case)[0], **kwargs)
   ylabel('dp [bar]')
   title('Pressure drop from toe to heel')
   legend(loc='best')
   xlabel('TIME [days]')
   grid(True)
   y_ = dp


def _comparison(varnm, cases, plot_kwargs=None):
   '''
   compares cases[n] to cases[n+1]
   '''
   global y_
   figure()
   plotted = {}
   ymin = Inf; ymax = -Inf
   red = linspace(0, 1, len(cases)/2)
   for n in range(len(cases)/2):
      s1 = get_summary(cases[2*n])
      s2 = get_summary(cases[2*n+1])
      y = interp(s2.time, s1.time, s1.get(varnm))
      y = (y - s2.get(varnm)) / s2.get(varnm) * 100. # want percent
      kwargs = {'color':(red[n], 0, 1-red[n])}
      plot(s2.time, y, label='%s/%s'%(UT.basename(s1.nm), UT.basename(s2.nm)), **kwargs)
   ylabel('rel. change %s [%%]' % (varnm))
   title('comparison of cases')
   legend(loc='best')
   xlabel('TIME [days]')
   grid(True)

def _xplot(varnm1, varnm2, cases, refcase=None, plot_kwargs=None, nolegend=False):
   global y_
   figure(); 
   ymin = Inf; ymax = -Inf
   n = -1
   if refcase: ref = get_summary(refcase)
   red = linspace(0, 1, len(cases))
   for case in cases:
      n += 1
      s = get_summary(case)
      varnm1 = _fix(s,varnm1)
      varnm2 = _fix(s,varnm2)
      if varnm1 == 'PVI':
         y1 = _pvi(UT.basename(case))
      else:
         y1 = s.get(varnm1)
      y2 = s.get(varnm2)
      if refcase:
         y2 = interp(ref.time, s.time, y2)
         y2 = (y2 - ref.get(varnm2)) / ref.get(varnm2) * 100. # want percent
         y1 = ref.time
      if   plot_kwargs              : kwargs = plot_kwargs(case, cases)
      elif len(cases) > len(COLOURS): kwargs = {'color':(red[n], 0, 1-red[n]), 'linestyle':'-', 'marker':None}
      else                          : kwargs = {'color':COLOURS[n], 'linestyle':'-', 'marker':None}
      if varnm2.startswith('WMCTL'):
         kwargs['linestyle'] = 'None'
         kwargs['marker']    = 'o'
         kwargs['markeredgewidth'] = 0
      if UI.accum: y2 = UT.cumulative(y1, y2)
      if UI.yscaler:
         if not UI.silent: print 'Warning: UI.yscaler = ', UI.yscaler
         y2 *= UI.yscaler
      if varnm1 == 'TIME' and UI.plot_dates:
         y1 += s.startdate
         plotfunc = plot_date
      else:
         plotfunc = plot
      plotfunc(y1, y2, label=os.path.splitext(case)[0], **kwargs)
   #
   if UI.plot_hist:
      varnmh = _fix(s, varnm2.split(s.separator)[0]+'H')
      yh = s.get(varnmh)
      plotfunc(y1, yh, 'k--', label=os.path.splitext(case)[0]+' (hist)')
   #
   if refcase:
      ylabel('rel. change %s [%%]' % (varnm2))
      titl = 'compare to case ' + UT.basename(refcase)
   else:
      ylabel('%s [%s]' % (varnm2, s.unit(varnm2)))
      titl = varnm2
   if UI.accum: titl += ' (accumulated)'
   if varnm1 == 'TIME' and UI.plot_dates:
      gcf().autofmt_xdate()   # beautify
   elif varnm1 == 'PVI':
      xlabel('PVI [-]')
   else:
      xlabel('%s [%s]' % (varnm1, s.unit(varnm1)))
   if not nolegend: legend(loc='best')
   title(titl)
   grid(True)
   if varnm2.startswith('WMCTL'):
      yticks(range(-7,14), ('bhp (grp)', 'thp (grp)', 'resv (grp)', 'lrat (grp)',   # grp = group control
                           'grat (grp)', 'wrat (grp)', 'orat (grp)',
                           'shut/stop',
                           'orat', 'wrat', 'grat',
                           'lrat', 'resv', 'thp', 'bhp',
                           ' ', 'crat', ' ', 'GOR Penalty', 'drawdown', 'ngl'))

def _meanplot(varnm, cases):
   figure(); 
   y = 0
   t = []
   for case in cases:
      s = get_summary(case)
      if len(t) == 0:
         t = s.get('TIME')   # use first case
         varnm = _fix(r,varnm)
      y += interp(t, s.time, s.get(varnm))
   if UI.plot_dates:
      t += s.startdate
      plotfunc = plot_date
   else:
      plotfunc = plot
   y /= len(cases)
   plotfunc(t, y, '-')
   ylabel('%s [%s]' % (varnm, s.unit(varnm)))
   title('Mean values')
   grid(True)
   DBG.t, DBG.y = t, y

def _wbarplot(varnm, wname_pattern, case, t0, shift=0, mark_highest=True, rel_height=False):
   '''
   cross-plot for wells. useful in case of having many wells. wname_pattern should be a reg-exp pattern'
   '''
   s = get_summary(case)
   ptrn = '%s%s%s' % (varnm, ECL.SEPARATOR, wname_pattern)
   varnms = _get_varnms(wname_pattern, varnm, case)
   figure()
   yvals = []
   maxval = -1.
   i = -1
   for varnm in varnms:
      i += 1
      val = s.get(varnm)[UT.find_index(s.time, t0, eps=t0/40.)]
      yvals.append(val)
      if abs(val) > maxval:
         maxval = abs(val)
         ind    = i
   colors = ['b' for i in arange(len(yvals))]
   if rel_height:
      yvals = array(yvals)
      yvals /= maxval
      yvals *= 100 # percent
   if mark_highest: colors[ind] = 'r'
   bar(arange(len(yvals))+0.5-shift, yvals, color=colors)
   xticks(arange(len(yvals))+1,[x.split(ECL.SEPARATOR)[1] for x in varnms], rotation=90)
   if rel_height:
      ylabel('relative [%]')
   else:
      ylabel('%s [%s]' % (varnm, s.unit(varnm)))
   title('%s @ %i days' % (varnm, int(t0)))

def _get_varnms(wname_pattern, varnm, case):
   s = get_summary(case)
   ptrn = '%s%s%s' % (varnm, ECL.SEPARATOR, wname_pattern)
   varnms = filter(lambda v: fnmatch.fnmatch(v, ptrn), s.varnames())
   varnms.sort()
   return varnms

def _wxplot(varnm1, varnm2, wname_pattern, case, plot_kwargs=None, nolegend=False):
   '''
   cross-plot for wells. useful in case of having many wells. wname_pattern should be a reg-exp pattern'
   '''
   global y_
   figure(); 
   plotted = {}
   ymin = Inf; ymax = -Inf
   s = get_summary(case)
   ptrn = '%s%s%s' % (varnm2, ECL.SEPARATOR, wname_pattern)
   varnms = _get_varnms(wname_pattern, varnm2, case)
   red = linspace(0, 1, len(varnms))
   n = -1
   for varnm in varnms:
      n += 1
      y1 = s.get(varnm1)
      y2 = s.get(varnm)
      if   plot_kwargs               : kwargs = plot_kwargs(case, cases)
      elif len(varnms) > len(COLOURS): kwargs = {'color':(red[n], 0, 1-red[n])}
      else                           : kwargs = {'color':COLOURS[n]}
      if varnm2.startswith('WMCTL'):
         kwargs['linestyle'] = 'None'
         kwargs['marker']    = 'o'
         kwargs['markeredgewidth'] = 0
      if UI.accum: y2 = UT.cumulative(y1, y2)
      plot(y1, y2, label=varnm, **kwargs)
   #
   ylabel('%s [%s]' % (varnm2, s.unit(varnms[0])))
   titl = case
   if UI.accum: titl += ' (accumulated)'
   xlabel('%s [%s]' % (varnm1, s.unit(varnms[0])))
   if not nolegend: legend(loc='best')
   title(titl)
   grid(True)
   if varnm2.startswith('WMCTL'):
      yticks(range(-7,8), ('bhp (gr)', 'thp (gr)', 'resv (gr)', 'lrat (gr)',   # gr = group control
                           'grat (gr)', 'wrat (gr)', 'orat (gr)',
                           'shut/stop',
                           'orat', 'wrat', 'grat',
                           'lrat', 'resv', 'thp', 'bhp'))

def _segm_length(case):
   p = get_parameters(case)
   segments = well_segments(case, WELLNM)
   return p._dy * arange(len(segments))

def _segment_data_at_given_time(varnm, cases, t0, plot_kwargs=None, adjust_to_zero=False, use_icd_segm=True, accum=False):
   global y_
   figure()
   red = linspace(0., 1., len(cases)) # color scaler
   n = 0
   for case in cases:
      case = UT.basename(case)
      if use_icd_segm: segments = icd_segments(case, WELLNM)
      else           : segments = well_segments(case, WELLNM)
      s = get_summary(case)
      if t0 > 0: tix = UT.find_closest_index(s.time, t0)
      else     : tix = len(s.time) # last one
      if accum:
         tixs = range(tix)   # time indices
         m = s.get_segm_data(varnm, WELLNM, segments)[:, tixs]
         y_ = array([UT.cumulative(s.time[tixs], m[i,:])[-1] for i in range(len(segments))])
      elif varnm == 'SGVF':
         y_ = _gvf_segmentdata(r, segments, tix)
      else:
         y_ = s.get_segm_data(varnm, WELLNM, segments)[:, tix]
      if   plot_kwargs              : kwargs = plot_kwargs(case, cases)
      elif len(cases) > len(COLOURS): kwargs = {'color':(red[n], 0, 1-red[n])}
      else                          : kwargs = {'color':COLOURS[n]}
      if adjust_to_zero: y_ -= min(y_)
      if UI.md_accum: y_ = cumsum(y_)
      if UI.plot_md:
         sl = UI.segm_length(case)
         x = cumsum(sl)
         y_ /= sl        # flow pr. meter
         xlbl = 'md [m]'
      else:
         x = segments + 1  # 26/4-17: removed adjustment of x-variable
         xlbl = 'segment# [-]'
      if UI.yscaler:
         if not UI.silent: print 'Warning: UI.yscaler = ', UI.yscaler
         y_ *= UI.yscaler
      plot(x, y_, label=os.path.splitext(case)[0], **kwargs)
      n += 1
   #
   legend(loc='best')
   ylabel('%s [%s]' % (varnm, s.unit('%s%s%s%s%i'%(varnm, ECL.SEPARATOR, WELLNM, ECL.SEPARATOR, segments[0]))))
   xlabel(xlbl)
   ylim(0, ylim()[1])
   grid(True)
   titl = '%s @ t= %i days' % (varnm, t0)
   if UI.md_accum and accum: titl += ' (accumulated time and space)'
   elif UI.md_accum        : titl += ' (accumulated in space)'
   elif accum              : titl += ' (accumulated in time)'
   title(titl)

def _connection_data_at_given_time(varnm, cases, t0, dx, plot_kwargs=None):
   if WELLNM == 'none':
      print 'Need to set WELLNM!'
      return
   global y_
   figure()
   red = linspace(0., 1., len(cases)) # color scaler
   n = 0
   for case in cases:
      dd = ECL.DataDeck(case)
      ijks = dd.raw2values(dd.get_raw_data('COMPDAT', get_all=True), 2,4, identif=WELLNM, match_col=1)
      case = UT.basename(case)
      s = get_summary(case)
      t0_indx = UT.find_closest_index(s.time, t0)
      if UI.accum: tixs = range(t0_indx)   # time indices
      else       : tixs = (t0_indx,)
      y_ = 0
      sep = ECL.SEPARATOR  # for convinience
      if varnm == 'CLFR':
         m  =  array([s.get('%s%s%s%s%i,%i,%i'%('COFR', sep, WELLNM, sep, ijk[0],ijk[1],ijk[2])) for ijk in ijks]).T
         m  += array([s.get('%s%s%s%s%i,%i,%i'%('CWFR', sep, WELLNM, sep, ijk[0],ijk[1],ijk[2])) for ijk in ijks]).T
      else:
         m = array([s.get('%s%s%s%s%i,%i,%i'%(varnm, sep, WELLNM, sep, ijk[0],ijk[1],ijk[2])) for ijk in ijks]).T
      for tix in tixs:
         if UI.accum:
            if tix > 0: dt = s.time[tix] - s.time[tix-1]
            else      : dt = s.time[tix]
         else: dt = 1    # dont scale it
         y_ += dt*m[tix,:]
      if   plot_kwargs              : kwargs = plot_kwargs(case, cases)
      elif len(cases) > len(COLOURS): kwargs = {'color':(red[n], 0, 1-red[n])}
      else                          : kwargs = {'color':COLOURS[n]}
      x = dx*arange(len(y_))
      step(x, y_, label=os.path.splitext(case)[0], **kwargs)
      n += 1
   #
   legend(loc='best')
   ylabel('%s for well %s'%(varnm,WELLNM))
   xlabel('md [m]')
   ylim(0, ylim()[1])
   grid(True)
   titl = '%s @ t= %i days' % (varnm, s.time[t0_indx])
   if UI.md_accum and UI.accum: titl += ' (accumulated time and space - to 1st order)'
   elif UI.md_accum           : titl += ' (accumulated in space)'
   elif UI.accum              : titl += ' (accumulated in time - to 1st order)'
   title(titl)

def _compiled(cases, segm1, segm2, plot_kwargs=None):
   if not UI.COMP_VARS:
      print 'You need to set UI.COMP_VARS'
      return
   os.system('rm -f xxx*.png')
   figsize = rcParams['figure.figsize']
   rcParams['figure.figsize'] = 7.5, 4
   n = 0
   for varnm in UI.COMP_VARS:
      if varnm.startswith('S'):
         _xplot('TIME', '%s%s%s%s%i'%(varnm, ECL.SEPARATOR, WELLNM, ECL.SEPARATOR, segm1), cases, plot_kwargs=UI.plot_kwargs, nolegend=False)
         n += 1
         savefig('xxx%03i_%s_%i.png' % (n, varnm, segm1))
         close()
         _xplot('TIME', '%s%s%s%s%i'%(varnm, ECL.SEPARATOR, WELLNM, ECL.SEPARATOR, segm2), cases, plot_kwargs=UI.plot_kwargs, nolegend=False)
         n += 1
         savefig('xxx%03i_%s_%i.png' % (n, varnm, segm2))
         close()
      elif varnm.startswith('W'):
         _xplot('TIME', '%s%s%s'%(varnm, ECL.SEPARATOR, WELLNM), cases, plot_kwargs=UI.plot_kwargs, nolegend=False)
         n += 1
         savefig('xxx%03i_%s.png' % (n, varnm))
         close()
      else:
         _xplot('TIME', varnm, cases, plot_kwargs=UI.plot_kwargs, nolegend=False)
         n += 1
         savefig('xxx%03i_%s.png' % (n, varnm))
         close()
   # put together to one figure
   prefix = cases[0][0] # typically 'A' (from A1.DATA)
   fname = '%s_collage.png' % ('_vs_'.join([UT.basename(case) for case in cases]))
   UT.make_collage(UT.glob('xxx*.png', sortit=True), len(UI.COMP_VARS)/2+1, 2, newfile=fname)
   os.system('gthumb %s &' % fname)
   rcParams['figure.figsize'] = figsize
   #os.system('rm -f xxx*.png')

def _prod_pr_part(varnm, case, pr_segment, accum):
   figure()
   titl = ('%s %s' % (case, varnm))
   s = get_summary(case)
   p = get_parameters(case)
   s1, s2 = [_j2segm(case, p.channel_j1, True), _j2segm(case, p.channel_j2, True)]
   segm_toe  = arange(s2+1, 2*p._nsegm+1)
   segm_chnl = arange(s1, s2+1)
   segm_heel = arange(p._nsegm+2, s1)
   if False: # set manually
      shift = 10
      segm_toe  = arange(s2+shift, s2+shift+1)
      segm_heel  = arange(s1-shift, s1-shift+1)
      titl += ' shift= %i'%shift
   y1 = 0; y2 = 0; y3 = 0;
   for segm in segm_toe:
      y1 += s.get_segm_data(varnm, WELLNM, [segm])
   for segm in segm_chnl:
      y2 += s.get_segm_data(varnm, WELLNM, [segm])
   for segm in segm_heel:
      y3 += s.get_segm_data(varnm, WELLNM, [segm])
   if pr_segment:
      y1 /= len(segm_toe)
      y2 /= len(segm_chnl)
      y3 /= len(segm_heel)
   if accum:
      y1 = UT.cumulative(s.time, y1)
      y2 = UT.cumulative(s.time, y2)
      y3 = UT.cumulative(s.time, y3)
   plot(s.time, y1, 'y', label='toe')
   plot(s.time, y2, 'r', label='channel')
   plot(s.time, y3, 'c', label='heel')
   legend(loc='best')
   if varnm == 'SLFR': ylbl = 'SLFR [SM3/DAY]'
   else: ylbl = '%s [%s]' % (varnm, s.unit('%s-%s-1'%(varnm, WELLNM)))
   ylabel(ylbl)
   xlabel('time [DAYS]')
   grid(True)
   if accum: titl += ' accumulated'
   if pr_segment: titl += ' (pr. segment)'
   title(titl)

def _gvf_segmentdata(r, segments, ind):
   o = s.get_segm_data('SOFR',  WELLNM, segments, scaler=UI.Bo)[:, ind]
   g = s.get_segm_data('SGFRF', WELLNM, segments, scaler=UI.Bg)[:, ind]
   w = s.get_segm_data('SWFR',  WELLNM, segments, scaler=UI.Bw)[:, ind]
   return g/(g+o+w)

def _segment_data_for_given_case(varnm, case, times, adjust_to_zero=False, onefig=True, use_icd_segm=True):
   global y_
   s = get_summary(case)
   if use_icd_segm: segments = icd_segments(case, WELLNM)
   else           : segments = well_segments(case, WELLNM)
   figs = []
   if onefig:
      f = figure()
      figs.append(f.number)
   red = linspace(0., 1., len(times)) # color scaler
   n = 0
   for t in times:
      # requested values
      ind = UT.find_closest_index(s.time, t)
      t0 = s.time[ind]
      if varnm == 'SGVF':
         y_ = _gvf_segmentdata(r, segments, ind)
         ylbl = 'SGVF [-]'
      else:
         y_ = s.get_segm_data(varnm, WELLNM, segments)[:, ind]
         ylbl = '%s [%s]' % (varnm, s.unit('%s-%s-%i'%(varnm, WELLNM, segments[0])))
      if UI.yscaler:
         if not UI.silent: print 'Warning: UI.yscaler = ', UI.yscaler
         y_ *= UI.yscaler
      if adjust_to_zero: y_ -= min(y_)
      if not onefig:
         f = figure()
         figs.append(f.number)
      plot(segments+1, y_, label='%i days'%t0, color=(red[n],0,1-red[n]))  # 26/4-17: removed adjustment of x-variable
      n += 1
   #
   titl = '%s - %s' % (varnm, os.path.splitext(case)[0])
   for figno in figs:
      figure(figno)
      legend(loc='best')
      ylabel(ylbl)
      xlabel('segment# [-]')
      grid(True)
      if adjust_to_zero: titl += ' (adjusted to 0-level)'
      title(titl)

def _contours(varnm, cases, accum, zmin=None, zmax=None, t_end=None, relative=False, use_icd_segm=True):
   for case in cases:
      if use_icd_segm: segments = icd_segments(case, WELLNM)
      else           : segments = well_segments(case, WELLNM)
      s = get_summary(case)
      if UI.plot_md: segm_length_func = lambda segments: cumsum(UI.segm_length(case))
      else         : segm_length_func = None
      ax = s.contour_plot(varnm, WELLNM, segments, zmin=zmin, zmax=zmax, accum=accum, relative=relative, Bg=UI.Bg, Rs=UI.Rs, segm_length_func=segm_length_func)
      titl = '%s %s %s' % (os.path.splitext(case)[0], WELLNM, varnm)
      if accum: titl += ' accum.'
      if varnm == 'SFFR': titl += ' (in situ)'
      ax.title.set_text(titl)
      if t_end: ax.set_ylim(0, t_end)

def _is_shut(r):
   is_shut = (s.get('WMCTL-X21')[-1] == 0)
   if is_shut: print s.nm, 'is closed'
   return is_shut

def _comparewells(varnm, cases, well1, well2):
   t = None
   n1 = n2 = 0
   y1 = y2 = 0.
   for case in cases:
      if well1 in case:
         s = get_summary(case)
         varnm = _fix(r,varnm)
         if t is None:
            t = s.time
            y1 = zeros(t.shape)
         if _is_shut(r): continue
         n1 += 1
         y1 += interp(t, s.time, s.get(varnm))
      elif well2 in case:
         s = get_summary(case)
         if t is None:
            t = s.time
            y2 = zeros(t.shape)
         n2 += 1
         y2 += interp(t, s.time, s.get(varnm))
   if n1 > 0: y1 /= float(n1)
   if n2 > 0: y2 /= float(n2)
   figure()
   plot(t, y1, label=well1)
   plot(t, y2, label=well2)
   xlabel('TIME [days]')
   ylabel(s.unit(varnm))
   grid(True)
   legend(loc='best')

def _barplot_vals(varnm, cases, t0=-1, shift=0, mark_highest=True, rel_height=False):
   '''
   if t0<=0, it will get the latest timestep
   '''
   figure()
   yvals = []
   maxval = -1.
   i = -1
   for case in cases:
      i += 1
      s = get_summary(case)
      varnm = _fix(r,varnm)
      indx = UT.find_index(s.time, t0, eps=t0/40.) if t0 > 0 else -1
      val = s.get(varnm)[indx]
      yvals.append(val)
      if abs(val) > maxval:
         maxval = abs(val)
         ind    = i
   colors = ['b' for i in arange(len(yvals))]
   if rel_height:
      yvals = array(yvals)
      yvals /= maxval
      yvals *= 100 # percent
   if mark_highest: colors[ind] = 'r'
   bar(arange(len(yvals))+0.5-shift, yvals, color=colors)
   xticks(arange(len(yvals))+0.5,[UT.basename(x) for x in cases], rotation=90)
   if rel_height:
      ylabel('%s [%%] @ %i days' % (varnm, int(t0)))
   else:
      ylabel('%s [%s] @ %i days' % (varnm, s.unit(varnm), int(t0)))

def _standard_deviation_along_well(varnm, wellnm, cases, tstep):
   figure()
   yvals = []
   for case in cases:
      s = get_summary(case)
      segments = icd_segments(case, WELLNM)
      yvals.append(std(s.get_segm_data(varnm, wellnm, segments, tsteps=(tstep,))))
   bar(arange(len(yvals))+0.5, yvals)
   xticks(arange(len(yvals))+1,[UT.basename(x) for x in cases], rotation=90)
   unit = s.unit(varnm+'-'+wellnm+'-'+str(segments[0]))
   ylabel('%s [%s]' % (varnm, unit))
   title('Standard deviation @ %i days' % s.time[tstep])

def _water_breakthrough_time(case, segment, wct_limit=0.05):
   s = get_summary(case)
   wct = s.get_segm_data('SWCT', WELLNM, [segment])
   try:
      indx = (wct>wct_limit).nonzero()[0][0]
   except:
      return Inf
   return s.time[indx]

def _special_plot1(varnm, cases, j):
   '''
   plot x=dp_sand/dp_well versus y=FOPT 
   need a lot of cases with varying well_diam to get
   a variation in x.
   '''
   cases.sort()
   x = []; y = []
   for case in cases:
      s = get_summary(case)
      varnm = _fix(r,varnm)
      dps = _dp_sand(UT.basename(case), y_indices)
      _welldp((case,))
      dpw = array(y_)
      ts  = _find_drawdown_sand(UT.basename(case), 0, j)[3]
      dps = interp(s.time, ts, dps)
      x.append(mean(dps[1:]/dpw[1:]))
      y.append(s.get(varnm)[-1])
   #close('all')
   figure()
   plot(x, y, '-*')
   ylabel(varnm)
   xlabel('dp_sand / dp_well [-]')
   grid(True)

def _variable_vs_param(varnm, parnm, cases, rel_changes, tstep=-1):
   global y_
   figure()
   x = []; y_ = []
   for case in cases:
      if parnm == 'licd':
         xval = JOINT_LENGTH/_get_param_value(UT.basename(case), 'valves_pr_joint')
      else:
         xval = _get_param_value(UT.basename(case), parnm)
      x.append(xval)
      s = get_summary(case)
      varnm = _fix(r,varnm)
      y_.append(s.get(varnm)[tstep])
   if rel_changes:
      y_ = array(y_)/max(y_) * 100
      ylabel(varnm + '[%]')
   else:
      ylabel(varnm)
   plot(x, y_, '-*')
   xlabel(parnm)
   grid(True)
   title('%s @ t = %i days' % (varnm, s.time[tstep]))

def _variable_vs_variable(varnm1, varnm2, t0, cases):
   global y_
   figure()
   x = []; y_ = []
   for case in cases:
      s = get_summary(case)
      varnm1 = _fix(r,varnm1)
      varnm2 = _fix(r,varnm2)
      x.append(s.get_interp_value(varnm1, t0))
      y_.append(s.get_interp_value(varnm2, t0))
   plot(x, y_, '-*')
   xlabel(varnm1)
   ylabel(varnm2)
   grid(True)
   title('@ t = %i days' % t0)

def _segments_diff(cases, varnm, segm1, segm2, plot_kwargs=None):
   figure()
   red = linspace(0, 1, len(cases))
   n = 0
   for case in cases:
      s = get_summary(case)
      y1 = s.get('%s-%s-%i' % (varnm, WELLNM, segm1))
      y2 = s.get('%s-%s-%i' % (varnm, WELLNM, segm2))
      if plot_kwargs: kwargs = plot_kwargs(case, cases)
      else          : kwargs = {'color':(red[n], 0, 1-red[n])}
      plot(s.time, y1-y2, label=os.path.splitext(case)[0], **kwargs)
      n += 1
   xlabel('TIME [days]')
   ylabel('%s [%s]' % (varnm, s.unit('%s-%s-%i' % (varnm, WELLNM, segm1))))
   grid(True)
   title('delta %s between segments %i and %i' % (varnm, segm1, segm2))
   legend(loc='best')

def _wbt(cases, parnm, y_indices, relative=False):
   figure()
   x = []
   # get absissa
   for case in cases:
      if parnm == 'dpvalve':
         x.append(_mean_dpvalve([case], -1)[0])
      else:
         x.append(_get_param_value(UT.basename(case), parnm))
   # plot one line for each j (or difference to first j, if relative is True)
   n = -1
   for j in y_indices:
      n += 1
      wbt = []
      for case in cases:
         s = get_summary(case)
         segm = _j2segm(UT.basename(case),j,True)
         wbt.append(_water_breakthrough_time(case, segm))
      if n == 0 and relative:
         wbt0 = array(wbt)
         continue
      if relative: y = array(wbt) - wbt0
      else:        y = wbt
      plot(x, y, '-*', label='j= %i'%j)
   legend(loc='best')
   if relative:
      ylabel('difference in water break through time [days]')
      title('comparing to j= %i' % y_indices[0])
   else :
      ylabel('water break through time [days]')
   xlabel(parnm)
   grid(True)

def _flowr_func(q_tot, nsegments, perm1, perm2, expon):
   return q_tot / float(nsegments) / (1. + (perm2/perm1)**expon)

def _find_drawdown_sand2(case, tstep, y_indx):
   '''
   this one takes away the static pressure.
   '''
   p   = get_parameters(case)
   r3  = CM2.add(case+'.UNRST', nx=p.nx, ny=p.ny, nz=p.nz)  # 3D data
   r3i = CM3.add(case+'.INIT', nx=p.nx, ny=p.ny, nz=p.nz)   # 3D data
   r   = get_summary(case)
   #
   pr = r3.get('PRESSURE', tstep)[:,y_indx-1,:]
   dz = r3i.get('DZ', 0)[0,y_indx-1,:]
   well_i, _, well_k = well_indices(p)
   rho_g = r3.get('GAS_DEN', tstep)[:,y_indx-1,:]
   rho_o = r3.get('OIL_DEN', tstep)[:,y_indx-1,:]
   rho_w = r3.get('WAT_DEN', tstep)[:,y_indx-1,:]
   sg = r3.get('SGAS', tstep)[:,y_indx-1,:]
   sw = r3.get('SWAT', tstep)[:,y_indx-1,:]
   so = 1 - sg - sw
   # calclulate static pressure
   static   = zeros(pr.shape)
   rho = rho_g*sg + rho_o*so + rho_w*sw
   dp_static = rho*GRAVITY*dz / 100000. # Pa -> bar
   for k in range(1, len(dz)):
      static[:,k] = static[:,k-1] + dp_static[:,k]
   # calculate effective drawdown
   tstep2  = UT.find_index(s.time, r3.time[tstep], 1.)
   p_segm  = s.get('SPR-%s-%i'%(WELLNM,well_segments(case,WELLNM)[y_indx-1]))[tstep2]
   dp_segm = s.get('SPRD-%s-%i'%(WELLNM,icd_segments(case,WELLNM)[y_indx-1]))[tstep2]
   drawdown = pr - p_segm - dp_segm + (static[well_i-1,well_k-1]-static)
   return drawdown

def _find_drawdown_sand(case, tstep, y_indx):
   '''
   finds dradown in x-section along the well. 
   finds the drawdown as the difference between reservoir pressure and the
   reservoir pressure in the cell where the well is located for each section.
   dont care about hydraulic pressure or the angle...
   y_indx is y-index, i.e. "j". (eclipse numbered, so starts with j=1)
   '''
   p = get_parameters(case)
   r3  = CM2.add(case+'.UNRST', nx=p.nx, ny=p.ny, nz=p.nz)  # 3D data
   r3i = CM3.add(case+'.INIT', nx=p.nx, ny=p.ny, nz=p.nz)   # 3D data
   #
   pr     = r3.get('PRESSURE', tstep)[:,y_indx-1,:]
   z     = r3i.get('DEPTH', 0)[0,y_indx-1,:]
   dx    = r3i.get('DX', 0)[:,y_indx-1,0]
   dz    = r3i.get('DZ', 0)[0,y_indx-1,:]
   x     = zeros(len(dx))
   for i in arange(1, len(dx)):
      x[i] = x[i-1] + 0.5*(dx[i-1]+dx[i])
   #drawdown = pr[:,::-1] # turn it
   well_i, _, well_k = well_indices(p)
   # subtract static pressure
   rho_g = r3.get('GAS_DEN', tstep)[:,y_indx-1,:]
   rho_o = r3.get('OIL_DEN', tstep)[:,y_indx-1,:]
   rho_w = r3.get('WAT_DEN', tstep)[:,y_indx-1,:]
   sg = r3.get('SGAS', tstep)[:,y_indx-1,:]
   sw = r3.get('SWAT', tstep)[:,y_indx-1,:]
   so = 1 - sg - sw
   static   = zeros(pr.shape)
   drawdown = zeros(pr.shape)
   g = GRAVITY/100000.  # include conversion to Bar
   for k in range(len(dz)-1, -1, -1):
      for i in range(len(x)):
         static[i,k] = rho_g[i,k]*sg[i,k]*dz[k]*g + rho_o[i,k]*so[i,k]*dz[k]*g + rho_w[i,k]*sw[i,k]*dz[k]*g
         if k < len(dz)-1: static[i,k] += static[i,k+1]
         drawdown[i,k] = pr[i,k] - pr[well_i-1, well_k-1] + static[i,k]
   #if tstep == 30: stop
   return (drawdown, x, z, r3.time, static)

def _show_xsection(case, z_indx, varnm, welldepth=0):
   '''
   show animation of xsection of the reservoir
   '''
   p = get_parameters(case)
   r3  = CM2.add(case+'.UNRST', nx=p.nx, ny=p.ny, nz=p.nz)  # 3D data
   r3i = CM3.add(case+'.INIT', nx=p.nx, ny=p.ny, nz=p.nz)   # 3D data
   fig = figure()
   for tstep in arange(len(r3.time)):
      dx = r3i.get('DX', 0)[:,0,z_indx]
      x  = zeros(len(dx))
      for i in arange(1, len(dx)): x[i] = x[i-1] + 0.5*(dx[i-1]+dx[i])
      dy = r3i.get('DY', 0)[0,:,z_indx]
      y  = zeros(len(dy))
      for i in arange(1, len(dy)): y[i] = y[i-1] + 0.5*(dy[i-1]+dy[i])
      only_one_step = False
      if varnm in ['PERMX','PERMY','PERMZ','PORO']:
         v = r3i.get(varnm, 0)[:,:,z_indx]
         only_one_step = True
      elif varnm == 'SOIL':
         sg = r3.get('SGAS', tstep)[:,:,z_indx]
         sw = r3.get('SWAT', tstep)[:,:,z_indx]
         v = 1. - sg - sw
      else:
         v = r3.get(varnm, tstep)[:,:,z_indx]
      v = v[::-1,:] # turn it
      ax = PU.contourf(y[::-1], x, v, fig=fig)
      #ax.axis('equal')
      ax.set_title('%s. %s. %i days. z = %i' % (case, varnm, r3.time[tstep], z_indx))
      ax.set_xlabel('y [m]')
      ax.set_ylabel('x [m]')
      fname = '%s_%s_%02i-%02i.png' % (case, varnm, z_indx, tstep)
      savefig(fname)
      print 'saving file', fname
      if only_one_step: break

def _drawdown_contours(case, tstep, segno):
   drawdown, x, z, t = _find_drawdown_sand(case, tstep, segno)
   ax = PU.contourf(x, z, drawdown.T, zmin=None, zmax=None)
   ax.axis('equal')
   ax.set_title('%s. Drawdown [bar]. %i days. Segment %i' % (case, t[tstep], segno))
   ax.set_xlabel('x-section length [m]')
   ax.set_ylabel('depth [m]')

def _grid_numbers(case):
   p = get_parameters(case)
   nx = int(p.nx)
   ny = int(p.ny)
   nz = int(p.nz)
   well_k = well_indices(p)[2]
   return (nx, ny, nz, well_k)

def _drawdown_profiles(case, tsteps, y_indx, relative_vals):
   figure()
   well_k = _grid_numbers(case)[3]
   for tstep in tsteps:
      drawdown, x, z, t = _find_drawdown_sand(case, tstep, y_indx)
      drdwn_x = drawdown[:, well_k-1]
      if relative_vals: drdwn_x = drdwn_x / max(drdwn_x)
      plot(x, drdwn_x, label='%i days' % t[tstep])
   xlabel('distance [m]')
   if relative_vals: ylabel('relative drawdown [-]')
   else            : ylabel('drawdown [bar]')
   grid(True)
   legend(loc='best')
   title('%s. Drawdown (sandface only). j= %i.' % (case, y_indx))

def _drawdown(case, y_indices):
   '''
   plots dp_valve and dp_sand (from outer boundary to well) against time
   for the given cells
   '''
   figure()
   p = get_parameters(case)
   well_k = _grid_numbers(case)[3]
   n = -1
   for j in y_indices:
      n += 1
      # sand drawdown and pressure support
      t = _find_drawdown_sand(case, 0, j)[3]
      dpsand = []
      p_support = []
      for tstep in range(len(t)):
         p_support.append(CM2.get(case+'.UNRST').get('PRESSURE', tstep)[p.nx-1,j-1,0])
         dpsand.append((_find_drawdown_sand(case, tstep, j)[0])[-1, well_k-1])
      plot(t, p_support, COLOURS[n]+'.-', label='press. supp. j= %i'%j)
      plot(t[1:], dpsand[1:], COLOURS[n]+'-*', label='dP sand j= %i'%j)  # skip first point which is 0 (only confusing)
      # valve drawdown
      s = get_summary(case)
      dp_valve = s.get_segm_data('SPRD', WELLNM, [_j2segm(case,j,True)])
      plot(s.time, dp_valve, COLOURS[n]+'--', label='dP valve j= %i'%j)
      # total drawdown
      dp_tot = dp_valve + interp(s.time, t, dpsand)
      plot(s.time, dp_tot, COLOURS[n]+'-.', label='dP tot j= %i'%j)
   # plot BHP as well
   plot(s.time, s.get('WBHP-'+WELLNM), 'k.', label='BHP')
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('%s. Drawdowns.' % (case))

def _dp_sand(case, y_indices, z_indx):
   '''
   plots dp_sand for the given sections. it is found by dp from well-cell
   to cell given by z_indx and y_indices
   '''
   global y_
   figure()
   p = get_parameters(case)
   well_i, _, well_k = well_indices(p)
   n = -1
   for j in y_indices:
      n += 1
      # sand drawdown and pressure support
      t = _find_drawdown_sand(case, 0, j)[3]
      dpsand = []
      for tstep in range(1,len(t)-1):  # for some reason - the last step is not always found for all variables
         dpsand.append(_find_drawdown_sand2(case, tstep, j)[well_i-1, z_indx-1])
      plot(t[1:-1], dpsand, COLOURS[n]+'-*', label='j= %i'%j)
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('%s. dP sand. k=%i' % (case, z_indx))
   y_ = dpsand
   return dpsand

def _dp_sand2(unrstfiles, y_indx, z_indx, plot_kwargs=None):
   '''
   plots dp_sand for the given cases. it is found by dp from well-cell
   to cell given by z_indx and y_indx
   '''
   figure()
   n = -1
   if not plot_kwargs: red = linspace(0., 1., len(unrstfiles)) # color scaler
   for fname in unrstfiles:
      n += 1
      case = UT.basename(fname)
      p = get_parameters(case)
      well_k = _grid_numbers(case)[3]
      # sand drawdown and pressure support
      t = _find_drawdown_sand(case, 0, y_indx)[3]
      dpsand = []
      if plot_kwargs: kwargs = plot_kwargs(fname, unrstfiles)
      else          : kwargs = {'color':(red[n], 0, 1-red[n])}
      for tstep in range(len(t)):
         dpsand.append((_find_drawdown_sand(case, tstep, y_indx)[0])[-1, z_indx-1])
      plot(t[1:], dpsand[1:], label=case, **kwargs)  # skip first point which is 0 (only confusing)
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('dP sand. j,k = %i,%i' % (y_indx,z_indx))

def _dp_sand3(cases, y_indx, plot_kwargs=None):
   '''
   plots dp_sand for the given cases
   '''
   global y_
   figure()
   n = -1
   if not plot_kwargs: red = linspace(0., 1., len(cases)) # color scaler
   for case in cases:
      n += 1
      icd_segm  = _j2segm(case,y_indx,True)
      well_segm = _j2segm(case,y_indx,False)
      p = get_parameters(case)
      s = get_summary(case)
      p_support = s.get('BPR-%i,1,1' % p.nx)   # assumes this block is representative for the pressure support
      dp_sand = p_support - s.get('SPR-%s-%i'%(WELLNM,well_segm)) - s.get('SPRD-%s-%i'%(WELLNM,icd_segm))
      if plot_kwargs: kwargs = plot_kwargs(case, cases)
      else          : kwargs = {'color':(red[n], 0, 1-red[n])}
      plot(s.time[1:], dp_sand[1:], label=case, **kwargs)  # skip first point which is 0 (only confusing)
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('dP sand. j = %i' % (y_indx))
   y_ = dp_sand

def _dp_valve(case, y_indices):
   '''
   plots dp_valve for the given section
   '''
   figure()
   s = get_summary(case+'.RSM')
   n = -1
   for j in y_indices:
      n += 1
      plot(s.time, s.get_segm_data('SPRD', WELLNM, [_j2segm(case,j,True)]), COLOURS[n]+'-', label='j=%i'%j)
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('%s. dP inflow valve' % case)

def _dp_valve2(y_indx, cases, plot_kwargs=None):
   '''
   plots dp_valve for the given section
   '''
   figure()
   if not plot_kwargs: red = linspace(0., 1., len(cases)) # color scaler
   n = -1
   for case in cases:
      n += 1
      s = get_summary(case)
      if plot_kwargs: kwargs = plot_kwargs(case, cases)
      else          : kwargs = {'color':(red[n], 0, 1-red[n])}
      plot(s.time, s.get_segm_data('SPRD', WELLNM, [_j2segm(case,y_indx,True)]), label=case, **kwargs)
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('dP inflow valve. j=%i' % y_indx)

def _press_support(case, y_indices):
   '''
   plots pressure support for the given cell
   '''
   figure()
   n = -1
   for j in y_indices:
      n += 1
      t = _find_drawdown_sand(case, 0, j)[3]
      p_support = []
      for tstep in range(len(t)):
         p_support.append(CM2.get(case+'.UNRST').get('PRESSURE', tstep)[-1,j-1,0])
      plot(t, p_support, COLOURS[n]+'.-', label='j=%i'%j)
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('%s. Pressure support' % case)

def _press_support2(unrstfiles, y_indx):
   '''
   plots pressure support for the given cell
   '''
   figure()
   n = -1
   for fname in unrstfiles:
      n += 1
      case = UT.basename(fname)
      t = _find_drawdown_sand(case, 0, y_indx)[3]
      p_support = []
      for tstep in range(len(t)):
         p_support.append(CM2.get(case+'.UNRST').get('PRESSURE', tstep)[-1,y_indx-1,0])
      plot(t, p_support, COLOURS[n]+'.-', label=case)
   xlabel('time [days]')
   ylabel('[bar]')
   grid(True)
   legend(loc='best')
   title('Pressure support. j=%i'%y_indx)

def _valves_pr_joint(datafile):
   '''
   assumes RCP valve
   '''
   p = get_parameters(UT.basename(datafile))
   if p.inflow_control == 'AR3' : 
      valve_type = 'ar3_peregrino'
   else:
      raise Exception("No such inflow control: %s" % p.inflow_control)
   l = VC.rcp_length(p.dp_valve_oil, p.qliq*UI.Bo/p._nsegm, p._dy, p.rho_oil, p.mu_oil, valve_type)
   vpj = JOINT_LENGTH / l # valves pr joint
   return vpj

def _nozzle_diameter(datafile):
   diam = float(UT.grep_column(datafile, '-- valve has diam =', 6, False)[0].split('mm')[0])
   return diam

def _valves_pr_joint2(datafile):
   '''
   assumes RCP valve. reads vpj based from actual datatfile.
   does not handle defaulted values ('*') or ranges of segments.
   returns vpj for each segment.
   '''
   f = open(datafile)
   in_wsegaicdsection = False
   licds = []
   for line in f:
      if line.startswith('--') : continue  # skip commented line
      if line.strip() == ''    : continue  # skip blank line
      if 'WSEGAICD' in line:
         in_wsegaicdsection = True
         continue
      if not in_wsegaicdsection: continue
      if line.strip() == '/'   : break     # marks end of WSEGAICD-section
      licds.append(float(line.split()[4])) # so - let's read the data
      if line.count('/') == 2  : break     # marks end of WSEGAICD-section
   f.close()
   return JOINT_LENGTH/array(licds)

def _plot_valves_pr_joint(datafiles, plot_kwargs=None):
   global y_
   figure()
   if not plot_kwargs: red = linspace(0., 1., len(datafiles)) # color scaler
   n = -1
   for datafile in datafiles:
      n += 1
      if plot_kwargs: kwargs = plot_kwargs(datafile, datafiles)
      else          : kwargs = {'color':(red[n], 0, 1-red[n])}
      kwargs['marker'] = '*'
      y_ = _valves_pr_joint2(datafile)
      plot(y_, label=UT.basename(datafile), **kwargs)
   xlabel('segment # [-]')
   ylabel('valves pr. joint [-]')
   grid(True)
   legend(loc='best')
   title('valves pr. joint')

def _nvalves(cases):
   '''
   finds number of valves for the complete well.
   assumes constant segment length (=_dy), but vpj may vary per segment.
   '''
   y = []
   for case in cases:
      case = UT.basename(case)
      vpj = _valves_pr_joint2(case+'.DATA')
      njoints_pr_segm = _get_param_value(case, '_dy') / JOINT_LENGTH
      y.append(sum(vpj)*njoints_pr_segm)
   return y

def _legend(cases):
   leg = []
   for case in cases:
      s = get_summary(case)
      txt = ''
      vpj = _valves_pr_joint(s.nm+'.DATA')
      txt += '%.1f vpj' % mean(vpj)
      if vpj[0] != vpj[-1]: txt += ' *' # non equi-spaced
      leg.append(txt)
   legend(leg, loc='best')

def _accum_pr_section(varnm, cases, y_indices, use_icd_segm):
   figure()
   for case in cases:
      s = get_summary(case)
      segments = [_j2segm(s.nm, j, use_icd_segm) for j in y_indices]
      y = []
      t = []
      i = -1
      for t0 in s.time:
         i += 1
         t.append(t0)
         y.append(sum(s.get_segm_data(varnm, WELLNM, segments)[:, i]))
      Y = [trapz(y[:i], t[:i]) for i in arange(1,len(t))]
      plot(t[1:], Y, label=s.nm)
   xlabel('time [DAYS]')
   ylabel('SM3')
   title('%s on segments %s' % (varnm, segments))
   grid(True)
   legend(loc='best')


def _accum_pr_part(varnm, cases, t_end, pr_segment=False, compl_nm_only=False):
   figure()
   col = ''
   nfiles = len(cases)
   m = zeros((nfiles, 3))
   names = []
   n = -1
   for case in cases:
      n += 1
      s = get_summary(case)
      p = get_parameters(s.nm)
      # indices are a mess...
      nsegm = p._nsegm                                     # for convinience
      s1, s2 = [_j2segm(s.nm, p.channel_j1, True), _j2segm(s.nm, p.channel_j2, True)]
      segm_heel = arange(nsegm+2, s1)
      segm_chnl = arange(s1, s2+1)
      segm_toe  = arange(s2+1, 2*nsegm+1)
      y1 = []; y2 = []; y3 = []; t = []
      i = -1
      for t_ in s.time:
         if t_ > t_end: break
         i += 1
         t.append(t_)
         y1.append(sum(s.get_segm_data(varnm, WELLNM, segm_heel)[:, i]))
         y2.append(sum(s.get_segm_data(varnm, WELLNM, segm_chnl)[:, i]))
         y3.append(sum(s.get_segm_data(varnm, WELLNM, segm_toe) [:, i]))
      if pr_segment: (n1, n2, n3) = (len(segm_heel), len(segm_chnl), len(segm_toe)) # num segm pr part
      else         : (n1, n2, n3) = (1., 1., 1.)
      m[n,0] = (trapz(y1, t)/n1)  # heel section
      m[n,1] = (trapz(y2, t)/n2)  # channel section
      m[n,2] = (trapz(y3, t)/n3)  # toe section
      col += 'kbrm'
      if compl_nm_only: names.append(_compl_nm(UT.basename(s.nm)))
      else:             names.append(UT.basename(s.nm))
   ind = arange(nfiles)
   p1 = bar(ind, m[:,0], color='y')
   p2 = bar(ind, m[:,1], color='r', bottom=m[:,0])
   p3 = bar(ind, m[:,2], color='c', bottom=m[:,0]+m[:,1])
   xticks(arange(nfiles)+0.4, names)
   xlim(-0.5, nfiles+1.2)
   if pr_segment: xtra = '(Pr. segment)'
   else         : xtra = ''
   title('Accumulated %s @ %i days. %s' % (varnm, t_end, xtra))
   legend((p1[0],p2[0],p3[0]), ('heel', 'channel', 'toe'), loc='upper right')

def _weighted_flow(q_segm, params):
   w = ones(q_segm.shape)
   q_tot = sum(q_segm)
   for pm in params:
      w[pm[0]-1:pm[1]] = pm[2]*w[pm[0]-1:pm[1]] # fra og med den forste, til og med den siste
   q = w * q_segm
   q *= q_tot/sum(q)
   return q

def _perturb(K, shift, scaler_lo, scaler_hi):
   Knew = K.copy()
   if shift > 0:
      Knew[shift:] = K[:-shift]
      Knew[:shift] = K[0]
   elif shift < 0:
      Knew[:shift] = K[-shift:]
      Knew[shift:]  = K[-1]
   m = mean(Knew)
   ind = (Knew <= m).nonzero()[0]
   Knew[ind] = scaler_lo * Knew[ind]
   ind = (Knew >= m).nonzero()[0]
   Knew[ind] = scaler_hi * Knew[ind]
   return Knew

def _add_timeseries(varnm, cases):
   s0 = get_summary(cases.pop(0))
   ysum = 0
   for case in cases:
      s = get_summary(case)
      varnm = _fix(r,varnm)
      ysum += interp(s0.time, s.time, s.get(varnm))
   return (s0.time, (ysum/len(cases) - s0.get(varnm))/s0.get(varnm)*100.)

def _read_relperm(swof_file):
   sfile = open(swof_file)
   header = True
   sw = []; kr_w = []; kr_o = []
   for line in sfile:
      if 'SWOF' in line:
         header = False
         continue
      if header: continue
      if line.startswith('/') : break     # done
      if line.startswith('--'): continue  # skip comments
      rec = line.split()
      if len(rec) == 4:                   # allow blank lines
         sw.append(float(rec[0]))
         kr_w.append(float(rec[1]))
         kr_o.append(float(rec[2]))
   sfile.close()
   return (array(sw), array(kr_w), array(kr_o))

def plot_relperms(swof_files):
   figure()
   i = -1
   for swof_file in swof_files:
      i += 1
      sw, kr_w, kr_o = _read_relperm(swof_file)
      plot(sw[:-1], kr_w[:-1], COLOURS[i], label=UT.basename(swof_file))
      plot(sw[:-1], kr_o[:-1], COLOURS[i], label='_nolegend_')
   xlabel('water saturation [-]')
   ylabel('relperm [-]')
   legend()
   grid(True)

def _ali_permx():
   fname = '/project/multiscale/ICD_Rotvoll/Brazil/runs/GENERIC1/TEST_ATLE.permx'
   dx = 50.
   dy = 50.
   dz = 1.
   m = ECL.read_property(fname, 'PERMX')

def _plot_swat_contours(i1,i2, j, k, fnames):
   for fname in fnames:
      case = UT.basename(fname)
      p = get_parameters(case)
      r3 = CM2.add(fname,nx=p.nx, ny=p.ny, nz=p.nz)
      r3.contour_plot('SWAT', 'X', i1,i2, j,99999, k,99999)
      savefig('%s.png' % case)

def _get_box_value(r3, varnm, tstep, i1,i2, j1,j2, k1,k2, averaged):
   '''
   must take wedge shape into account when calculating mean value inside box
   '''
   m = r3.get(varnm, tstep)
   p = get_parameters(r3.casenm)
   dz_scaler = linspace(p.wedge_factor, 1., p.nx+1) # wedge weighting
   mm    = 0                                        # mean value to be calculated
   for i in arange(i1,i2):
      for j in arange(j1,j2):
         for k in arange(k1,k2):
            mm += m[i,j,k]
            if averaged: mm *= dz_scaler[i]
   if averaged:
      w_tot = (j2-j1)*(k2-k1) * sum(dz_scaler[i1:i2])  # tot weight
      return mm/w_tot
   else:
      return mm

def _plot_box_profile(varnm, i1,i2, j1,j2, k1,k2, fnames, averaged):
   figure()
   i = -1
   for fname in fnames:
      i += 1
      case = UT.basename(fname)
      p = get_parameters(case)
      r3 = CM2.add(fname,nx=int(p.nx), ny=int(p.ny), nz=int(p.nz))
      y = []
      t = [] # dont use s.time since some data is not for all timesteps
      for tstep in range(len(r3.data[varnm])):
         y.append(_get_box_value(r3, varnm, tstep, i1,i2, j1,j2, k1,k2, averaged))
         t.append(r3.time[tstep])
      plot(t, y, COLOURS[i], label=UT.basename(fname))
   xlabel('time [days]')
   ylabel(varnm)
   legend(loc='best')
   titl = '%s in box (%i-%i, %i-%i, %i-%i)' % (varnm, i1+1,i2, j1+1,j2, k1+1,k2)
   if averaged: titl = 'Avaraged ' + titl
   title(titl)
   grid(True)
   return (t, y)

def _plot_box_profile2(varnm, i1,i2, j_pairs, k1,k2, fname, averaged):
   figure()
   i = -1
   case = UT.basename(fname)
   p = get_parameters(case)
   r3 = CM2.add(fname,nx=p.nx, ny=p.ny, nz=p.nz)
   for j_pair in j_pairs:
      j1 = j_pair[0] - 1
      j2 = j_pair[1]
      i += 1
      y = []
      for tstep in range(r3.ntimesteps):
         y.append(_get_box_value(r3, varnm, tstep, i1,i2, j1,j2, k1,k2, averaged))
      plot(r3.time, y, COLOURS[i], label='j%i - j%i' % (j1+1,j2))
   xlabel('time [days]')
   ylabel(varnm)
   legend(loc='best')
   title('Average %s in box (i%i-i%i, k%i-k%i) for case %s' % (varnm, i1+1,i2, k1+1,k2, fname))
   grid(True)

def _fluids_in_place(prtfiles, fluid, plotit=True):
   '''
   for some reason, FWIP is incorrect for some cases (does not match FWPT)
   '''
   if plotit: figure()
   i = -1
   for prtfile in prtfiles:
      i += 1
      f = open(prtfile)
      y = []; t = []
      for line in f:
         if 'BALANCE  AT' in line: # time
            rec = line.split()
            t.append(float(rec[2]))
            continue
         if 'CURRENTLY IN PLACE' in line:
            if len(y) == len(t): continue   # redundant lines
            rec = line.split()
            if   fluid.lower() == 'oil'  : val = float(rec[4])
            elif fluid.lower() == 'water': val = float(rec[6])
            else                         : val = float(rec[4]) + float(rec[6]) # total
            y.append(val)
            continue
      f.close()
      if plotit: plot(t, y, COLOURS[i], label=UT.basename(prtfile))
   if plotit:
      xlabel('time [days]')
      ylabel('[m3]')
      legend(loc='best')
      title('Fluid In Place: ' + fluid)
      grid(True)
   return t, y

def _pvi(case, swi=0.05):
   '''
   find PVI - Pore Volume Injected.
   since we use MULTPV for the aquifier and dont actually inject anything,
   we calculate pore volume of the oil zone only, and use the produced liquid to
   find PVI.
   '''
   porevol = (1. + swi)*_fluids_in_place(('%s.PRT'%case,), 'OIL', False)[1][0]
   s = get_summary(case+'.RSM')
   return s.get('FLPT') / porevol

def _effectivedp_vs_flow(cases, j1, j2, pres=237., skip=7):
   '''
   plot effective dp vs flow for two different j-positions (valve characteristics including reservoir).
   pres: assumes constant pressure in water zone.
   skip: possible to skip first timesteps since the behavior the first few
   timesteps is strange (BHP decreases before increasing)
   '''
   figure()
   colors = 'krymcgb'
   i = -1
   for case in cases:
      i += 1
      s = get_summary(case)
      dp1 = pres - s.get('SPR-%s-%i' % (WELLNM, _j2segm(case, j1)))
      dp2 = pres - s.get('SPR-%s-%i' % (WELLNM, _j2segm(case, j2)))
      q1  = s.get('SLFR-%s-%i'%(WELLNM, _j2segm(case, j1, True)))[skip:]
      q2  = s.get('SLFR-%s-%i'%(WELLNM, _j2segm(case, j2, True)))[skip:]
      plot(q1, dp1[skip:], '%s*'%colors[i], label='%s j=%i'%(UT.basename(case), j1))
      plot(q2, dp2[skip:], '%so'%colors[i], label='%s j=%i'%(UT.basename(case), j2))
   grid(True)
   legend(loc='best')
   xlabel('flowrate [m3/d/segm]')
   ylabel('total drawdown [bar]')

def _cases(datafiles):
   cases = []
   for datafile in datafiles:
      #if os.path.exists(datafile): cases.append(UT.basename(datafile))
      if os.path.exists(datafile): cases.append(datafile)
   if not cases: raise Exception("WARNING: no DATA-files found")
   if UI.sortit: cases.sort()
   return cases

def _robustness(cases, varnm, ngroups=3):
   groups = _groups(cases, ngroups)
   figure()
   i = 0
   for group in groups:
      i +=1
      c = COLOURS[i]
      vals = get_val(varnm, group)
      mval = mean(vals)
      plot(i*ones(len(vals)), vals, c+'*', markersize=12)
      plot(i, mval, c+'o', markersize=20)
   grid(True)
   legend(loc='best')
   ylabel(get_summary(cases[0]).unit(varnm))
   xlim(0.5, ngroups+0.5)
   names = ['none', 'rcp', 'icd', 'vfp'][:ngroups]
   xticks(arange(ngroups)+1, names, fontsize=20)
   title(varnm)

def _valve_characterstics(casenm, segments, tsteps, plot_perf_curves=True, ax_limits=(0,700, -1, 75)):
   global y_
   p = get_parameters(casenm)
   segs = icd_segments(casenm, WELLNM)
   id = ECL.InflowData(casenm, 0, Inf, 12./p.valves_pr_joint, segs[0], segs[-1], 950, 1050, False)
   qmax = max(id.q[:,1:].flatten()) # avoid inf for first timestep
   qv   = linspace(0,qmax, 50)
   dr = 1/float(len(id.t))
   figs = []
   for segm in segments: figs.append(figure().number)
   #
   for tstep in arange(len(id.t))[tsteps]:
      red = dr*tstep
      n = -1
      for segm in segments:
         n += 1
         figure(figs[n])
         segm_ = segm - segs[0]
         q = id.q[segm_,tstep]
         if plot_perf_curves:
            if 'ar' in p.inflow_control.lower():
               dpv = VC.rcp_dp2(p.inflow_control, id.rho[segm_,tstep], id.mu[segm_,tstep], qv)
            else:
               # icd: dp = c*q**2
               c = id.dp[segm_,tstep] / q**2
               dpv = c*qv**2
            plot(qv, dpv, '--', color=(red,0,1-red),label='%i days'%id.t[tstep])
         plot(q, id.dp[segm_,tstep], 'o', color=(red,0,1-red),label='_nolabel_', markersize=12)
   n = -1
   for segm in segments:
      n += 1
      figure(figs[n])
      title('%s - segm %i' % (casenm, segm))
      xlabel('flow rate [l/h]')
      ylabel('dp [bar]')
      grid(True)
      axis(ax_limits)
      legend(loc='best')
   y_ = id

def _fluidheight_contours(case, fluid, tsteps, hmin, hmax, rel_h):
   '''
   rel_h: if True, will give relative height (useful for wedge shapes, for example)
   '''
   global y_
   nx, ny, nz = ECL.get_dimensions(case+'.DATA')
   r3 = CM2.add(case+'.UNRST', nx=nx, ny=ny, nz=nz)  # 3D data
   ri = CM3.add(case+'.INIT', nx=nx, ny=ny, nz=nz)   # 3D data
   x = cumsum(ri.get('DX',0)[:,0,0])
   y = cumsum(ri.get('DY',0)[0,:,0])
   fh = zeros((nx,ny)) # fluid height
   if not tsteps: tsteps = range(len(r3.time))
   for tstep in tsteps:
      sw = r3.get('SWAT',tstep)
      for i in range(nx):
         for j in range(ny):
            if fluid == 'oil': sat = 1. - sw[i,j,:]
            else             : sat = sw[i,j,:]
            fh[nx-i-1,j] = sum(sat*ri.get('DZ',0)[i,j,:])
            if rel_h:
               fh[nx-i-1,j] /= sum(ri.get('DZ',0)[i,j,:]) / 100.  # in %
      ax = PU.contourf(y,x, fh, zmin=hmin, zmax=hmax)
      ax.axis('equal')
      titl = '%s. %s height @ %i days' % (case, fluid, r3.time[tstep])
      if rel_h: titl += ' (relative height)'
      ax.set_title(titl)
      ax.set_xlabel('x-coord [m]')
      ax.set_ylabel('y-coord [m]')
      ax.set_yticks((0,max(x)))
      ax.set_yticklabels(('0', '%i' % (int((max(x)/100*1.1))*100) ) )
   show()
   y_ = fh

def help(mode=None):
   '''
   prints some useful info about the analyze-function.
   if 'mode' is given, only that functionality is explained
   '''
   import sys, inspect
   fname = inspect.getabsfile(sys.modules[__name__])  # *this* file
   f = open(fname)
   inside_analyze_func = False
   docs = {}
   while True:
      line = f.readline()
      if line.startswith('def analyze'):
         inside_analyze_func = True
         continue
      if not inside_analyze_func : continue
      if line.startswith('def')  : break    # we passed the analyze function
      line = line.rstrip()
      if 'mode ==' in line:
         rec = line.split()
         mode_ = rec[3].strip().replace("'","").replace(":","")
         docs[mode_] = []
         line = f.readline()
         # include some lines below the 'mode ==' line
         while '##' in line or 'args[' in line or 'enough_args' in line:
            if not 'enough_args' in line: docs[mode_].append(line)
            line = f.readline()
         continue
   f.close()
   if not mode:
      print "for analyze(mode, *args), the following modes are available:"
      keys = docs.keys(); keys.sort()
      for mode_ in keys: print mode_
   else:
      if not docs.has_key(mode):
         print 'No such mode: %s' % mode
         return
      for line in docs[mode]:
         print line.strip()

def _statistics(varnm, t0, plot_all, ref_value, descriptions, case_groups):
   if not len(descriptions) == len(case_groups):
      print "length of descriptions-list must match number of case_groups"
      return
   # handle ref_value:
   relative = False
   try:
      len(ref_value)
      relative = True
   except Exception:
      if ref_value > 0:
         relative = True
         ref_value = ref_value*ones(len(case_groups))
   figure()
   i = 0
   for case_group in case_groups:
      y = []
      for case in case_group:
         s = get_summary(case)
         if t0 < 0:
            ind = -1
         else:
            ind = UT.find_index(s.time, t0, max(1., t0/10.))
         y.append(s.get(varnm)[ind])
      y  = array(y)
      if relative: y = (y-ref_value[i])/ref_value[i] * 100 # want percent
      m  = mean(y)
      sd = std(y)
      plot(i, m, 'mo', markersize=15)
      plot((i,i), (m-sd, m+sd), 'k-', linewidth=2)
      plot(i, m+sd, 'k_', markersize=10)
      plot(i, m-sd, 'k_', markersize=10)
      if plot_all:
         plot(i*ones(len(y)), y, 'kx', markersize=6)
      print '%s: mean=%.1f stddev=%.1f' % (descriptions[i], m, sd)
      i += 1
   xticks(range(i), descriptions)
   if relative > 0:
      ylabel('rel. change in %s [%%]'%varnm)
   else:
      ylabel('%s [%s]'%(varnm, s.unit(varnm)))
   title(os.getcwd().split('/')[-1])
   xlim(-0.2, i-0.8)
   legend(('mean', 'stddev'), loc='best')

def _mean_dpvalve(cases, tstep):
   '''
   note: tstep = -1 gives mean dp for the entire simulation.
   '''
   dps = []
   for case in cases:
      s = get_summary(case)
      dp = []
      for segmno in s.icd_segments(WELLNM):
         dp_ = s.get('SPRD-%s-%i'%(WELLNM, segmno))[1:]
         dt = diff(s.time)
         dp.append(cumsum(dt*dp_)[tstep] / s.time[tstep])
      dps.append(mean(dp))
   return dps

def _plot_block_data_at_given_time(varnm, blocks, t0, cases):
   figure()
   for case in cases:
      s = get_summary(case)
      ind = UT.find_closest_index(s.time, t0)
      print 'using time', s.time[ind]
      plot(arange(len(blocks))+1, s.get_block_data(varnm, blocks, tsteps=[ind])[:,0], label=s.shortnm)
   grid(True)
   xlabel('blocks [-]')
   ylabel(varnm)
   legend(loc='best')
   title('%s @ %.1f days'%(varnm, t0))
   show()

def _plot_single_block_data(varnm, block, cases):
   figure()
   for case in cases:
      s = get_summary(case)
      plot(s.time, s.get_block_data(varnm, [block]), label=s.shortnm)
   grid(True)
   xlabel('TIME [days]')
   ylabel(varnm)
   title('Block = %i, %i, %i' % block)
   legend(loc='best')
   show()

def _var_vs_dpvalve(patterns, descriptions,
                    vpj_indxs=None, varnm='FOPT', tstep=-1, ref_case=None, plotit=True, colours=None, relative=False):
   dp = []
   y  = []
   n_groups = len(patterns)
   for pattern in patterns:
      cases = _cases(UT.glob(pattern))
      y.append(get_val(varnm, cases, tstep))
      dp.append(_mean_dpvalve(cases, tstep))
   if plotit:
      if not colours:
         if n_groups <= 4: colours = ['r', 'g', 'm', 'c']
         else            : colours = COLOURS
      c = colours # for convinience
      figure()
      if ref_case is not None: ref_value = get_val(varnm, ref_case, tstep)[0]
      if relative: scaler = 100/ref_value  # percent
      else       : scaler = 1.             # convinient
      for n in range(n_groups):
         plot(dp[n], array(y[n])*scaler, c[n]+'--s', label=descriptions[n])
      #if ref_case is not None:
      #   plot(0, ref_value*scaler, 'k*', markersize=18, label='Ref.')
      if vpj_indxs is not None:
         # mark 1 vpj and 2 vpj with triangles
         for n in range(n_groups):
            i1 = (array(vpj_indxs[n]) == 1).nonzero()[0][0]
            i2 = (array(vpj_indxs[n]) == 2).nonzero()[0][0]
            plot(dp[n][i1], y[n][i1], 'kv', markersize=8)
            plot(dp[n][i2], y[n][i2], 'k^', markersize=8)
      legend(loc='lower center')
      xlabel('Mean averaged dp-valve [bar]')
      if relative: ylabel('Relative change in %s [%%]' % varnm)
      else       : ylabel('%s [%s]' % (varnm, get_summary(cases[0]).unit(varnm)))
      grid(True)
      a = axis()
      axis((-5, 1.02*a[1],  0.98*a[2], 1.02*a[3]))
   return dp, y

def enough_args(args, n_args, mode):
   # go to help if number of  args is too low
   if not len(args) >= n_args:
      print 'number of arguments = %i < %i. showing help...' % (len(args), n_args)
      help(mode)
      return False
   else:
      return True

def analyze(mode, *args):
   global y_
   if not WELLNM:
      print 'You need to initialize some functions and variables (typically run PP_setup.py)'
      return
   if mode == 'run'  : 
      ## runs the given cases.
      ## example: analyze('run', 'A*.DATA')
      if not enough_args(args, 1, mode): return
      cases = UT.glob(args[0:])
      ECL.lsf_run2(cases)
   elif mode == 'check'  : 
      ## check if simulation had problems/warnings/errors
      ## example: analyze('check', 'A*.LOG')
      if not enough_args(args, 1, mode): return
      logfiles = UT.glob(args[0:])
      logfiles.sort()
      figure(1); clf()
      figure(2); clf()
      red = 0.
      for logfile in logfiles:
         os.system('grep -H Errors %s' % logfile)
         os.system('grep -H Problems %s' % logfile)
         os.system('grep -H Warnings %s' % logfile)
         os.system('grep elapsed %s' % logfile)
         figure(1)
         ECL.analyze_run((logfile,), newfig=False, use_basename=False, color=(red, 0, 1-red))
         figure(2)
         ECL.analyze_run((logfile,), newfig=False, use_basename=False, messagetype='WARNING', color=(red, 0, 1-red))
         red += 1. / len(logfiles)
   elif mode == 'get':
      ## just report last value of the given variable
      ## result is put into global y_ variable (and is also returned)
      ## example: analyze('get', 'FOPT', 'A*.DATA')
      if not enough_args(args, 2, mode): return
      varnm = args[0]
      cases = _cases(UT.glob(args[1:]))
      y_ = get_val(varnm, cases)
      return y_
   elif mode == 'get_par':
      ## just get the requeste parameter for given case(s)
      ## example: analyze('get_par', 'valve_pr_joint', 'A*.DATA')
      if not enough_args(args, 2, mode): return
      parnm = args[0]
      cases = _cases(UT.glob(args[1:]))
      y_ = _get_params(parnm, cases)
      return y_
   elif mode == 'xplot':
      ## cross plot any two variables against each other, for a group of cases.
      ## example: analyze('xplot', 'TIME', 'FOPT', 'A*.DATA')
      ## example: analyze('xplot', 'FWPT', 'FOPT', 'A*.DATA')
      if not enough_args(args, 3, mode): return
      varnm1 = args[0]
      varnm2 = args[1]
      cases  = _cases(UT.glob(args[2:]))
      _xplot(varnm1, varnm2, cases, plot_kwargs=UI.plot_kwargs)
      if UI.legend: UI.legend(cases)
   elif mode == 'meanplot':
      ## plot mean value of a given variable for a bunch of cases
      ## example: analyze('xplot', 'FWPT', 'FOPT', 'A*.DATA')
      if not enough_args(args, 1, mode): return
      varnm = args[0]
      cases = _cases(UT.glob(args[1:]))
      _meanplot(varnm, cases)
   elif mode == 'wxplot':
      ## cross plot any two variables against each other, for a group of wells for one case.
      ## example: analyze('wxplot', 'FWPT', 'WOPT', 'P*', 'A1.DATA')
      ## note that wellname-pattern ('P*') uses wildcard notation ala unix filenames.
      if not enough_args(args, 4, mode): return
      varnm1        = args[0]
      varnm2        = args[1]
      wname_pattern = args[2]
      case          = _cases(UT.glob(args[3]))[0]
      _wxplot(varnm1, varnm2, wname_pattern, case, plot_kwargs=UI.plot_kwargs)
      if UI.legend: UI.legend(cases)
   elif mode == 'comparison':
      ## compares a group of cases to a reference case. plots relative difference in %.
      ## example: analyze('comparison', 'FOPT', 'A1', 'A?.DATA')
      if not enough_args(args, 3, mode): return
      varnm   = args[0]
      refcase = args[1]
      cases   = _cases(UT.glob(args[2:]))
      _xplot('TIME', varnm, cases, refcase=refcase,plot_kwargs=None)
      if UI.legend: UI.legend(cases)
   elif mode == 'comparison2':
      ## compares a group of cases, two and two (meaning that case 1 is compared to 2, 3 to 4 and so on.
      ## plots relative difference in %.
      ## example: analyze('comparison2', 'FOPT', 'A?.DATA')
      if not enough_args(args, 2, mode): return
      varnm = args[0]
      cases = _cases(UT.glob(args[1:]))
      _comparison(varnm, cases, plot_kwargs=UI.plot_kwargs)
   elif mode == 'block':
      ## plot any variable for some given block (e.g. BPR)
      ## example: analyze('block', 'BPR', (10,18,3), 'A?.DATA')
      if not enough_args(args, 3, mode): return
      varnm          = args[0]          # typically BPR
      block          = args[1]          # (i,j,k)
      cases          = _cases(UT.glob(args[2:]))
      _plot_single_block_data(varnm, block, cases)
      if UI.legend: UI.legend(cases)
   elif mode == 'blocks':
      ## plot any variable for some given blocks (e.g. BPR)
      ## example: analyze('blocks', 'BPR', blocks, tstep, 'A?.DATA')
      if not enough_args(args, 4, mode): return
      varnm          = args[0]          # typically BPR
      blocks         = args[1]          # list of (i,j,k)'s
      t0             = args[2]          # at this given time [days]
      cases          = _cases(UT.glob(args[3:]))
      _plot_block_data_at_given_time(varnm, blocks, t0, cases)
      if UI.legend: UI.legend(cases)
   elif mode == 'segments':
      ## plot any variable along the well (e.g. SOFR) at a given time, for a group of cases.
      ## note: UI.plot_md is used.
      ## note: if UI.accum is True, it will accumulate in time, but this is only done
      ##       to first order accuracy (doesn't use trapz or similair techniques).
      ## example: analyze('segments', 'SOFR', 365, False, 'A?.DATA')
      if not enough_args(args, 5, mode): return
      varnm          = args[0]          # typically SOFR
      t0             = args[1]          # in days (-1 for last timestep does not work --> fix!)
      adjust_to_zero = args[2]          # boolean
      use_icd_segm   = args[3]          # boolean. if false, it uses well_segments
      accum          = args[4]          # boolean. if true, it will accumulate rates etc.
      cases          = _cases(UT.glob(args[5:]))
      _segment_data_at_given_time(
            varnm, cases, t0, plot_kwargs=UI.plot_kwargs, adjust_to_zero=adjust_to_zero, use_icd_segm=use_icd_segm, accum=accum)
      if UI.legend: UI.legend(cases)
   elif mode == 'segments2':
      ## plot any variable along the well (e.g. SOFR) at a given case, for different times.
      ## example: analyze('segments2', 'SOFR', 'A1.DATA', [365, 730, 1095], False)
      if not enough_args(args, 5, mode): return
      varnm          = args[0]           # typically SOFR
      case           = args[1]
      times          = args[2]          # time in days
      adjust_to_zero = args[3]          # boolean
      use_icd_segm   = args[4]          # boolean. if false, it uses well_segments
      _segment_data_for_given_case(varnm, case, times, adjust_to_zero=adjust_to_zero, use_icd_segm=use_icd_segm)
      if UI.legend: UI.legend(cases)
   elif mode == 'connections':
      ## plot any variable along the well (e.g. COFR) at a given time, for a group of cases.
      ## note: dx is used for the x-axis. assumed to be constant...
      ## note: if UI.accum is True, it will accumulate in time, but this is only done
      ##       to first order accuracy (doesn't use trapz or similair techniques).
      ## example: analyze('connections', 'COFR', 365, 50, 'A?.DATA')
      if not enough_args(args, 4, mode): return
      varnm          = args[0]          # typically COFR
      t0             = args[1]          # in days
      dx             = args[2]          # [m]
      cases          = _cases(UT.glob(args[3:]))
      _connection_data_at_given_time( varnm, cases, t0, dx, plot_kwargs=UI.plot_kwargs)
      if UI.legend: UI.legend(cases)
   elif mode == 'contours':
      ## makes the 'einstein-plots'.
      ## example: analyze('contours', 'SOFR', None, None, True, 'A?.DATA')
      if not enough_args(args, 5, mode): return
      varnm        = args[0]   # typically SOFR
      zmin         = args[1]   # if None, minimum value is used
      zmax         = args[2]   # if None, maximum value is used
      use_icd_segm = args[3]   # boolean. if false, it uses well_segments
      cases        = _cases(UT.glob(args[4:]))
      _contours(varnm, cases, UI.accum, zmin, zmax, use_icd_segm=use_icd_segm)
   elif mode == 'barplot':
      ## creates bar plot for any variable at given time, for a group of cases.
      ## example: analyze('barplot', 'FOPT', 3650, True, True, 'A?.DATA')
      if not enough_args(args, 5, mode): return
      varnm        = args[0]         # typically FOPT
      t0           = float(args[1])
      rel_height   = args[2]         # boolean
      mark_highest = args[3]         # boolean
      cases        = _cases(UT.glob(args[4:]))
      _barplot_vals(varnm, cases, t0, mark_highest=mark_highest, rel_height=rel_height)
   elif mode == 'wbarplot':
      ## creates bar plot for any variable at given time, for a group of wells (for one case)
      ## example: analyze('wbarplot', 'FOPT', 3650, True, True, '*', 'A.DATA')
      if not enough_args(args, 5, mode): return
      varnm         = args[0]         # typically FOPT
      t0            = float(args[1])
      rel_height    = args[2]         # boolean
      mark_highest  = args[3]         # boolean
      wname_pattern = args[4]
      case          = _cases(UT.glob(args[5]))[0]
      _wbarplot(varnm, wname_pattern, case, t0, 0, mark_highest, rel_height)
   elif mode == 'flatness':
      ## not much in use...
      if not enough_args(args, 3, mode): return
      varnm = args[0]   # typically SOFR
      tstep = args[1]
      cases = _cases(UT.glob(args[2:]))
      _standard_deviation_along_well(varnm, WELLNM, cases, tstep-1)
      if UI.legend: UI.legend(cases)
   elif mode == 'comparewells':
      if not enough_args(args, 4, mode): return
      varnm = args[0]   # typically FOPT
      well1 = args[1]
      well2 = args[2]
      cases = _cases(UT.glob(args[3:]))
      _comparewells(varnm, cases, well1, well2)
      titl = ''
      for str in args[4:]: titl += str + ' '
      title(titl)
   elif mode == 'drawdownc':
      ## not so useful, use drawdownp in stead
      if not enough_args(args, 3, mode): return
      case  = args[0]
      tstep = args[1]
      segno = args[2]
      _drawdown_contours(case, tstep, segno)
   elif mode == 'drawdownp':
      ## not much in use...
      if not enough_args(args, 3, mode): return
      case = args[0]
      y_indx = args[1]
      tsteps = args[2]
      _drawdown_profiles(case, tsteps, y_indx, False)
   elif mode == 'drawdown':
      ## this one gives a lot of info, showing both dp_sand and dp_valve + press_support
      if not enough_args(args, 2, mode): return
      case      = args[0]
      y_indices = args[1]
      _drawdown(case, y_indices)
   elif mode == 'press_support':
      ## example: analyze('press_support', 'A1', [20, 35])
      if not enough_args(args, 2, mode): return
      case      = args[0]
      y_indices = args[1]
      _press_support(case, y_indices)
   elif mode == 'press_support2':
      ## example: analyze('press_support2', 20, 'A?.UNRST')
      if not enough_args(args, 2, mode): return
      y_indx     = args[0]
      unrstfiles = UT.glob(args[1:])
      if UI.sortit: unrstfiles.sort()
      _press_support2(unrstfiles, y_indx)
   elif mode == 'dp_sand':
      ## example: analyze('dp_sand', 'A1', [20, 35], 1)
      if not enough_args(args, 3, mode): return
      case      = args[0]
      y_indices = args[1]
      z_indx     = int(args[2])
      _dp_sand(case, y_indices, z_indx)
   elif mode == 'dp_sand2':
      ## example: analyze('dp_sand2', 20, 1, 'A?.UNRST')
      if not enough_args(args, 3, mode): return
      y_indx     = int(args[0])
      z_indx     = int(args[1])
      unrstfiles = UT.glob(args[2:])
      if UI.sortit: unrstfiles.sort()
      _dp_sand2(unrstfiles, y_indx, z_indx, UI.plot_kwargs)
      if UI.legend: UI.legend(unrstfiles)
   elif mode == 'dp_sand3':
      ## doesnt use UNRST-files. assumes constant pressure support
      ## example: analyze('dp_sand3', 20, 'A?.DATA')
      if not enough_args(args, 2, mode): return
      y_indx = int(args[0])
      cases  = _cases(UT.glob(args[1:]))
      _dp_sand3(cases, y_indx, UI.plot_kwargs)
      if UI.legend: UI.legend(cases)
   elif mode == 'dp_valve':
      ## example: analyze('dp_valve', 'A1', [20, 35])
      if not enough_args(args, 2, mode): return
      case      = args[0]
      y_indices = args[1]
      _dp_valve(case, y_indices)
   elif mode == 'dp_valve2':
      ## example: analyze('dp_valve2', 20, 'A?.DATA')
      if not enough_args(args, 2, mode): return
      y_indx = int(args[0])
      cases  = _cases(UT.glob(args[1:]))
      _dp_valve2(y_indx, cases, UI.plot_kwargs)
      if UI.legend: UI.legend(cases)
   elif mode == 'xsection':
      ## shows a contour map of one layer of cells
      ## example: analyze('xsection', 'A1.DATA', 3, 'SWAT')
      if not enough_args(args, 3, mode): return
      case   = args[0]
      z_indx = args[1] - 1
      varnm  = args[2]             # typically SWAT
      _show_xsection(case, z_indx, varnm)
   elif mode == 'section':
      ## not much in use...
      if not enough_args(args, 3, mode): return
      varnm     = args[0]   # typically SOFR
      y_indices = args[1]
      cases     = _cases(UT.glob(args[2:]))
      _accum_pr_section(varnm, cases, y_indices)
   elif mode == 'vpj':
      ## plot valves pr joint for each case.
      ## example: analyze('vpj', 'A?.DATA')
      if not enough_args(args, 1, mode): return
      datafiles = UT.glob(args[0:])
      if UI.sortit: datafiles.sort()
      _plot_valves_pr_joint(datafiles, plot_kwargs=UI.plot_kwargs)
      if UI.legend: UI.legend(datafiles)
   elif mode == 'dp_well':
      ## plots pressure drop along the well (tubing)
      ## example: analyze('dp_well', 'A?.DATA')
      if not enough_args(args, 1, mode): return
      cases = _cases(UT.glob(args[0:]))
      _welldp(cases, plot_kwargs=UI.plot_kwargs)
   elif mode == 'swat_contours': # typically in the channel
      ## not much in use...
      if not enough_args(args, 5, mode): return
      i1         = args[0] - 1
      i2         = args[1]
      j          = args[2] - 1
      k          = args[3] - 1
      unrstfiles = UT.glob(args[4:])
      if UI.sortit: unrstfiles.sort()
      _plot_swat_contours(i1,i2, j, k, unrstfiles)
   elif mode == 'box_profile': # typically in the channel
      ## not much in use...
      if not enough_args(args, 9, mode): return
      varnm      = args[0]        # typically SWAT
      i1         = args[1] - 1
      i2         = args[2]
      j1         = args[3] - 1
      j2         = args[4]
      k1         = args[5] - 1
      k2         = args[6]
      averaged   = args[7]        # boolean
      unrstfiles = UT.glob(args[8:])
      if UI.sortit: unrstfiles.sort()
      _plot_box_profile(varnm, i1,i2, j1,j2, k1,k2, unrstfiles, averaged)
   elif mode == 'box_profile2':   # typically in the channel
      ## not much in use...
      if not enough_args(args, 8, mode): return
      varnm     = args[0]      # typically SWAT
      i1        = args[1] - 1
      i2        = args[2]
      j_pairs   = args[3]      # typically [(14,14), (24,24)]
      k1        = args[4] - 1
      k2        = args[5]
      unrstfile = args[6]
      averaged  = args[7]      # boolean
      _plot_box_profile2(varnm, i1,i2, j_pairs, k1,k2, unrstfile, averaged)
   elif mode == 'relperm':
      ## not much in use...
      if not enough_args(args, 1, mode): return
      swof_files = UT.glob(args[0:])
      plot_relperms(swof_files)
   elif mode == 'fluid_in_place':
      ## not much in use...
      if not enough_args(args, 2, mode): return
      fluid    = args[0]
      prtfiles = UT.glob(args[1:])
      if UI.sortit: prtfiles.sort()
      _fluids_in_place(prtfiles, fluid)
   elif mode == 'prod_pr_part':
      ## sees the reservoir as 3 different parts and plots
      ## the production for each part as a function of time
      if not enough_args(args, 4, mode): return
      varnm    = args[0]   # typically SOFR
      t_end       = int(args[1])
      pr_segment  = int(args[2]) # boolean
      cases = _cases(UT.glob(args[3:]))
      _accum_pr_part(varnm, cases, t_end, pr_segment)
   elif mode == 'prod_pr_part2':
      ## sees the reservoir as 3 different parts and plots
      ## the production for each case as barplots
      if not enough_args(args, 3, mode): return
      varnm   = args[0] # typically SOFR
      case = args[1]
      pr_segment  = int(args[2]) # boolean
      _prod_pr_part(varnm, case, pr_segment, UI.accum)
   elif mode == 'wbt':
      ## plots water break through time at various positions for
      ## different cases a function of a given paramater
      ## if parnm = 'dpvalve', it will plot it against averaged and mean valve dp
      ## example: analyze('wbt', [20, 35], 'dpvalve', 'A?.DATA')
      if not enough_args(args, 3, mode): return
      y_indices = args[0]
      parnm     = args[1]                  # typically valves_pr_joint
      relative  = args[2]                  # bool. relative to first index or not
      cases     = _cases(UT.glob(args[3:]))
      _wbt(cases, parnm, y_indices, relative)
   elif mode == 'var_vs_par':
      ## plot a (result) variable vs. a input parameter
      ## example: analyze('var_vs_par', 'FOPT', 'valves_pr_joint', True, 'A?.DATA')
      if not enough_args(args, 4, mode): return
      varnm  = args[0]                   # typically FOPT
      parnm  = args[1]                   # typically dp_valve_oil
      rel_ch = args[2]                   # boolean
      cases  = _cases(UT.glob(args[3:]))
      _variable_vs_param(varnm, parnm, cases, rel_ch)
   elif mode == 'var_vs_var':
      ## plot a (result) variable vs. another variable.
      ## example: analyze('var_vs_var', 'FWPT', 'FOPT', 3650 , 'A?.DATA')
      if not enough_args(args, 4, mode): return
      varnm1 = args[0]                   # typically FWPT
      varnm2 = args[1]                   # typically FOPT
      t0     = args[2]                   # at given time (in days)
      cases  = _cases(UT.glob(args[3:]))
      _variable_vs_variable(varnm1, varnm2, t0, cases)
   elif mode == 'var_vs_dpvalve':
      ## plot a (result) variable vs. averaged mean dp-valve
      ## very useful when comparing f.ex. AICD to ICD with varying number of vpj
      ## example: analyze('var_vs_dpvalve', ('A2?.DATA', 'B2?.DATA'), ['AICD', 'AICV'], vpjs, 'FOPT', -1, 137500, None)
      if not enough_args(args, 7, mode): return
      patterns     = args[0]             # typically ('A2?.DATA', 'A3?.DATA')
      descriptions = args[1]             # typically ('AICD', 'ICD')
      vpj_indxs    = args[2]             # vpjs-list (or None - if you dont need 1vpj and 2vpj to be marked)
      varnm        = args[3]             # typically FOPT
      tstep        = args[4]             # typically -1 (last time step)
      ref_case     = args[5]             # typically case with no inflow control (or None)
      colours      = args[6]             # typically None or ('g', 'c') etc.    (or None)
      relative     = args[7]             # boolean. print relative change compared to ref_case?
      y_ = _var_vs_dpvalve(patterns, descriptions, vpj_indxs, varnm, tstep, ref_case, True, colours, relative)
   elif mode == 'segm_diff':
      ## plot differences in variable-value between two different segments
      ## example: analyze('segm_diff', 'SOFR', 58, 75, 'A?.DATA')
      if not enough_args(args, 3, mode): return
      varnm = args[0]                    # typically SOFR
      segm1 = args[1]
      segm2 = args[2]
      cases = _cases(UT.glob(args[3:]))
      _segments_diff(cases, varnm, segm1, segm2, plot_kwargs=UI.plot_kwargs)
      if UI.legend: UI.legend(cases)
   elif mode == 'dp_vs_flow':
      ## plot effective dp vs flow
      ## example: analyze('dp_vs_flow', 20, 35, 'A?.DATA')
      if not enough_args(args, 3, mode): return
      j1 = args[0]
      j2 = args[1]
      cases = _cases(UT.glob(args[2:]))
      _effectivedp_vs_flow(cases, j1, j2)
   elif mode == 'compiled':
      ## creates a "1-pager" with a lot of plots
      ## example: analyze('compiled', 58, 75 'A?.DATA')
      ## note: must make sure UI.COMP_VARS is set first
      if not enough_args(args, 3, mode): return
      segm1  = args[0]
      segm2  = args[1]
      cases  = _cases(UT.glob(args[2:]))
      _compiled(cases, segm1, segm2, plot_kwargs=UI.plot_kwargs)
   elif mode == 'robustness':
      ## plots a lot of simulations in a graph as requested by sigurd,
      ## i.e. it plots the mean value and and all values for each group.
      ## to be used when a lot of different paramateres have been changed
      ## to check the robustness for each valve type, typically.
      ## example: analyze('robustness', 'FOPT', 3, 'A?.DATA')
      if not enough_args(args, 3, mode): return
      varnm   = args[0]                   # typically FOPT
      ngroups = args[1]                   # typically 3 (none, icd, rcp)
      cases   = _cases(UT.glob(args[2:]))
      _robustness(cases, varnm, ngroups)
   elif mode == 'valvechar':
      ## plots valve characteristics from the simulation.
      ## must have well info (using WELDEBUG) in DBG-file.
      ## example: analyze('valvechar', 'A1.DATA', [58, 75], [1, 10, 20])
      if not enough_args(args, 3, mode): return
      casenm   = args[0]
      segments = args[1]          # typically one inside channel and one outside: (57, 75)
      tsteps   = args[2]
      _valve_characterstics(casenm, segments, tsteps)
   elif mode == 'fluidheight':
      ## plots a contour-map showing 'effective' fluid height
      ## example: analyze('fluidheight', 'A1', 'oil', [1, 20], None, None, True)
      if not enough_args(args, 6, mode): return
      casenm = args[0]
      fluid  = args[1]   # oil or water
      tsteps = args[2]
      hmin   = args[3]   # if None, minimum value is used
      hmax   = args[4]   # if None, maximum value is used
      rel_h  = args[5]   # if True, use relative height (in %)
      _fluidheight_contours(casenm, fluid, tsteps, hmin, hmax, rel_h)
   elif mode == 'statistics':
      ## plots some statistics for the given cases.
      ## example: analyze('statistics', 'FWPT', -1, True, -1, ['OH', 'TR7', 'EQS'], 'B1?', 'B2?', 'B3?')
      if not enough_args(args, 6, mode): return
      varnm     = args[0]    # typically FOPT
      t0        = args[1]    # @ time, in days. if < 0, last tstep is used
      plot_all  = args[2]    # if true, all results will show as small dots (together with mean etc)
      ref_value = args[3]    # if > 0, all values will be relative (in percent) to this. could be vector or scalar.
      descr     = args[4]    # description of each case-group
      case_groups = []
      for pattern in args[5:]:
         case_groups.append(_cases(UT.glob(pattern)))
      _statistics(varnm, t0, plot_all, ref_value, descr, case_groups)
   else:
      print "no such mode: " + mode
   if UI.show_fig: show()

def xplots(varnms, *args):
   cases = UT.glob(args[0:])
   for varnm in varnms:
      _xplot('TIME', varnm, cases)

def contours(varnms, zmin, zmax, *args):
   cases = UT.glob(args[0:])
   for varnm in varnms:
      _contours(varnm, cases, UI.accum, zmin, zmax)

def create_tanklevel(fluid, rel_height, bg_pic, prefix='./'):
   '''
   creates a picture of a tank with some fluid in it.
   fluid is 'oil' or 'water'
   rel_height is in %
   bg_pic: the background tank picture. there is one on /project/RCP/active/fluent/Atle_Resources/Pics/tank.png
   '''
   height = (515-40)*rel_height/100.
   ystart = 515 - height
   if fluid == 'oil': colour = 'black'
   else             : colour = 'blue'
   fluid_pic = '_fluid_.png'
   cmd = 'convert -size 300x630 xc:white -fill %s -stroke black -draw "rectangle 25,%i  275,515" %s' \
      % (colour, ystart, fluid_pic)
   #print cmd
   os.system(cmd)
   fname = '%s%s_%05i.png' % (prefix, fluid, rel_height*100)
   os.system('composite -blend 80%%X20%% %s %s %s' % (bg_pic, fluid_pic,fname) )
   os.remove(fluid_pic)
   return fname

def create_setup_file(wellnm, fname='PP_setup.py'):
   f = open(fname, 'w')
   txt = '''
PP.icd_segments   = PP._icd_segments
PP.well_segments  = PP._well_segments
PP.WELLNM         = '%s'
PP.UI.plot_md     = True
PP.UI.plot_kwargs = None
PP.UI.legend      = None
''' % wellnm
   f.write(txt)
   f.close()
   print '%s has been created' % fname

def create_files(templfile, keyw, vals, prefix=None, fmt='%.2f', overwrite=False, file_ext='DATA'):
   '''
   easy way of creating a bunch of files with just one value that is changing.
   the tempfile is typically an input file (Eclipse: Data-file) with one value
   replaced by a special keyword (e.g. _LICD_).
   if prefix is not provided, it will use the basename of the templfile (upper case).
   '''
   if prefix is None:
      prefix = UT.basename(templfile).upper()
   ndigits = int(log10(len(vals))) + 1
   fnm_fmt = '%s%%0%ii.%s' % (prefix, ndigits, file_ext)
   i = 0
   fnms = []
   for val in vals:
      i += 1
      fnm = fnm_fmt % i
      if os.path.exists(fnm) and not overwrite:
         print fnm, 'exists. skipped'
         continue
      repl = [(keyw, fmt%val)]
      UT.replace_in_file(repl, templfile, fnm)
      print 'created file', fnm
      fnms.append(fnm)
   return fnms

def vpj2ratio(vpj):
   '''
   A function most useful for Troll... They use ratios like 1:1 instead of vpj=0.5
   n1:n2 means n1 screens followed by n2 blanks
   Note: only valid for vpj <= 1
   '''
   if vpj > 1: raise Exception('cannot handle vpj > 1')
   n1, a2 = float.as_integer_ratio(vpj)
   n2 = a2 - n1
   return '%i:%i'%(n1,n2), n1, n2



