import sys, re, glob, os, tempfile, pdb, time
import subprocess
import pylab as pl
import math
import shutil
import AtlejgTools.Utils as UT
import AtlejgTools.SimulationTools.UnitConversion as U

def _acosd(x):
   return math.degrees(math.acos(x))

def read_exported_file(fname, varnms):
   '''
   reads a file created by 'file export ...' in fluent.
   returns a hash with the requested data and each one of varnms as keys.
   varnms is a list of fluid varible names (like cell-volume)

   .cdat files: assumes the requested varnames comes first in the file. then it breaks.

   '''
   d = {}
   ################
   # CDAT files
   ################
   if '.cdat' in fname:
      # CFD-Post compatible file format
      f = open(fname)
      while True:
         line = f.readline()
         if not line: break
         if 'domain' in line:
            rec = line.split()
            varnm = rec[1][1:-1]
            print varnm
            if not varnm in varnms: break   # often the cdat-files has more variables than we asked for
            ncells = int(rec[7])
            print UT.basename(fname), varnm, ncells
            y = pl.zeros(ncells)
            # special treatment of first value
            f.readline()
            line = f.readline()
            y[0] = float(line.split()[1])
            # now read the values
            for i in pl.arange(1, ncells): y[i] = float(f.readline())
            # if case has more than 1 fluid zone, cdat-file will report them separately.
            # i will concatenate them.
            if d.has_key(varnm):
               d[varnm] = pl.concatenate((d[varnm], y))
            else:
               d[varnm] = y
      f.close()
   ################
   #
   ################
   return d

def read_force_report(filenm) :
   '''read summary of a fluent force report.
   not thoroughly tested...'''
   file = open(filenm , "r")

   for line in file:
      if 'net' in line:
         rec = line.split()
         f = {}
         f['pres'] = float(rec[1])
         f['visc'] = float(rec[2])
         f['tot']  = float(rec[3])
         file.close()
         return f

def read_report(filenm, zn_name) :
   '''"greps" out values that have been written for specified
   zone to a fluent report file
   not thoroughly tested...'''
   file = open(filenm , "r")
   # loop file
   vals = []
   for line in file:
      # split line and search for given zn_name
      rec = line.split()
      if len(rec) == 2 and zn_name == rec[0]:
         vals.append(float(rec[1]))
   file.close()
   return pl.array(vals)

def read_xyfile(fname, skip=0, last_entry_only=False, use_pylab=True):
   '''
   note: if multiple zones are reported in the file, then use_pylab should be False
   '''
   if last_entry_only: use_pylab = False
   if use_pylab:
      return pl.loadtxt(fname, skiprows=(4+skip), comments=')')
   file = open(fname , "r")
   headerlines = 4
   x = []; y = []
   i = 0
   for line in file:
      i += 1
      if i <= headerlines + skip                      : continue
      if len(line) <= 1                               : continue
      if line.startswith(')') or line.startswith('(') : continue
      if not last_entry_only :
         rec = line.split()
         x.append(float(rec[0]))
         y.append(float(rec[1]))
   if last_entry_only :
      rec = line.split()
      x = float(rec[0])
      y = float(rec[1])
   return (pl.array(x),pl.array(y))

def read_monitorfile(fname, skip=0, use_pylab=True):
   '''
   about use_pylab: see read_xyfile
   '''
   return read_xyfile(fname, skip=skip - 2, use_pylab=use_pylab)

def plot_monitorfiles(fnames, skip=0, xscaler=1, yscaler=1, use_pylab=True) :
   '''
   about use_pylab: see read_xyfile
   '''
   plot_xyfiles(fnames, skip - 2, use_pylab=use_pylab, yscaler=yscaler, xscaler=xscaler)

def plot_xyfiles(fnames, skip=0, xscaler=1, yscaler=1, use_pylab=True, vary_style=True):
   '''
   about use_pylab: see read_xyfile
   vary_style: if many files, try to make each plot-style unique
   '''
   pl.figure()
   if vary_style: ps = UT.PlotStyle()
   for fname in fnames:
      d = read_xyfile(fname, skip, use_pylab=use_pylab)
      ps.next()
      if vary_style:
         pl.plot(d[:,0]*xscaler, d[:,1]*yscaler, color=ps.colour(), marker=ps.marker(), label=fname)
      else:
         pl.plot(d[:,0]*xscaler, d[:,1]*yscaler, label=fname)
   pl.legend(loc='best')
   pl.gca().grid()
   #pl.show() # probably need this for windows

def read_dpm_report_file(fname, skip=0):
   '''
   reads a fluent dpm report file, returns 'matrix'.
   it reads the file, cleans it, and write numbers to a tempfile, which then is read by loadtxt
   '''
   file    = open(fname , "r")
   tmpfile = open(tempfile.mktemp(), "w")
   headerlines = 2
   i = 0
   for line in file:
      i += 1
      if i <= headerlines + skip : continue
      tmpfile.write(line[2:line.find(')')] + "\n") # trim line, get rid of paranthesis
   file.close()
   tmpfile.close()
   d = pl.loadtxt(tmpfile.name)
   os.unlink(tmpfile.name)
   return d
def read_flowtimes(fname):
   '''Reads flowtimes from transcript file'''
   trafile = open(fname)
   flowtimes = []
   for line in trafile:
      if not 'Flow time =' in line: continue
      rec = line.split()
      ftime = rec[3]
      flowtimes.append(float(ftime[:-2]))
   return pl.array(flowtimes)

def get_timesteps(t1, t2):
   '''
   given two simulations with diferent timestepping,
   this gives a map to find timesteps where simulation time are (more or less) the same
   typically, t1 is courser than t2
   note: note very thoroghly tested
   ex: t1 = FL.read_flowtimes('27.tra')
       t2 = FL.read_flowtimes('20.tra')
       ind = FL.get_timesteps(t1, t2)
       then t1[23] == t2[ind[23]]
   '''
   ind = pl.interp(t1, t2, range(len(t2)))
   tsteps = [int(i) for i in ind]
   return pl.array(tsteps)

def create_pipeline_geometry(coord, diams, journalfile, overlap=None):
   '''
   creates gambit journal-file which makes a pipeline geometry.
   coord      : the x- and y-coordinates of the junctions between segments.
   diams      : diameter of sections
   journalfile: the gambit journal file that will be created
   note1: journal starts with a reset - which is necesarry.
   note2: an extra segments is added (dont remember why - and do not dare to
          remove it...)
   NOTE!  consider using create_pipeline_geometry2 in stead...
   '''
   if not overlap: overlap = max(diams)
   # note: we are adding a horizontal first section
   # so that it is easy to rotate
   coord = pl.concatenate(([[coord[0,0] - 2*overlap, coord[0,1]]], coord))
   npoints   = len(coord)                       # adding section
   nsections = npoints - 1
   #
   # finding pipe-sections, lengths, and radius
   ps = pl.zeros((npoints, 2))                     # pipe-sections (incremental)
   l  = pl.zeros(nsections)
   r  = pl.zeros(nsections)
   for i in range(nsections):
      ps[i,:] = coord[i+1,:] - coord[i,:]
      l[i] = pl.norm(ps[i,:])
      if i > 0: r[i] = 0.5*diams[i-1]           # diam -> radius
   r[0] = r[1]                                  # added section
   #
   # finding (relative) angles between pipe-sections
   da = pl.zeros(nsections-1)
   for i in range(nsections-1):
      da[i] = 0.5*_acosd(pl.dot(-ps[i,:],ps[i+1,:])/(l[i]*l[i+1]))
      # need cross product to decide whether we go up or down from here
      v1 = pl.concatenate((ps[i], [0]))
      v2 = pl.concatenate((ps[i+1], [0]))
      c = pl.cross(v1,v2)
      if c[2] >= 0: da[i] = 90 - da[i]
      else:         da[i] = da[i] - 90
   #
   f = open(journalfile,'w')
   #
   f.write('reset\n')
   f.write('face create "for_split" width 1 height 1 yzplane rectangle\n')
   #
   # first two pieces must be handled separetely
   f.write('volume create height %.3f radius1 %.3f radius3 %.3f offset %.3f 0 0 xaxis frustum\n' % (l[0]+overlap,r[0],r[0],(l[0]+overlap)/2.))
   f.write('/\n')
   f.write('volume move "volume.1" offset %.3f 0 0 connected\n' % -l[0])
   f.write('/ #1\n')
   #
   f.write('volume create height %.3f radius1 %.3f radius3 %.3f offset %.3f 0 0 xaxis frustum\n' % (l[1]+2*overlap,r[1],r[1],(l[1]+2*overlap)/2.))
   f.write('volume move "volume.2" offset %.3f 0 0 connected\n' % -overlap)
   f.write('volume move "volume.2" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n' % -da[0])
   f.write('volume move "volume.1" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n' % da[0])
   f.write('volume split "volume.1" faces "for_split" connected keeptool\n')
   f.write('volume split "volume.2" faces "for_split" connected keeptool\n')
   f.write('volume delete "volume.2" "volume.3" lowertopology\n')
   f.write('volume split "volume.4" volumes "volume.1" connected bientity\n')
   f.write('volume move "volume.4" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n' % da[0])
   f.write('volume move "volume.4" offset %.3f 0 0 connected\n' % -l[1])
   f.write('/ #2\n')
   # this one is needed for final rotation
   f.write('edge split "edge.19" percentarclength 0.5 connected\n')
   #
   # then we can generalize
   for i in range(2, nsections):
      f.write('volume create height %.3f radius1 %.3f radius3 %.3f offset %.3f 0 0 xaxis frustum\n' % (l[i]+2*overlap,r[i],r[i],(l[i]+2*overlap)/2.))
      f.write('volume move "volume.%d" offset %.3f 0 0 connected\n' % (4*(i-1)+2,-overlap))
      if abs(da[i-1]) > 1e-4: # gets error message if angle is 0
         f.write('volume move "volume.%d" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n' % (4*(i-1)+2,-da[i-1]))
         f.write('volume move "volume.%d" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n' % (4*(i-1),da[i-1]))
      f.write('volume split "volume.%d" faces "for_split" connected keeptool\n' % (4*(i-1)))
      f.write('volume split "volume.%d" faces "for_split" connected keeptool\n' % (4*(i-1)+2))
      f.write('volume delete "volume.%d" "volume.%d" lowertopology\n' % (4*(i-1)+3,4*(i-1)+2))
      f.write('volume split "volume.%d" volumes "volume.%d" connected bientity\n' % (4*(i-1)+4,4*(i-1)))
      f.write('volume move "volume.%d" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n' % (4*(i-1)+4,da[i-1]))
      f.write('volume move "volume.%d" offset %.3f 0 0 connected\n' % (4*(i-1)+4,-l[i]))
      f.write('/i= %i\n' % i)
   #
   ## rotate and move to correct position
   s = '''
vertex create "origio" coordinates 0 0 0
volume align "volume.1" translation "vertex.5" "origio" connected
vertex create "x1" coordinates 1 0 0
vertex create "x2" coordinates 0 0 1
volume align "volume.1" translation "origio" "origio" rotation1 "vertex.21" \
  "x1" connected
'''
   f.write(s)
   f.close()
   #
   return ps, da, coord

def create_pipeline_geometry2(depths, length, diam, journalfile, overlap=None):
   '''
   creates gambit journal-file which makes a pipeline geometry.
   depths : depths of pipe-junctions
   lenght : length of pipe sections (or really: lateral distance between
            junctions)
   diams      : diameter of sections
   journalfile: the gambit journal file that will be created
   note1: journal starts with a reset - which is necesarry.
   '''
   coord = pl.zeros((len(depths),2))
   for i in range(len(depths)):
      coord[i,0] = i*length
      coord[i,1] = depths[i]
   diams = diam*pl.ones(len(depths))
   return create_pipeline_geometry(coord, diams, journalfile, overlap=overlap)

def read_convergence(trafile):
   f = open(trafile)
   iterations = []
   t          = []
   iter       = []
   first      = pl.Inf
   rightplace = False
   i = 0
   for line in f:
      i += 1
      if line.startswith('Flow time ='):
         t.append(float(line.split()[3].replace('s,','')))
         iter.append(curr-first)
         first = pl.Inf
         rightplace = False
      if line.startswith('  iter  continuity'): rightplace = True
      if rightplace:
         rec = line.split()
         if rec and rec[0].isdigit() and len(rec) > 2:
            curr = int(rec[0])
            if curr < first: first = curr
   f.close()
   return pl.array(t), pl.array(iter)

def plot_convergence(trafiles, as_function_of_time=True):
   if type(trafiles) is str: trafiles =[trafiles] # handle string / list
   pl.figure()
   for trafile in trafiles:
      t, iter = read_convergence(trafile)
      if as_function_of_time:
         pl.plot(t, iter, 'o', markersize=6, label=UT.basename(trafile))
      else:
         pl.plot(iter, 'o', markersize=6, label=UT.basename(trafile))
   if as_function_of_time: pl.xlabel('time [sec]')
   else:                   pl.xlabel('time step [-]')
   pl.ylabel('# iterations [-]')
   pl.grid(True)
   pl.legend(loc='best')
   pl.show()

def get_density(bcfile, mat_nm):
   import pyparsing
   line = filter(lambda x: x.startswith('(materials'), open(bcfile).readlines())[0]
   if not line: return -1
   l = pyparsing.OneOrMore(pyparsing.nestedExpr()).parseString(line).asList()
   return float([x[3][1][2] for x in l[0][1] if x[0]==mat_nm][0])

def n_running(joufiles, global_count, verbose):
   '''
   is based on the job being started with option -driver null
   global_count: if True, it will count all fluent jobs - not just the provided cases
   '''
   ps = [x for x in subprocess.check_output(['ps', '-efww']).split('\n') if '-driver' in x]
   if global_count: n = len(ps)
   else:
      n = 0
      for joufile in joufiles:
         for psline in ps:
            if joufile in psline:
               n += 1
               break
   if verbose: print 'number of running cases is', n
   return n

def run(joufiles, dims, max_running=1, version=162, sleep_time=5, verbose=True, global_count=False):
   '''
   joufiles may be just a string (will be converted to list)
   dims is '2d' or '3d'
   max_running is useful when you have a limited number of licenses.
   sleep_time: time to wait until re-checking how many jobs are running [seconds]
   verbose: write some info (boolean)
   global_count: if True, it will count all fluent jobs - not just the provided cases
   note1: will not return until all jobs have been started
   note2: does not support parallel
   '''
   if type(joufiles) is str: joufiles = [joufiles]
   for joufile in joufiles:
      runid = UT.basename(joufile)
      cmd = '/prog/Fluent/ansys_inc/v%i/fluent/bin/fluent -cc %s -i %s.jou -gu -driver null > %s.tra 2>&1 &' % (version, dims, runid, runid)
      if len(joufiles) > max_running:
         while True:
            if n_running(joufiles, global_count, verbose) >= max_running: time.sleep(sleep_time)
            else: break
      print cmd
      os.system(cmd)

class CaseAnalyzer(object):
   '''
   Useful for analyzing (plotting) results from simulations
   Usage:
      ca = FL.CaseAnalyzer(_read_data)
      ca.barplot('time', 'pab00/*.tra', force_read=True)
      ca.xplot('time', 'vf2_hi', 'pab00/f*.tra')
   Note: A 'case' is a tuple like  ('aca10/1', 'aca10', '1')
   '''
#
   def __init__(self, read_func, label_func=None, linest_func=None):
      '''
      read_func   = read_func(case). must return a struct with data (like s.q, s.time etc.)
      label_func  = label_func(case). must return string
      linest_func = linest_func(case). must return string
      '''
      self.read_func   = read_func
      self.label_func  = label_func
      self.linest_func = linest_func
      self._data       = {}
#
   def _cases(self, patterns):
      '''
      gives a list like this:
           label      dir   runid
      [ ('aca10/1', 'aca10', '1'),
        ('aca10/2', 'aca10', '2'),
        ('aca10/3', 'aca10', '3')]
      '''
      files = UT.glob(patterns, sortit=True)
      cases = []
      for fnm in files:
         rec = fnm.split('/')
         runid = UT.basename(rec[-1])
         if len(rec) > 1:
            direct = rec[0]
            lbl = '%s/%s' % (direct,runid)
         else           :
            direct = '.'
            lbl = '%s/%s' % (UT.basename(os.getcwd()),runid)
         cases.append((lbl, direct, runid))
      return cases
#
   def get_data(self, case, force_read):
      if (not force_read) and (self._data.has_key(case[0])): return self._data[case[0]]
      s = self.read_func(case)
      self._data[case[0]] = s
      return s
#
   def _colourid(case, cases, group_mode):
      i = int(array([case == x for x in cases]).nonzero()[0][0]) # postion in list
      if   group_mode == 0:
         return i
      elif group_mode == 1:
         # grouping 1&4, 2&5, 3&6
         n = len(cases)/2
         if i >= n: i = i - n
         return i
      else:
         # grouping 1&2, 3&4, 5&6
         return i/2
#
   def xplot(self, varnm1, varnm2, patterns, group_mode=0, force_read=False): 
      '''
      patterns  : typically 'pab00/*.tra' or ['pab00/*.tra', 'pab01/*.tra']
      group_mode: 0 (no grouping),1 or 2. try changing if grouping is wrong
      '''
      cases = self._cases(patterns)
      pl.figure()
      i = -1
      for case in cases:
         i += 1
         s = self.get_data(case, force_read)
         cid = self._colourid(case, cases, group_mode) if group_mode else i
         x = s.__dict__[varnm1]
         y = s.__dict__[varnm2]
         lbl = self.label_func(case) if self.label_func else s.lbl
         lnst = self.linest_func(case) if self.label_func else '--'
         pl.plot(x, y, color=UT.COLOURS[cid], linestyle=lnst, label=lbl)
      pl.xlabel(varnm1)
      pl.ylabel(varnm2)
      pl.grid(True)
      pl.legend(loc='best')
      pl.show()
#
   def barplot(self, varnm, patterns, tstep=-1, force_read=False):
      '''
      patterns: typically 'pab00/*.tra' or ['pab00/*.tra', 'pab01/*.tra']
      '''
      cases = self._cases(patterns)
      labels = []; y = [] 
      for case in cases:
         s = self.get_data(case, force_read)
         lbl = self.label_func(case) if self.label_func else s.lbl
         labels.append(lbl)
         y.append(s.__dict__[varnm][tstep])
      pl.figure()
      pl.bar(pl.arange(len(cases))+0.2, y)
      pl.xticks(pl.arange(len(cases))+.6, labels)
      pl.ylabel(varnm)
      pl.xlim(0, len(cases)+.3)
      pl.show()

class BC_data(object):
   '''
   reads bc-file into a list of list - using Lisp functionality.
   use method get() to obtain values
   usage:
      bc = BC_data(bcfile)
      density = bc.get(['materials','air','density','constant'])[2]
   '''
   def __init__(self, bcfile):
      import Tools.lisp as lisp
      kws = ('rp', 'dv', 'cx1', 'bc') # these are the sections of a bc-file. skip ni as it causes problems...
      self._data = []
      s = ''
      for line in open(bcfile):
         line = line.strip()
         for kw in kws:
            if line.startswith('(%s'%kw):
               if s: self._data.extend(lisp.parse(s)[1])
               s = ''
         s += line
      self._data.extend(lisp.parse(s)[1])
#
   def get(self, keywords, indx=0, data=None):
      '''
      try to find key (=keywords[indx]) in nested list.
      recursive.
      density = bc.get(['materials','air','density','constant'])[2]
      '''
      if not data: data = self._data
      for x in data:
         if x[0] == keywords[indx]:
            if indx == len(keywords)-1: return x
            else                      : return self.get(keywords, indx+1, x) # recursive - next key
         if type(x[0]) == list:
            r = self.get(keywords, indx, x)                                  # recursive - go deeper
            if r: return r
      return None  # didnt find it

def regular_tsteps_animation(caseid, dname1, pattern, npics, ext='.png', dname2='Anim_regular', skip=0):
   '''
   If simulation is run with varying time-step, this one copies the hardcopies obtained at regular
   time intervals.
   NB! Must have a transcript-file to get time-steps from
   Input:
   caseid           : typically 'a02'
   dname1           : typically 'Anim'
   pattern          : typically 'vof' => Anim/a02_vof_0125.png
   npics            : how many pics to end up with
   dname2           : this is where copies are put
   ext              : hardcopy file extension
   '''
   if not os.path.exists(dname2): os.mkdir(dname2)
   t2 = read_flowtimes('%s.tra'%caseid)
   t = pl.linspace(0,t2[-1], npics)
   n = int(round(float(len(t2)) / len(UT.glob('%s/%s_%s_*%s'%(dname1,caseid,pattern,ext))))) # how often hardcopies are made
   indx = [(i/n)*n for i in get_timesteps(t,t2)]   # only have a hardcopy for every n timestep (so 181 => 180 with n=10)
   for i in range(skip,npics):
      fname1 = '%s/%s_%s_%04i%s'%(dname1,caseid,pattern,indx[i],ext)
      fname2 = '%s/%s_%05i%s'%(dname2,caseid,i,ext) 
      if not os.path.exists(fname2):
         print 'copying %s to %s' % (fname1, fname2)
         shutil.copy(fname1, fname2)
      else:
         print fname2, 'already in', dname2

def simtime(caseid, dname, pattern, ext='.dat', scaler=3600.):
   '''
   assumes file names are like a01_0005.png, which gives tstep = 5
   '''
   t = read_flowtimes('%s.tra'%caseid)
   pattern2 = '%s/%s%s*%s'%(dname,caseid,pattern,ext)
   files =  UT.glob(pattern2, sortit=True)
   if not files:
      print 'found no files matching', pattern2
      return [], []
   stime = []
   cputime = []
   t0 = os.path.getmtime(files[0])
   for f in files:
      stime.append(t[int(UT.basename(f).split('_')[-1])-1])
      cputime.append(int(os.path.getmtime(f)-t0))
   return pl.array(stime), pl.array(cputime)/scaler

def create_dpm_injection_file(fname, points, velocity, diam, massflow=1, temp=U.DEG_ZERO+20, verbose=True):
   '''
   creates a file of named particles. to be used for injections of type 'file'.
   points.
   input:
    points: list of (x,y,z)
    velocity: function s.t. vx,vy,vz = velocity(x,y,z)
   '''
   f = open(fname, 'w')
   np = 0
   for x,y,z in points:
      np += 1
      vx, vy, vz = velocity(x,y,z)
      f.write('(( %g %g %g %g %g %g %g %g %g) P%05i)\n' % (x,y,z, vx,vy,vz, diam, temp, massflow, np))
   f.close()
   if verbose: print fname, 'was created'

# testing code
if __name__ == '__main__':
   f = read_force_report(sys.argv[1])
   print f

   '''
   d = read_bcfile(sys.argv[1])
   kws = ['residuals/settings', 'x-velocity']
   kws = ['adapt/region/zmax']
   kws = ['read-fan-curve']
   kws = ['les-wale-sgs-on?']
   kws = ['unit-table']
   kws = ['materials', 'air', 'density', 'constant']
   x = bc_lookup(d, kws, 0)
   print x
   '''
