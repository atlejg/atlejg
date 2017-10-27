import os, sys
import numpy
import glob as glob_
import pylab as pl
import re
import tempfile
import string

class Struct: pass

COLOURS = ('b', 'g', 'r', 'c', 'm', 'y', 'k', '0.3', '0.6', '0.9') # '0.x' gives grayscale

class PlotStyle(object):
   colours = COLOURS
   markers = ('-', '*', '*-', '.', '--*')

   def __init__(self, ci=0, mi=0):
      self.ci  = ci # color index
      self.mi  = mi # marker index

   def next_marker(self):
      self.mi = (self.mi + 1) % len(PlotStyle.markers)

   def next_color(self):
      self.ci = (self.ci + 1) % len(PlotStyle.colours)

   def prev_marker(self):
      self.mi = (self.mi - 1) % len(PlotStyle.markers)

   def prev_color(self):
      self.ci = (self.ci - 1) % len(PlotStyle.colours)

   def fmt(self, shift_color=True, shift_marker=False):
      if shift_color : self.next_color()
      if shift_marker: self.next_marker()
      str = PlotStyle.colours[self.ci] + PlotStyle.markers[self.mi]
      return str

   def reset(self):
      self.ci = 0
      self.mi = 0

def replace_in_file(replace,origfile,newfile):
   '''
   replace is a list of 2-tuples (or lists) where the first element is replaced by the second.
   '''
   of = open(origfile)
   nf = open(newfile,'w')
   for line in of:
      for r in replace:
         line = line.replace(r[0],r[1])
      nf.write(line)
   of.close()
   nf.close()
   return newfile

def run_in_new_thread(): pass

def tcsh(cmd) :
   cmd = "tcsh -c '%s'" % cmd
   #print cmd
   os.system(cmd)

# some useful maths for standard python
def interpolate(x, y, x0, index=None) :
   '''
   assumes y is a linear function of x on the intervals between points.
   x must be monotonically increasing.
   extrapolates linearly from the last two points in the series.
   using plain python - so lists are allowed. but not very effective...
   '''
   # if index was not provided, we need to find it ourselves
   if index == None :
      for i in range(len(x)) :
         if x[i] > x0 : break
   index = i-1
   #
   # find interpolated value
   a = (y[index+1]-y[index]) / (x[index+1]-x[index])
   return y[index] + a*(x0-x[index])

def find_roots(x, y, last_only=False) :
   roots   = []
   indices = []
   for i in range(0, len(x)-1) :
      if y[i] * y[i+1] < 0 :
         a = (y[i+1]-y[i]) / (x[i+1]-x[i])
         b = y[i] - a*x[i]
         roots.append(-b/a)
         indices.append(i)
   if last_only and roots : return (roots[-1], indices[-1])
   else                   : return (roots, indices)

def column_indx(col_nm):
   '''col_nm : excel-name of column (a,b, ... aa, ab, ... etc)'''
   n = len(col_nm)
   indx = 0
   for x in col_nm:
      n -= 1
      indx += (ord(x.lower()[-1])-96) * 26**n
   return indx - 1

def read_excel_column(sheet, col_nm, line_from, line_to, file=None, book=None, return_refs=False, return_array=False):
   '''
   if file and book is None, it will assume sheet is book.sheet (performance gain)
   if sheet is integer (int) it will use sheet_by_index, else sheet_by_name
   '''
   import xlrd
   if file != None:
      book = xlrd.open_workbook(file)
   if book != None:
      if type(sheet) == int:
         sh = book.sheet_by_index(sheet)
      else:
         sh = book.sheet_by_name(sheet)
   else: sh = sheet
   vals =  sh.col_values(column_indx(col_nm), line_from-1, line_to)
   if return_array :
      import numpy
      vals = numpy.array(vals)
   if return_refs: return vals, book, sheet
   else          : return vals


def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    found on www.scipy.org/Cookbook/SignalSmooth. -atle,
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def glob(patterns, sortit=False, suffix=None):
   '''
   Just a wrapper around glob to make it accept lists. Useful when running scripts under
   ipython which doesnt expand filenames from wildcards.
   note: when using suffix, include the '.' (like suffix='.txt')
   '''
   # if is input is a string we'll make it a list
   if type(patterns) is str: patterns = (patterns, )
   files = []
   for pattern in patterns :
      if not suffix is None: pattern += suffix
      files.extend(glob_.glob(pattern))
   if sortit: files.sort()
   return files

def clean_string_for_fname(str):
   # replace space, brackets, and slashes with underscores
   str = str.replace(' ','_')
   str = str.replace('(','_')
   str = str.replace(')','_')
   str = str.replace('[','_')
   str = str.replace(']','_')
   str = str.replace('/','_')
   str = str.replace('?','_')
   str = str.replace(':','_')
   return str

def hardcopies(fignos, path='./', prefix='', figno_only=False):
   '''
   fignos could be int, in which case it will be interpreted as 1, 2, .., fignos.
   else it should be a list of figure numbers.
   '''
   if type(fignos) is int: fignos = range(1, fignos+1)
   fnames = []
   fignos = pl.array(fignos)
   for figno in fignos:
      pl.figure(figno)
      if figno_only: titl = ''
      else         : titl = pl.get(pl.gca(), 'title')
      titl = clean_string_for_fname(titl)
      fname = ('%s/%s%02i_%s.png' % (path, prefix, figno, titl))
      fnames.append(fname)
      pl.savefig(fname)
      print 'saving hardcopy %s' % fname
   return fnames

def close_figs(n1, n2=None):
   if n2: figs = range(n1, n2+1)
   else : figs = range(1, n2+1)
   for figno in figs:
      pl.figure(figno)
      pl.close()

def basename(fnames):
   '''
   ut.basename('/project/RCP/active/test.txt') = 'test'
   '''
   if not (type(fnames) is list or type(fnames) is tuple): fnames = (fnames,)
   basenames = []
   for fname in fnames:
      basenames.append(os.path.splitext(os.path.basename(fname))[0])
   if len(basenames) == 1: return basenames[0]
   return basenames

def savetxt(fname, data, **kwargs):
   '''saves data as a columns to text file.
   data could be matrix or list of vectors'''
   if isinstance(data, pl.ndarray):
      # assuming its a matrix
      m = data.reshape(data.shape[-1], -1)
   else:
      # assuming its a list of vectors
      i = 0
      for v in data:
         if i == 0: m = v
         else     : m = pl.vstack((m, v))
         i += 1
      m = m.T
   pl.savetxt(fname, m, **kwargs)

def savecsv(fname, data, fmt='%.18e'):
   savetxt(fname, data, fmt=fmt, delimiter=';')

def split_by_length(s, length, stripit=False):
   '''split given string into a list of strings by chopping the string in regular parts'''
   strings = []
   n1 = 0
   while True:
      n2 = min(len(s), n1 + length)
      if stripit: strings.append(s[n1:n2].strip())
      else      : strings.append(s[n1:n2])
      if n2 == len(s): return strings
      n1 += length

def grep_column(fname, expr, column_no, is_float=True, skip=0):
   '''returns an array with all numbers sitting in the given column number in the given file where the
   given expression matches as a regular expression'''
   file = open(fname)
   v = []
   i = -1
   for line in file:
      i += 1
      if i < skip: continue
      match = re.search(expr, line)
      if match: 
         if is_float: v.append(float((line.split())[column_no-1]))
         else       : v.append((line.split())[column_no-1])
   if is_float: return pl.array(v)
   else       : return v

def read_inifile(fname):
   '''
   very basic inifile reader. assumes only one section and
   a flat structure. useful for replacing sys.argv when calling scripts.
   returns a Struct with the keys as members.
   eg.
   [params]
   root = .
   gives ini.root = '.'
   '''
   import ConfigParser
   cf = ConfigParser.ConfigParser()
   cf.read(fname)
   sect_nm = cf._sections.keys()[0]
   pairs = cf._sections[sect_nm]
   s = Struct()
   for key in pairs.keys():
      s.__dict__[key] = pairs[key]
   return s

def cumulative_old(t, y):
   '''
   bug: first 2 values returned are 0. the first one should be disregarded,
        and instead one value is "added" to the end of the array
   2015-03-31: fixed bug.
   '''
   return pl.array([pl.trapz(y[:n], t[:n]) for n in pl.arange(1, len(t)+1)])

def cumulative(t, y):
   '''
   *much* more effective than cumulative_old,
   and it uses the same principle as trapz
   '''
   cy = pl.zeros(len(t))
   for i in pl.arange(1,len(t)):
      cy[i] = cy[i-1] + (y[i-1]+y[i])/2.*(t[i]-t[i-1])
   return cy

def grep(pattern, items):
   '''
   search for pattern among items (which is iterable list) and return those that
   match (as regular expression)
   '''
   returnlist = []
   for item in items:
      if re.search(pattern, item): returnlist.append(item)
   return returnlist

class CacheManager(object):
   """
   simple class for caching data objects read from files.
   useful for things like Eclipse summary-readers.
   """
   def __init__(self, read_function, verbose=False):
      """
      must supply a function for reading data object from file.
      """
      self.coll          = {}              # collection of data-objects
      self.read_function = read_function   # function for obtaining data-objects from file
      self.verbose       = verbose         # silent or not?
   def add(self, fname, reread=True, **kwargs):
      '''
      assumes that the supplied read_function handles non-existing file etc.
      reread: if False, will not re-read file even if it has changed. useful for ongoing simulations
      '''
      fname = fullpath(fname)
      timestamp = os.path.getmtime(fname)
      if (not self.coll.has_key(fname)) or (reread and self.coll[fname]['timestamp'] < timestamp):
         if self.verbose: print "adding '%s'" % fname
         data = self.read_function(fname, **kwargs)
         self.coll[fname] = {'data':data, 'timestamp':timestamp}
      elif self.verbose: print "'%s' already here" % fname 
      return self.coll[fname]['data']
   def get(self, fname, **kwargs):
      """
      just an alias for add...
      """
      return self.add(fname, **kwargs)
   def delete(self, fname):
      if not self.coll.has_key(fname):
         print "delete(): did not find key '%s'" % fname
         return
      self.coll.pop(fname)['data']
   def reset(self):
      for key in self.coll.keys():
         self.delete(key)

def tmpfile(dir='.'):
   return tempfile.mkstemp(dir=dir)[1]

def make_collage(pic_files, nrows, ncolumns, newfile='collage.png'):
   '''
   makes a collage of nrows and ncolumns of the given picture files.
   uses ImageMagick convert
   '''
   fnames = []
   for n in range(nrows):
      pics = pic_files[n*ncolumns:(n+1)*ncolumns]
      fname = tempfile.mkstemp()[1]
      cmd = '/usr/bin/convert +append ' + string.join(pics, ' ') + ' ' + fname
      os.system(cmd)
      fnames.append(fname)
   cmd = '/usr/bin/convert -append ' + string.join(fnames, ' ') + ' ' + newfile
   os.system(cmd)
   print newfile, 'created'
   [os.unlink(fname) for fname in fnames] # delete intermediate files

def each_figure_do(fignos, cmd):
   for figno in fignos:
      pl.figure(figno)
      eval(cmd)

def combine_plots2(fig_numbers):
   '''
   given a list of figure numbers, it will plot
   all data in a new figure, using attributes from existing plots.
   '''
   newfig = pl.figure().number
   for fign in fig_numbers:
      pl.figure(fign)
      lines = pl.gca().get_lines()
      pl.figure(newfig)
      for line in lines:
         xy = line.get_xydata()
         pl.plot(xy[:,0], xy[:,1],
                 label=line.get_label(),
                 marker=line.get_marker(),
                 linestyle=line.get_linestyle(),
                 color=line.get_color(),
                 )
   a = pl.figure(fig_numbers[0]).get_axes()[0]
   pl.figure(newfig)
   pl.legend(loc='best')
   pl.grid(True)
   pl.title(a.get_title())
   pl.xlabel(a.get_xlabel())
   pl.ylabel(a.get_ylabel())
   pl.show()

def combine_plots(fig_numbers, return_data=False, marker='None', linestyle='-'):
   '''
   given a list of figure numbers, it will plot
   all data in a new figure.
   if asked for, it will return the data (and labels)
   if no line is wanted - use linestyle=''
   '''
   xydata = []
   labels = []
   for fign in fig_numbers:
      pl.figure(fign)
      lines = pl.gca().get_lines()
      for line in lines:
         xydata.append(line.get_xydata())
         labels.append(line.get_label())
   pl.figure()
   for d, lbl in zip(xydata, labels):
      pl.plot(d[:,0], d[:,1], label=lbl, marker=marker, linestyle=linestyle)
   pl.legend(loc='best')
   pl.grid(True)
   if return_data: return xydata, labels

def replot(fig_number, attrib, value):
   '''
   will plot all data in a figure into a new figure, but now with a new attribute - 
   like replot(4, 'marker', '-')
   '''
   pl.figure(fig_number)
   lines = pl.gca().get_lines()
   newfig = pl.figure().number
   for line in lines:
      xy = line.get_xydata()
      pl.plot(xy[:,0], xy[:,1], label=line.get_label(), **{attrib:value})
   a = pl.figure(fig_number).get_axes()[0]
   pl.figure(newfig)
   pl.legend(loc='best')
   pl.title(a.get_title())
   pl.xlabel(a.get_xlabel())
   pl.grid(True)
   pl.show()

def replot2(fig_number, **kwargs):
   '''
   will plot all data in a figure into a new figure, but now with a new attributes
   like UT.replot2(5, **{'linestyle':'', 'marker':'*'})
   '''
   pl.figure(fig_number)
   lines = pl.gca().get_lines()
   newfig = pl.figure().number
   for line in lines:
      xy = line.get_xydata()
      pl.plot(xy[:,0], xy[:,1], label=line.get_label(), **kwargs)
   a = pl.figure(fig_number).get_axes()[0]
   pl.figure(newfig)
   pl.legend(loc='best')
   pl.title(a.get_title())
   pl.xlabel(a.get_xlabel())
   pl.grid(True)
   pl.show()

class InputValues(object):
   '''
   reads input variables from a simple text file, which is often more useful than reading
   from command line or STDIN.
   text file must look like
   # case 1
   float length = 3.                    # [m]
   float width  = 1.                    # [m]
   int   index  = 2
   str   direct = ../..                 # note: this is a string
   bool  force  = True
   eval  vals1  = [1,3,5,6]             # list
   eval  vals2  = range(5)              # list
   eval  q      = pl.linspace(0,1, 11)  # use pl to get pylab-things
   ...
   note1: variable names cannot start or end with white-space characters.
          same goes with values.
   note2: comments are allowed (see example)
   note3: could use add_variable to create variables and create_input_file to make the file.
   note4: valid types are: float, str, int, bool and eval
   note5: eval could be used for lists (or arrays)
   note5: consider using 'import * from constants' instead... (put variables in file constants.py)
   '''
   def __init__(self, fname='', verbose=False, import_statement=None):
      self.varnms = []
      self.types  = {}
      if fname: self.read_input_file(fname, verbose, import_statement)
#
   def __str__(self, fname=''):
      s = ''
      for varnm in self.varnms:
         s += '%s %s = %s\n' % (self.types[varnm], varnm, str(self.__dict__[varnm]))
      return s
#
   def add_variable(self, varnm, value, typ):
      self.__dict__[varnm] = value
      self.types[varnm] = typ
      if not varnm in self.varnms: self.varnms.append(varnm)
#
   def create_input_file(self, fname):
      f = open(fname, 'w')
      for varnm in self.varnms:
         f.write('%s %s = %s\n' % (self.types[varnm], varnm, str(self.__dict__[varnm])))
      f.close()
      print 'creating input file %s' % fname
#
   def read_input_file(self, fname, verbose=False, import_statement=None):
      if verbose: print 'reading input file %s' % fname
      if import_statement: exec(import_statement)
      f = open(fname)
      for line in f:
         # handle comments
         line = line.strip().split('#')[0]
         if not line: continue  # there was nothing in front of the '#'
         # handle type
         typ = line.split()[0].strip()
         if typ in ('str', 'eval'): # want to allow spaces in strings. not very elegant code...
            rec = line.split('=')
            varnm = rec[0].split()[1]
            value = string.join(rec[1:], '=').strip() # reconstruct
            if typ == 'eval': value = eval(value)
         else:
            rec = line.split()
            line = string.join(rec[1:], '') # reconstruct
            rec = line.split('=')
            varnm = rec[0].strip()
            value = rec[1].strip()
            if verbose: print varnm, '=', value
            if   typ == 'float': value = float(value)
            elif typ == 'bool' : value = bool(eval(value)) # note: bool('0') is True...
            elif typ == 'int'  : value = int(value)
            elif typ == 'eval' : value = eval(value)
         self.add_variable(varnm, value, typ)
      f.close()

def perturb(v):
   '''
   perturb the values of the list v into a new list.
   '''
   pv = []
   while len(v) > 0:
      i = pl.randint(len(v))
      pv.append(v.pop(i))
   return pv

def elbow_func(p1, p2, p3, x):
   '''
   linear elbow function: assuming linear between
   three points p1, p2, p3 (pN = (x,y)).
   returns value at point x.
   assumes p1[0] <= x <= p3[0] and p1[0] <= p2[0] <= p3[0]
   '''
   if x > p2[0]: p1, p2 = (p2, p3)   # smart trick
   a = (p2[1]-p1[1]) / (p2[0]-p1[0]) # calculate stigningstallet
   return p1[1] + a*(x-p1[0])

_phi    = (1 + pl.sqrt(5)) / 2.  # 1.6180339887498949
_goldenr = 2 - _phi              # 0.3819660112501051
def find_max(f, a, b, c, dx, *args):
   ''' 
   based on: http://en.wikipedia.org/wiki/Golden_section_search
   a and c are the current bounds; the maximum is between them.
   b is a center point (a < b < c). b could be None
   note that the function is not necessarily evaluating at a and c
   so make sure the maximum is not an endpoint.
   f(x) is some mathematical function - it must be unimodal
   '''
   if b == None: 
      b = a + _goldenr*(c-a)
   if c-b > b-a:
      x = b + _goldenr*(c-b)
      righty = True            # b is closer to c
   else:
      x = b - _goldenr*(b-a)
      righty = False
   if c-a < dx:
      x = (c+a) / 2. 
      return (x, f(x, *args))
   #assert(f(x) != f(b)) # !!! should do something here..
   if f(x, *args) > f(b, *args):
      if righty: return find_max(f, b, x, c, dx, *args)
      else:      return find_max(f, a, x, b, dx, *args)
   else:
      if righty: return find_max(f, a, b, x, dx, *args)
      else:      return find_max(f, x, b, c, dx, *args)

def find_min(f, a, b, c, dx, *args):
   find_max(-1*f, a, b, c, dx, *args)

def npv_scale(y, t, intr_rate):
   '''
   t in days
   intr_rate in %
   '''
   y_npv = [y[n]/(1+intr_rate/100.)**(t[n]/365.25) for n in pl.arange(len(y))]
   return pl.array(y_npv)

'''
better to use xkcdify below!!
def xkcd_plot(*args, **kwargs):
   #for some reason i have a hard time getting fonts right (want to use 'Humor Sans').
   #so i do it here...
   import matplotlib
   import matplotlib.pyplot as plt
   import matplotlib.font_manager as font_manager
   fig = plt.figure()
   path = '/project/RCP/active/fluent//Atle_Resources/Humor-Sans.ttf'
   fp = font_manager.FontProperties(fname=path)
   plt.xkcd()
   ax = plt.plot(*args, **kwargs)
   plt.xticks(plt.xticks()[0], fontproperties=fp)
   plt.yticks(plt.yticks()[0], fontproperties=fp)
   plt.show()
   plt.rcdefaults()
   return fig.get_axes()[0]
'''

def xkcdify_fonts(fig):
   '''
   for some reason i have a hard time getting fonts right (want to use 'Humor Sans').
   so i do it here...
   2017-09-21: seems like problem is gone. try to just use xkcd()
   remember doing xkcd() before creating the plot.
   example:
   xkcd()
   rcParams['lines.linewidth'] = 5
   fig = figure()
   plot(x,ysim, label='simulated')
   plot(x,yexp, label='exp')
   xkcdify_fonts(fig)
   '''
   import matplotlib
   import matplotlib.pyplot as plt
   import matplotlib.font_manager as font_manager
   path = '/project/RCP/active/fluent//Atle_Resources/Humor-Sans.ttf'
   fp = font_manager.FontProperties(fname=path)
   plt.figure(fig.number) # sets gca()
   #
   # loop all (sub)plots
   for ax in fig.get_axes():
      # titles:
      ax.set_title(ax.get_title(), fontproperties=fp)
      # xlabels:
      lbls = [txt.get_text() for txt in ax.get_xticklabels()]
      ax.set_xticklabels(lbls, fontproperties=fp)
      ax.set_xlabel(ax.get_xlabel(), fontproperties=fp)
      # ylabels:
      lbls = [txt.get_text() for txt in ax.get_yticklabels()]
      ax.set_yticklabels(lbls, fontproperties=fp)
      ax.set_ylabel(ax.get_ylabel(), fontproperties=fp)
      #
      # legends
      l = ax.get_legend()
      if l:
         txts = l.get_texts()
         lbl = []
         for txt in txts:
            lbl.append(txt.get_text())
         ax.legend(lbl, prop=fp, loc=l._loc_real)
   #
   # done
   plt.show()
   return fp # in case someone needs it

def rand_unif(min_val, max_val, nvals, sorted=True):
   '''
   gives a vector of random variables in the given interval. uniformly distributed.
   '''
   vals = pl.rand(nvals)
   vals *= (max_val - min_val)
   vals += min_val
   if sorted: vals.sort()
   return vals

def rand_norm(meanval, stddev, nvals, sorted=True, low_lim=None, hi_lim=None):
   '''
   gives a vector of normal distributed values.
   if asked for, it will force all values to be within the interval [low_lim, hi_him] (or just one of these)
   '''
   vals = pl.normal(meanval, stddev, nvals)
   if not low_lim is None: vals[vals <= low_lim] = low_lim
   if not hi_lim is None:  vals[vals >= hi_lim]  = hi_lim
   if sorted: vals.sort()
   return vals

def fullpath(fname):
   '''
   makes sure file-name has full path
   '''
   if fname.startswith('/'): return fname
   else:
      return '%s/%s' % (os.getcwd(), fname)

def find_index(t, t0, eps=1e-3):
   ''' find index where the array/list t is equal (or very close) to t0'''
   indx = 0
   for t_ in t:
      if abs(t_ - t0) < eps: return indx
      indx += 1
   raise Exception('did not find index') # should not get here

def find_closest_index(t, t0):
   ''' find index where the array/list t is closest to t0'''
   indx = 0
   dt_min = pl.Inf
   i = -1
   for t_ in t:
      i += 1
      dt = abs(t_ - t0)
      if dt < dt_min:
         dt_min = dt
         indx = i
   return indx

def clone(obj):
   '''
   doesnt work - why??
   '''
   import pickle
   f = open('t', 'w')
   pickle.dump(obj,f)
   f.close()
   f = open('t')
   cloned = pickle.load(f)
   f.close()
   return cloned
