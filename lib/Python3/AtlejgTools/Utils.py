import os, sys
import numpy as np
import glob as glob_
import pylab as pl
import re
import tempfile
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
from matplotlib import rcParams
import yaml
import pprint
from pathlib import Path


class Struct(object):
    '''
    a rather flat implempentation of a struct for python.
    note that attribute _keys is special and *cannot* be used.
    '''
    def __init__(self):
        self._keys = []
#
    def has(self, key):
        return hasattr(self, key)
#
    def get(self, key):
        return self.__dict__[key]
#
    def __setattr__(self, key, value):
        if hasattr(self, '_keys') and len(self._keys) > 0 and key == '_keys':
            raise Exception('Cannot use attribute "_keys" for Struct')
        self.__dict__[key] = value
        self._keys.append(key)
#
    def __repr__(self):   # used for __str__ as well
        s = ''
        for key in self._keys:
            if key == '_keys': continue
            s += f'{key:10s} : {pprint.pformat(self.get(key))}\n'
        return s

    
    pass

COLOURS = [x['color'] for x in rcParams['axes.prop_cycle']]    # wanna use same colours as default for 'plot'
COLOURS_OLD = ('b', 'g', 'r', 'c', 'm', 'y', '0.8', '0.65', '0.5', '0.3', '0.1', 'k') # '0.x' gives grayscale
MARKERS = ('', '*', 'o', 's')

# replacing string.join
def join(x, sep=' '):
    return sep.join(x)

class PlotStyle(object):
#
    def __init__(self, ci=0, mi=0):
        self.ci  = ci # color index
        self.mi  = mi # marker index
#
    def __next__(self):
        self.ci += 1
        if self.ci == len(COLOURS):
            self.ci = 0
            self.mi += 1
        if self.mi == len(MARKERS):
            self.mi = 0
#
    def colour(self):
        return COLOURS[self.ci]
#
    def marker(self):
        return MARKERS[self.mi]
#
    def format(self):
        return self.colour() + self.marker()
#
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

def interp_nan(x, y):
    '''
    interpolates nan-values according to neigbours values.
    very basic - uses piecewise linear interpolation
    '''
    ix = pl.isfinite(y)
    return pl.interp(x, x[ix], y[ix])

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
    note: to convert excel-dates to numpy-dates, add t0 = int(1899*365.24) + 4
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
        vals = np.array(vals)
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
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def glob(patterns, sortit=True, suffix=None):
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
    str = str.replace('|','_')
    return str

def hardcopies(fignos, path='./', prefix='', figno_only=False):
    '''
    fignos could be int, in which case it will be interpreted as 1, 2, .., fignos.
    else it should be a list of figure numbers.
    '''
    if type(fignos) is int: fignos = list(range(1, fignos+1))
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
        print('saving hardcopy %s' % fname)
    return fnames

def close_figs(n1, n2=None):
    if n2: figs = list(range(n1, n2+1))
    else : figs = list(range(1, n2+1))
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
    import configparser
    cf = configparser.ConfigParser()
    cf.read(fname)
    sect_nm = list(cf._sections.keys())[0]
    pairs = cf._sections[sect_nm]
    s = Struct()
    for key in list(pairs.keys()):
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
    and it uses the same principle as trapz.
    it handles NaNs by making a copy and replacing them with 0
    '''
    ixs = pl.where(pl.isnan(y))[0]
    if len(ixs):
        y = y.copy()
        y[ixs] = 0.
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
#
    def holds(self, fname):
        return fname in self.coll
#
    def add(self, fname, reread=True, **kwargs):
        '''
        assumes that the supplied read_function handles non-existing file etc.
        reread: if False, will not re-read file even if it has changed. useful for ongoing simulations
        '''
        fname = fullpath(fname)
        timestamp = os.path.getmtime(fname)
        if (not self.holds(fname)) or (reread and self.coll[fname]['timestamp'] < timestamp):
            if self.verbose: print("adding '%s'" % fname)
            data = self.read_function(fname, **kwargs)
            self.coll[fname] = {'data':data, 'timestamp':timestamp}
        elif self.verbose: print("'%s' already here" % fname)
        return self.coll[fname]['data']
#
    def get(self, fname, **kwargs):
        """
        just an alias for add...
        """
        return self.add(fname, **kwargs)
#
    def delete(self, fname):
        if fname not in self.coll:
            print("delete(): did not find key '%s'" % fname)
            return
        self.coll.pop(fname)['data']
#
    def reset(self):
        for key in list(self.coll.keys()):
            self.delete(key)

def tmpfile(dir='.'):
    return tempfile.mkstemp(dir=dir)[1]

def make_collage(pic_files, nrows, ncolumns, newfile='collage.png', convert_path=r'C:\Appl\Cygwin\bin\magick.exe', tmpdir='c:/TMP'):
    '''
    makes a collage of nrows and ncolumns of the given picture files.
    uses ImageMagick convert
    '''
    fnames = []
    for n in range(nrows):
        pics = pic_files[n*ncolumns:(n+1)*ncolumns]
        fname = tempfile.mkstemp(dir=tmpdir)[1]
        cmd = convert_path + ' ' + ' '.join(pics) + ' +append ' + fname
        print(cmd)
        os.system(cmd)
        fnames.append(fname)
    cmd = convert_path + ' ' + ' '.join(fnames) + ' -append ' + newfile
    print(cmd)
    os.system(cmd)
    print(newfile, 'created')
    #[os.unlink(fname) for fname in fnames] # delete intermediate files

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
        self.fname = fname
        if fname: self.read_input_file(verbose, import_statement)
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
        print('creating input file %s' % fname)
#
    def read_input_file(self, verbose=False, import_statement=None):
        if verbose: print('reading input file %s' % self.fname)
        if import_statement: exec(import_statement)
        f = open(self.fname)
        for line in f:
            # handle comments
            line = line.strip().split('#')[0]
            if not line: continue  # there was nothing in front of the '#'
            # handle type
            typ = line.split()[0].strip()
            if typ in ('str', 'eval'): # want to allow spaces in strings. not very elegant code...
                rec = line.split('=')
                varnm = rec[0].split()[1]
                value = join(rec[1:], '=').strip() # reconstruct
                if typ == 'eval': value = eval(value)
            else:
                rec = line.split()
                line = join(rec[1:], '') # reconstruct
                rec = line.split('=')
                varnm = rec[0].strip()
                value = rec[1].strip()
                if verbose: print(varnm, '=', value)
                if   typ == 'float': value = float(value)
                elif typ == 'bool' : value = True if 'true' in value.lower() else False
                elif typ == 'int'  : value = int(value)
                elif typ == 'eval' : value = eval(value)
            self.add_variable(varnm, value, typ)
        f.close()
#
    def has_key(self, key):
        return key in self.__dict__

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

def xkcdify_fonts(fig, ttf='/private/agy/Tools/atlejg/Resources/Misc/Humor-Sans-1.0.ttf'):
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
    #path = '/project/RCP/active/fluent//Atle_Resources/Humor-Sans.ttf'
    fp = font_manager.FontProperties(fname=ttf)
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


def peakdetect(y_axis, x_axis = None, lookahead = 500, delta = 0):
    """
    Atle: Copied from https://gist.github.com/gcalmettes/1784428

    Converted from/based on a MATLAB script at http://billauer.co.il/peakdet.html

    Algorithm for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively

    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- A x-axis whose values correspond to the 'y_axis' list and is used
        in the return to specify the postion of the peaks. If omitted the index
        of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 500)
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the algorithm from picking up false peaks towards to end of
        the signal. To work well delta should be set to 'delta >= RMSnoise * 5'.
        (default: 0)
            Delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the algorithm

    return -- two lists [maxtab, mintab] containing the positive and negative
        peaks respectively. Each cell of the lists contains a tupple of:
        (position, peak_value)
        to get the average peak value do 'np.mean(maxtab, 0)[1]' on the results
    """
    maxtab = []
    mintab = []
    dump = []   #Used to pop the first hit which always if false

    length = len(y_axis)
    if x_axis is None:
        x_axis = list(range(length))

    #perform some checks
    if length != len(x_axis):
        raise ValueError("Input vectors y_axis and x_axis must have same length")
    if lookahead < 1:
        raise ValueError("Lookahead must be above '1' in value")
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError("delta must be a positive number")

    #needs to be a numpy array
    y_axis = np.asarray(y_axis)

    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = np.Inf, -np.Inf

    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x

        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                maxtab.append((mxpos, mx))
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf

        ####look for min####
        if y > mn+delta and mn != -np.Inf:
            #Minima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                mintab.append((mnpos, mn))
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf


    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            maxtab.pop(0)
            #print "pop max"
        else:
            mintab.pop(0)
            #print "pop min"
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass

    return maxtab, mintab

def find_maxima(y, lookahead, delta=0):
    '''
    uses peakdetect to obtain maxima. only returns indices of maxima
    input:
     lookahead: distance to look ahead from a peak candidate to
                determine if it is the actual peak (default: 500)
                '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    '''
    r = peakdetect(y, lookahead=lookahead, delta=delta)[0]
    return [x[0] for x in r]

def find_minima(y, lookahead, delta=0):
    '''
    uses peakdetect to obtain minima. only returns indices of minima
     lookahead: distance to look ahead from a peak candidate to
                determine if it is the actual peak (default: 500)
                '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    '''
    r = peakdetect(y, lookahead=lookahead, delta=delta)[1]
    return [x[0] for x in r]

def plot_multiple(x, ys, ylabels, units, ylims, xlabel='', title='', xoffset=60, loc='best', figsz=(12,7), showit=True):
    '''
    x     : common x-axis variable
    ys    : list of y-axis variables
    units : list of y-labels
    ylims : could be list of []'s
    xlabel: x-label
    showit: True if figure is to be shown (not just a hardcopy - which must be done outside this routine)
    '''
    offset = r_[xoffset, 0]
    fig = figure(figsize=figsz)
    host = HostAxes(fig, [0.05, 0.1, 0.65, 0.8])
    host.axis["right"].set_visible(False)
    host.set_yticks([])
    host.set_xlabel(xlabel)
    if title: host.set_title(title)
    n = 0
    for y, ylbl, unit, ylim in zip(ys, ylabels, units, ylims):
        n += 1
        paras = ParasiteAxes(host, sharex=host)
        host.parasites.append(paras)
        pl, = paras.plot(x, y, label=ylbl)
        paras.set_ylabel('%s [%s]'%(ylbl, unit))
        if ylim: paras.set_ylim(ylim)
        if n == 1:
            ax = paras.axis["right"]
            paras0 = paras               # need it later...
        else:
            new_axisline = paras._grid_helper.new_fixed_axis
            paras.axis["right2"] = new_axisline(loc="right", axes=paras, offset=(n-1)*offset)
            ax = paras.axis["right2"]
        ax.label.set_color(pl.get_color())
        ax.major_ticklabels.set_visible(True)
        ax.set_visible(True)
    #
    # dont know why i need to do this to show ylabel for the first variable, but i need to...
    paras0.set_ylabel('%s [%s]'%(ylabels[0], units[0]))
    paras0.axis["right"].label.set_visible(True)
    #
    fig.add_axes(host)
    host.legend(loc=loc)
    if showit: show()

def get_yaml(fnm):
    '''
    read yaml-file into a Struct for easy access.
    and also avoid the open()
    '''
    yml = yaml.load(open(fnm), Loader=yaml.SafeLoader)
    s = Struct()
    for k,v in yml.items():
        s.__dict__[k] = v
    return s




def _get_size(bytes):
    """
    Returns size of bytes in a nice format
    """
    for unit in ['', 'K', 'M', 'G', 'T', 'P']:
        if bytes < 1024:
            return f"{bytes:.2f} {unit}B"
        bytes /= 1024


def get_processes_info(sort_by, descending=True, columns=None):
    '''
    copied / based on
    https://www.thepythoncode.com/code/make-process-monitor-python
    # the data-frame that contain all process info
    '''
    import psutil
    from datetime import datetime
    import pandas as pd
    import time
    #
    if columns is None:
        columns = "name,cpu_usage,memory_usage,read_bytes,write_bytes,status,create_time,nice,n_threads,cores"
    #
    processes = []
    for process in psutil.process_iter():
        # get all process info in one shot
        with process.oneshot():
            # get the process id
            pid = process.pid
            if pid == 0:
                # System Idle Process for Windows NT, useless to see anyways
                continue
            # get the name of the file executed
            name = process.name()
            # get the time the process was spawned
            try:
                create_time = datetime.fromtimestamp(process.create_time())
            except OSError:
                # system processes, using boot time instead
                create_time = datetime.fromtimestamp(psutil.boot_time())
            try:
                # get the number of CPU cores that can execute this process
                cores = len(process.cpu_affinity())
            except psutil.AccessDenied:
                cores = 0
            # get the CPU usage percentage
            cpu_usage = process.cpu_percent()
            # get the status of the process (running, idle, etc.)
            status = process.status()
            try:
                # get the process priority (a lower value means a more prioritized process)
                nice = int(process.nice())
            except psutil.AccessDenied:
                nice = 0
            try:
                # get the memory usage in bytes
                memory_usage = process.memory_full_info().uss
            except psutil.AccessDenied:
                memory_usage = 0
            # total process read and written bytes
            io_counters = process.io_counters()
            read_bytes = io_counters.read_bytes
            write_bytes = io_counters.write_bytes
            # get the number of total threads spawned by this process
            n_threads = process.num_threads()
            # get the username of user spawned the process
            try:
                username = process.username()
            except psutil.AccessDenied:
                username = "N/A"

        processes.append({
            'pid': pid, 'name': name, 'create_time': create_time,
            'cores': cores, 'cpu_usage': cpu_usage, 'status': status, 'nice': nice,
            'memory_usage': memory_usage, 'read_bytes': read_bytes, 'write_bytes': write_bytes,
            'n_threads': n_threads, 'username': username,
        })
    #
    # convert to pandas dataframe
    pcs = pd.DataFrame(processes)
    # set the process id as index of a process
    pcs.set_index('pid', inplace=True)
    # sort rows by the column passed as argument
    pcs.sort_values(sort_by, inplace=True, ascending=not descending)
    # pretty printing bytes
    pcs['memory_usage'] = pcs['memory_usage'].apply(_get_size)
    pcs['write_bytes'] = pcs['write_bytes'].apply(_get_size)
    pcs['read_bytes'] = pcs['read_bytes'].apply(_get_size)
    # convert to proper date format
    pcs['create_time'] = pcs['create_time'].apply(datetime.strftime, args=("%Y-%m-%d %H:%M:%S",))
    # reorder and define used columns
    pcs = pcs[columns.split(",")]
    return pcs

def read_graph_grabber_data(fnm):
    '''
    reads csv-file from GraphGrabber, with multiple data series, into a Struct
    '''
    lines = open(fnm).readlines()
    #
    series = []
    #
    for no, line in enumerate(lines):
        line = line.strip()
        if not line: continue
        if 'Series' in line:
            if no > 0: s.x, s.y = np.array(xs), np.array(ys)
            lbl = ' '.join(line.split()[2:])
            s = Struct()
            s.label = lbl
            series.append(s)
            xs, ys = [], []
            continue
        x, y = line.split(';')
        xs.append(float(x))
        ys.append(float(y))
    #
    s.x, s.y = np.array(xs), np.array(ys)
    return series

