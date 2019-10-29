import sys, tempfile, os, shutil, time

class Figure(object) :

    def __init__(self) :
        '''
        creator
        '''
        self._data = []
        self._legs = []
        #self._tmpdir = tempfile.mkdtemp(prefix='gnuplot_') # not availabel on 2.2.3 :-(
        self._tmpdir = tempfile.mktemp('.gnuplot'); os.mkdir(self._tmpdir)
        self._figtype = 'png'
        self._cwd     = os.getcwd()
        self._title   = 'notitle'
        self._xlabel  = 'xvar'
        self._ylabel  = 'yvar'
        self._xmin = None
        self._xmax = None
        self._ymin = None
        self._ymax = None
        self._location = 'top' # legend location

    def title(self, str=None) :
        '''
        title
        '''
        if not str == None : self._title = str
        return self._title

    def location(self, str=None) :
        '''
        location. could be top, bottom, left, right etc
        '''
        if not str == None : self._location = str
        return self._location

    def xlabel(self, str=None) :
        '''
        xlabel
        '''
        if not str == None : self._xlabel = str
        return self._xlabel

    def ylabel(self, str=None) :
        '''
        ylabel
        '''
        if not str == None : self._ylabel = str
        return self._ylabel


    def filename(self, relative_path) :
        '''
        get filename for _this_ figure
        '''
        if self._title.find('/') >= 0 :
        #if self._title.find(os.path.sep) >= 0 : # not available on 2.2.3 :-(
            filename = 'notitle'
            print("warning : cannot use title (%s) as filename. using '%s' instead" % (self._title, filename))
        else : filename = self._title
        return os.path.join(self._cwd, relative_path, filename.replace(' ', '_'))

    def axis(self, xmin, xmax, ymin, ymax) :
        '''
        set xmin, xmax, ymin,ymax
        '''
        self._xmin = xmin
        self._xmax = xmax
        self._ymin = ymin
        self._ymax = ymax

    def load(self, x, y, leg) :
        '''
        add data to be plotted in _this_ figure
        '''
        self._data.append((x, y))
        self._legs.append(leg)

    def cleanup(self) :
        '''
        deleting files and clearing data
        '''
        time.sleep(1) # make sure gnuplot has finished its business in the directory
        shutil.rmtree(self._tmpdir)
        self._data = []
        self._legs = []

    def _make_plot(self, hardcopy, relative_path='.') :
        # create the plot
        #
        if hardcopy :
            term = self._figtype
            options = ''
        else :
            term = 'x11'
            options = '-persist'
        gpfile = open(os.path.join(self._tmpdir, 'run.gnuplot'), 'w')
        print("Gnuplot::_make_plot : gpfile= " + gpfile.name)
        if self._xmin != None :
            gpfile.write("""\
   set xrange [%e:%e]
   set yrange [%e:%e]
            """ % (self._xmin, self._xmax, self._ymin, self._ymax))
        gpfile.write("""
  set grid
  set title  "%s"
  set xlabel "%s"
  set ylabel "%s"
  set key %s
  set term %s
  set output '%s.%s';
  plot \\\n""" \
        % (self._title, self._xlabel, self._ylabel, self._location, term, self.filename(relative_path), self._figtype))
        i = 0
        for (x, y) in self._data :
            fname = os.path.join(self._tmpdir,'%i.data' % i)
            _write_data(x, y, fname)
            gpfile.write("""\
   '%s' title '%s' with linespoints pointtype 5 \
            """ % (fname, self._legs[i]))
            i += 1
            if i < len(self._data) : gpfile.write(', \\\n  ')
        gpfile.close()
        #
        cmd = 'gnuplot %s %s' % (options, gpfile.name)
        failure = os.system(cmd)
        if failure    : print("Gnuplot.hardcopy(): running cmd = %s failed" % cmd)
        elif hardcopy : print("Gnuplot.hardcopy(): created file %s.%s" % \
              (self.filename(relative_path), self._figtype))

    def hardcopy(self, relative_path='.') :
        '''
        obtain a hardcopy of _this_ figure
        '''
        self._make_plot(True, relative_path)

    def show(self) :
        '''
        obtain a gnuplot of _this_ figure
        '''
        self._make_plot(False)

def _write_data(x, y, fname) :
    file = open(fname, 'w')
    for i in range(0, len(x)) :
        file.write('%e %e\n' % (x[i], y[i]))
    file.close()


if __name__ == '__main__' :

    x = list(range(0,5))
    y = list(range(1,6))

    gp = Figure()
    gp.load(x,y,'test1')
    gp.load(y,x,'test2')
    gp.title('test title')
    gp.hardcopy()
    gp.cleanup()
