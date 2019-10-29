from scipy.interpolate import interp1d
import pylab as pl
import codecs
import AtlejgTools.Utils                          as UT
import AtlejgTools.SimulationTools.UnitConversion as U

# i think these are constants. or must be found from mb.segments.casing_liner._option.__dict__.keys()
ID_BLANK = 7.
ID_ICD   = 30.

class Container(object):
    '''
    slightly intelligent struct
    '''
    def __init__(self):
        self.val = None
#
    def get_array(self):
        return pl.array([float(x.strip()) for x in self.val.split()])
#
    def get_branches(self, bnames):
        '''
        for ease of access
         note:
          may only be used on the top node!
         input:
          bnames: name of brances. typically ['mainbore', 'lateral 1', 'lateral 2']
         output:
          list of brances with measured depths, trajectories, screens etc.
        '''
        if not 'well' in list(self.__dict__.keys()):
            print('get_branches may only be called on the top node!')
            return None
        bs = []
        for bname in bnames:
            b = self.well.__dict__[bname]
            # md
            b.mds = b.nodes.md.get_array()[:-1]    # md-vector is one too long...
            b.ns = len(b.mds)                      # number of segments. for convinience
            # permeability [mD] that Netool uses
            b.ks = b.segments.kh.get_array() / U.DARCY * 1000.
            # path of branch
            xs = b.nodes.x.get_array()[:-1]
            ys = b.nodes.y.get_array()[:-1]
            zs = b.nodes.z.get_array()[:-1]
            b.path = pl.array(list(zip(xs, ys, zs)))
            # - from log (high resolution)
            b.md_log = self.well.mainbore.segments.kh._log.x.get_array()
            b.k_log  = self.well.mainbore.segments.kh._log.y.get_array() / U.DARCY * 1000.  # mD
            # AICD settings
            b.xs = self.well.mainbore.segments.aicd_x.get_array()
            # which segment is a screen? the others can be seen as blanks..
            mask = b.segments.inflow_control.get_array()==ID_ICD
            b.screens = pl.zeros(b.ns)                 # screens (0 or 1)
            b.screens[mask] = 1
            bs.append(b)
        return bs

def closest_indx(b, pos):
    '''
    finds index where the given position is closest to
    the path of the given branch
    b  : branch
    pos: position
    '''
    dists = norm(b.path - pos, axis=1)
    return argmin(dists)

def block_data(xs, ys, x1, x2, nx):
    '''
    note: it uses an average value for a segment covering before
          and after the x-values. this means x1 > min(xs) and x2 < max(xs)
    '''
    dx = (x2-x1) / nx
    ycs = UT.cumulative(xs, ys)
    yc_func = interp1d(xs, ycs, fill_value='extrapolate')
    xx = pl.linspace(x1-dx/2., x2+dx/2., nx+2)
    yy = yc_func(xx)
    return xx[:-1], pl.diff(yy) / dx

def _fix_varnm(varnm):
    varnm = varnm.strip().replace(' ', '_').replace('/', '_')  # cannot have spaces or '/'
    if varnm[0] in '0123456789': varnm = '_' + varnm           # cannot start with digit
    # treat special characters..
    varnm = varnm.replace('\u03bc', 'mu')   # viscosity -> mu
    varnm = varnm.replace('\u03c1', 'rho')  # density -> rho
    return varnm  # default

def read_ntlcase(ntlcase):
    lines = codecs.open(ntlcase, 'r', 'utf-8').readlines()[1:]  # first line is header: skip
    # puts all values into a struct like nc.well.lateral_1.segments.aicd_a.val
    nc = Container()
    for line in lines:
        line = line.strip()
        fullnm, val = line.split('=', 1)
        varnames = [_fix_varnm(x) for x in fullnm.split('.')]
        obj = nc
        for varnm in varnames:
            if varnm not in obj.__dict__:
                obj.__dict__[varnm] = Container()
            obj = obj.__dict__[varnm]
        obj.val = val.strip()
    return nc
