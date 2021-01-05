from scipy.signal import medfilt
import pylab as pl
import h5py
import pandas as pd
from datetime import datetime, date
from scipy.interpolate import interp1d
import AtlejgTools.SimulationTools.WellData as WellData
import AtlejgTools.Utils                    as UT

PPM_M        = 300e3 # PPM_M: ppm in mother-solution
MAX_VISC     = 50.
MAX_INJ_RATE = 1250.
MONTH_MAP    = {'JAN':1, 'FEB':2, 'MAR':3, 'APR':4, 'MAY':5, 'JUN':6, 'JUL':7, 'AUG':8, 'SEP':9, 'OCT':10, 'NOV':11, 'DEC':12} # useful when converting dates
TR_LABELS    = {'T-141e':'WT1', 'T-144c':'WT0', 'T-146m':'WT2'}
t0_excel = pl.date2num(datetime(1899, 12, 30))

def date2num1(date):
    '''
    DATE            WWCTH A-18
    TIME            FRACTION
    2/1/2014 0.005
    2/2/2014 0.006
    '''
    m, d, y = [int(x) for x in date.split('/')]
    return pl.date2num(datetime(y, m, d))

def date2num2(date):
    '''
    DATE            WWCTH A-18
    TIME            FRACTION
    15-Aug-2014 0.005
    16-Aug-2014 0.006
    '''
    d, m, y = date.split('-')
    m = MONTH_MAP[m.upper()]
    return pl.date2num(datetime(2000+int(y), m, int(d)))

def date2num3(date):
    '''
    handle date format like this: 05.05.2015 
    '''
    d, m, y = date.split('.')
    return pl.date2num(datetime(int(y), int(m), int(d)))

def read_txt_file1(fnm, skiprows, date2num_func=date2num1, delim=None):
    '''
    reads input data in text format that is typical for the peregrino pilot - like this:
    '''
    d   = pl.loadtxt(fnm, skiprows=skiprows, delimiter=delim, converters={0:date2num_func}, encoding='latin1')
    t   = d[:,0]
    ixs = pl.argsort(t)                                                 # make sure time is increasing. we do that by sorting
    return t[ixs], d[ixs,1]

def get_tracer_file(tracernm, wellnm, templ_fnm):
    if   tracernm == 'WT0':
        tracerfile = templ_fnm % ('114c', wellnm)
    elif tracernm == 'WT1':
        tracerfile = templ_fnm % ('141e', wellnm)
    elif tracernm == 'WT2':
        tracerfile = templ_fnm % ('146m', wellnm)
    return tracerfile

def get_shutins(datadir='/project/peregrino/users/agy/InputData/'):
    '''
    gets shutins from files provided by yngve.
    filenames are like:
    A-11_SHUT_2020.txt  A-11_TARs.txt  A-12_SHUT_2020_JFM.txt  A-18_SHUT_2020_JFM.txt 
    A-18_TARs.txt  A-22_SHUT_2020.txt  A-22_TARs.txt  A-27_SHUT_2020_JFM.txt
    '''
    shutins = {}
    for well in ['A-11', 'A-22', 'A-12', 'A-18', 'A-27']:
        dates, bhp = [], []
        for fn in UT.glob('%s/%s_*.txt'%(datadir,well)):   # could be more than one file per well
            d = read_txt_file1(fn, 1, date2num_func=date2num3, delim=';')
            dates.extend(d[0])
            bhp.extend(d[1])
        shutin = UT.Struct()
        shutin.nm = well
        dates = pl.array(dates)
        bhp = pl.array(bhp)
        # make sure data is sorted according to dates
        ixs = pl.argsort(dates)
        shutin.dates = dates[ixs]
        shutin.bhp = bhp[ixs]
        shutins[well] = shutin
    return shutins

def get_tracer(tracernm, wellnm, templ_fnm='/project/peregrino/users/agy/InputData/Tracers/T-%s_%s.txt', date2num_func=date2num1, skiprows=2):
    '''
    returns times, concentrations and unit of tracer
    '''
    fnm = get_tracer_file(tracernm, wellnm, templ_fnm)
    unit = UT.grep_column(fnm, 'TIME', 2, is_float=False)[0].lower()
    d = read_txt_file1(fnm, skiprows, date2num_func)
    return d[0], d[1], unit

def _read_tracer_data(excelfnm, pwell, iwell):
    #
    sheet    = 'Tracer Analysis Producer'
    #
    wellnms    = pl.array(UT.read_excel_column(sheet, 'A', 2, 99999, excelfnm))
    ixs        = (wellnms == pwell)
    dates      = pl.array(UT.read_excel_column(sheet, 'B', 2, 99999, excelfnm))[ixs]
    tracernms  = pl.array(UT.read_excel_column(sheet, 'C', 2, 99999, excelfnm))[ixs]
    conc       = pl.array(UT.read_excel_column(sheet, 'D', 2, 99999, excelfnm))[ixs]
    ts         = pl.array([float(x) + t0_excel for x in dates])
    #
    sheet    = 'Tracer Injection'
    #
    wellnms    = pl.array(UT.read_excel_column(sheet, 'A', 2, 99999, excelfnm))
    ixs        = (wellnms == iwell)
    dates      = pl.array(UT.read_excel_column(sheet, 'B', 2, 99999, excelfnm))[ixs]
    tracernms_ = pl.array(UT.read_excel_column(sheet, 'C', 2, 99999, excelfnm))[ixs]
    mass       = pl.array(UT.read_excel_column(sheet, 'D', 2, 99999, excelfnm))[ixs]
    #
    inj = {}
    for tracernm in pl.unique(tracernms_):
        ix = (tracernms_ == tracernm)
        t  = float(dates[ix]) + t0_excel
        inj[tracernm] = (t, float(mass[ix][0]))
    return tracernms, ts, conc, inj

def get_tracers(excelfnm='/project/peregrino/users/agy/InputData/Tracers/Tracer DB.xlsx', pwell='A-22', iwell='A-11'):
    '''
    read tracers from gulnar's excel sheet
    '''
    tracernms, ts, conc, inj = _read_tracer_data(excelfnm, pwell, iwell)
    #
    tracers = {}
    for tracernm in pl.unique(tracernms):
        ixs = (tracernms == tracernm)
        tr = UT.Struct()
        tr.nm = tracernm
        tr.label = TR_LABELS[tracernm]
        tr.t  = ts[ixs]
        tr.conc = conc[ixs] / 1000.
        tr.inj_t = inj[tracernm][0]
        tr.inj_m = inj[tracernm][1]
        tracers[tr.label] = tr
    return tracers

def get_bsw_wct(fnm, winsz=31, wc_func_only=True, skiprows=3, date2num_func=date2num1):
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
    t, wc = read_txt_file1(fnm, date2num_func, skiprows)
    wcf = medfilt(wc, winsz)
    wc_func = interp1d(t, wcf)
    if wc_func_only: return wc_func
    else           : return wc_func, t, wc, wcf

def visc_func_KE(ppm):
    return 0.7*(4e-6*ppm**2 - 0.0029*ppm + 0.6)   # excel trendline from kjetil E spread-sheet. scaled to match measured viscosities.

def visc_func_JFM(ppm):
    return 2.6e-6*ppm**2 + 0.4    # found by using fit_viscosity_function.py

def visc_func(ppm):
    return 0.9*2.6e-6*ppm**2 + 0.4    # based on visc_func_JFM, scaled to better match measured viscosities

def get_phi_F(t, c, q, m_inj, rho_inj):
    ''' 
    ref: Tracer interpretations... , Shook + Forsman
    calculates the storage capacity (Phi) and flow capacity (F) functions
    t must start at 0 [d]
    c has unit [g/m3]
    '''
    Et = c*rho_inj*q / m_inj
    Ett = Et*t
    int_Et = UT.cumulative(t, Et) 
    int_Ett = UT.cumulative(t, Ett)
    #   
    phi = int_Ett / int_Ett[-1]
    F = int_Et / int_Et[-1]
    return phi, F

def get_tracers(fnm='/project/peregrino/users/agy/InputData/Tracers/Tracer_analysis_Idaho_T144c_T141e_T146m_T145e_data-up-to-Mar15-2020_interpolation_JFMandDashboardRates.xls'):
    '''
    reads tracer-data (raw data and derived data) from yngve's excel
    '''
    detect_lvl = 0.01     # for defining when tracer break-through is happening
    inj_nms   = ['144c', '141e', '146m', '145e', '801']
    inj_dates = pl.array([pl.date2num(date(*x)) for x in \
                [(2014,7,30), (2017,1,14), (2017,11,26), (2019,1,3), (2019,7,25)]])
    #
    def _get_data(fnm, sheetnm, colnm, inj_nm):
        startln = 5 if sheetnm.startswith('Raw') else 10 
        y = UT.read_excel_column('%s_T%s'%(sheetnm,inj_nm), colnm, startln, 99999, fnm)
        y = [x if x else '0' for x in y]
        return pl.array([float(x) for x in y])
    #
    tracers = []
    for i, inj_nm in enumerate(inj_nms):
        s = UT.Struct()
        s.rdates = _get_data(fnm, 'Raw Input Data', 'B', inj_nm) + t0_excel   # raw dates
        s.rconcs = _get_data(fnm, 'Raw Input Data', 'C', inj_nm)              # raw concentration [pbb]
        s.dates = _get_data(fnm, 'Raw Input Data', 'G', inj_nm) + t0_excel    # interpol. dates
        s.concs = _get_data(fnm, 'Raw Input Data', 'H', inj_nm)               # concentration [g/sm3]
        # injection mass (scalar)
        s.inj_m = UT.read_excel_column('Input Data_T%s'%inj_nm, 'B', 2, 2, fnm)[0]
        # tracer break-through (scalar)
        ix = pl.where(s.rconcs > detect_lvl*max(s.rconcs))[0][0]
        s.breaktr = s.rdates[ix] - inj_dates[i]
        # remove data prior to tracer-start
        ixs = s.rdates >= inj_dates[i]
        s.rdates = s.rdates[ixs]
        s.rconcs = s.rconcs[ixs]
        # Ø-F etc.
        s.phis = _get_data(fnm, 'Results', 'M', inj_nm)
        s.fs   = _get_data(fnm, 'Results', 'N', inj_nm)
        s.tds = _get_data(fnm, 'Results', 'U', inj_nm)
        s.evs = _get_data(fnm, 'Results', 'V', inj_nm)
        s.name = inj_nm
        s.sim_name = 'WTPCWT%i:A-22'%i                    # useful when comparing to yngve's simulations
        s.label = pl.num2date(inj_dates[i]).strftime('%b:%-y')
        s.inj_date = inj_dates[i]
        tracers.append(s)
    return tracers

def get_a11_and_a22(dirname='/project/peregrino/users/agy/InputData/'):
    df = pd.read_csv('%s/A11.csv'%dirname, delimiter=';')
    a11 = UT.Struct()
    a22 = UT.Struct()
    a11.dates = df['DATE']
    a11.wir   = df['WWIRH']
    a11.bhp   = df['WBHPH']
    a11.cic   = df['WCIC']
    a11.cir   = df['WCIR']
    a11.cit   = df['WCIT']
    a11.time  = a11.dates - a11.dates[0]
    # also include start/stop of polymer-injection
    a11.pdates = [ (datetime(2018,1,6),  datetime(2018,4,15)),
                   (datetime(2019,1,13), datetime(2019,10,7)) ]
    # and some key shutins
    a11.shutins = [
        pl.date2num(datetime(2019,3,1,13,30)),
        pl.date2num(datetime(2019,5,9,14)),
        pl.date2num(datetime(2019,9,8,17)) ]
    a22.shutins = [pl.date2num(datetime(2020,4,8))]        # 2020 shutin of the field
    # and the ILT-dates
    a11.ilt_dates = [
        pl.date2num(datetime(2017,2,1)),
        pl.date2num(datetime(2019,7,28)) ]
    #
    # tars:
    a11.tars = [
        pl.date2num(datetime(2015,5,5)),
        pl.date2num(datetime(2016,4,22)),
        pl.date2num(datetime(2017,4,4)),
        pl.date2num(datetime(2018,3,29)),
        pl.date2num(datetime(2019,4,4)) ]
    a11.itars = [pl.where(a11.dates == x)[0][0] for x in a11.tars]   # useful indices
    a22.tars = a11.tars
    a22.itars = a11.itars
    df = pd.read_csv('%s/A22.csv'%dirname, delimiter=';')
    a22.dates = df['DATE']
    a22.opr   = df['WOPRH']
    a22.wpr   = df['WWPRH']
    a22.wct   = df['WWCTH']
    a22.bhp   = df['WBHPH']
    a22.time  = a22.dates - a22.dates[0]
    #
    # tracers
    a22.tracers = get_tracers()
    return a11, a22

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
    # q-skid is found in two ways
    #    skid-1: outflow
    tqs = h5['A-11']['WWIRS'][:,0]       # time - inj-rate series
    qs1  = h5['A-11']['WWIRS'][:,1] * 24  # m3/h -> m3/d
    #    resampling and summing. note: we WILL extrapolate if some tag is missing data in either direction
    t1 = min(tqc[0], tqi[0], tqm[0], tqd[0], tqs[0])
    t2 = max(tqc[-1], tqi[-1], tqm[-1], tqd[-1], tqs[-1])
    n  = max(len(tqc), len(tqi), len(tqm), len(tqd), len(tqs))
    dt = (t2-t1) / (n-1)
    t  = pl.linspace(t1,t2, n)
    print('last date', pl.num2date(t[-1]))
    qc = pl.interp(t, tqc, qc)
    qi = pl.interp(t, tqi, qi)
    qm = pl.interp(t, tqm, qm)
    qd = pl.interp(t, tqd, qd)
    qs1 = pl.interp(t, tqs, qs1)
    p  = pl.interp(t, tp, p)
    #    skid-2: adding inflows
    qs2 = qi + qd + qm
    #
    #
    # do some modifications since data is bad..
    #
    # must handle polymer-period 2 separately :-(
    ix = (t < pl.date2num(datetime(2019,2,13))).nonzero()[0][-1]
    qw = pl.concatenate((qc[:ix], qs1[ix:]))
    # also: remove highest values (not realistic)
    qw = pl.minimum(qw, MAX_INJ_RATE)
    # for some reason it has inj-rate before well is opened...
    ixs = t < pl.date2num(datetime(2014,6,1))
    qw[ixs] = 0.0
    # adjust inj-rate for first polymer-phase
    t1, t2 = [pl.date2num(datetime(2018,1,6)), pl.date2num(datetime(2018,4,15))]
    ixs = pl.logical_and(t>=t1, t<=t2)
    qw[ixs] = qs2[ixs]       # for the first polymer injection period, we must use qs2
    qw[ixs] *= 1.6           # and scale it
    #
    # create WellData object
    a11 = WellData.WellData('A11', welltype='inj', t=t, qw=qw, p=p, dt=dt)
    # add some properties
    ppm = PPM_M * qm / qs1
    ppm[ppm>PPM_M/60.]  = 0.
    ppm[ppm<0.]         = 0.
    visc = visc_func(ppm)
    visc[visc<0]        = pl.NaN
    visc[visc>MAX_VISC] = pl.NaN
    #    other useful stuff
    a11.ppm  = ppm
    a11.visc = visc
    a11.qi   = qi
    a11.qd   = qd
    a11.qm   = qm
    a11.qs1  = qs1
    a11.qs2  = qs2
    a11.qc   = qc
    #
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
    #
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
    typical usage:
    rs = concs / max(concs)
    [axvspan(ts[i], ts[i+1], facecolor='r', alpha=rs[i]) for i in range(len(concs)-1)]
    '''
    lines = open(fnm).readlines()
    dates, concs = [],[]
    for line in lines[3:]:
        if line.startswith('--'): continue
        r = line.strip().split()
        if len(r) == 1: concs.append(float(r[0]))
        else:           dates.append(datetime(int(r[2]), MONTH_MAP[r[1].replace("'","")], int(r[0])))
    return pl.array(dates), pl.array(concs)

