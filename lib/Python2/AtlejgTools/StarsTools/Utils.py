import sys
import AtlejgTools.Utils as UT
import AtlejgTools.EclipseTools.Utils as ECL
import pylab as pl

def _val(line):
   return float(line.split('+')[1].strip())

def read_outfile(fnm):
   fnm = UT.basename(fnm) + '.out'
   i = -pl.Inf
   t,bhp,dd,o,g,w,gor,wc = [],[],[],[],[],[],[],[]
   skip = True
   for line in open(fnm):
      if line.startswith(' Well parameters')              : skip = False
      if line.startswith(' Surface cumulative production'): skip = True
      if skip: continue
      if 'Cumulative Open days' in line: t.append(_val(line))
      if 'Bottom Hole'          in line: bhp.append(_val(line))
      if 'Drawdown'             in line: dd.append(_val(line))
      if 'Oil'                  in line: o.append(_val(line))
      if 'Gas'                  in line: g.append(_val(line))
      if 'Water'                in line: w.append(_val(line))
      if 'GOR'                  in line: gor.append(_val(line))
      if 'WCUT'                 in line: wc.append(_val(line))
   t,bhp,dd,o,g,w,gor,wc = pl.array(t),pl.array(bhp),pl.array(dd),pl.array(o),pl.array(g),pl.array(w),pl.array(gor),pl.array(wc)
   m = pl.array([t,o,g,w,o+w,o+w+g,wc,gor,bhp/100.,dd/100.])
   varnms = ['TIME','FOPR','FGPR','FWPR','FLPR','FTOT','FWCT','FGOR','BHP', 'DD']
   units  = ['days','m3/d','m3/d','m3/d','m3/d','m3/d','-',   '-',   'Bar', 'Bar']
   return ECL.TimeSeries(m.T, varnms, units, UT.basename(fnm))

def read_logfile(fnm):
   '''
   note: logfile dosnt have values per well. in particular - no BHP
   '''
   i = -pl.Inf
   m = []
   for line in open(fnm):
      if 'Stopping end time reached' in line: break
      if line.startswith(' ---Time Step'): i = -5
      i += 1
      if 1 <= i <= 20 :
         rec = line.split()
         if len(rec) != 14: continue
         t = float(rec[3])
         o = float(rec[5])
         g = float(rec[6])
         w = float(rec[7])
         l = o + w
         q = l + g
         wc = w / l
         gvf = g / q
         m.append(pl.r_[t,o,g,w,l,q,wc,gvf])
      if i == 20: i = -pl.Inf
   m = pl.array(m)
   varnms = ['TIME','FOPR','FGPR','FWPR','FLPR','FTOT','FWCT','FGVF']
   units  = ['days','m3/d','m3/d','m3/d','m3/d','m3/d','-',   '-']
   return ECL.TimeSeries(m, varnms, units, UT.basename(fnm))

if '__main__' in __name__:
   varnm = sys.argv[1]
   outfiles = UT.glob(sys.argv[2:])
   #
   figure()
   for outfile in outfiles:
      s = read_outfile(outfile)
      pl.plot(s.time, s.get(varnm), label=s.nm)
   #
   pl.xlabel('time [days]')
   pl.ylabel('%s [%s]' % (varnm, s.unit(varnm)))
   pl.legend(loc='best')
   pl.grid(True)
   #
   pl.show()

