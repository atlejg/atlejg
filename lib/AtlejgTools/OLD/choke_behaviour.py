'''
plots some result vector (like FOPR) vs time toghether
with a graph showing number of chokes being active
'''

varnm     = sys.argv[1]
datafiles = sys.argv[2:]
datafiles.sort()

def _resample(t, n_choked, t2):
   if len(n_choked) == 1 or max(n_choked) == 0: return zeros(len(t2))
   n_choked2 = []
   i = 0
   for t0 in t2[:-1]:
      if i < len(t)-1 and t0 > t[i+1]: i = i+1
      n_choked2.append(n_choked[i])
   n_choked2.append(n_choked[i])
   return n_choked2

def cooling(casenm, segm, cpw=4200., cpo=2100.):
    s  = PP.get_summary(casenm)
    wc = s.get('SWCT-%s-%i'%(PP.WELLNM,segm))
    q  = s.get('SLFR-%s-%i'%(PP.WELLNM,segm))
    return q * (wc*cpw + (1-wc)*cpo)

fig = figure()
caseno = 0
for datafile in datafiles:
   caseid = UT.basename(datafile)
   caseno += 1
   n_choked = [0]
   t = [0]
   lineno = -1
   j = -Inf
   for line in open('%s.PRT'%caseid).readlines():
      lineno += 1
      if line.startswith(' @--MESSAGE  AT TIME'):
         j = lineno
         t0 = float(line.split()[3])
      if 'CHOKE_' in line and lineno == j+2:
         n_choked.append(n_choked[-1]+1)
         t.append(t0)
      if 'UNCHOKE' in line and lineno == j+1:
         n_choked.append(0)
         t.append(t0)
   if caseno == 1:
      ax1 = fig.add_subplot(111)
      ax2 = ax1.twinx()
   s = PP.get_summary(caseid)
   ax1.plot(s.time, _resample(t,n_choked,s.time), PP.COLOURS[caseno-1], marker='o', markersize=5)
   ax2.plot(s.time, s.get(varnm), PP.COLOURS[caseno-1], label=caseid)

ax1.set_xlabel('time [days]')
ax1.set_ylabel('# valves choked')
ax2.set_ylabel(varnm)
ax2.legend(loc='best')
show()
