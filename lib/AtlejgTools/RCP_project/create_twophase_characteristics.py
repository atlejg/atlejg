def wcf(wc, wc2): return UT.elbow_func((0,0), (0.8,wc2), (1,1), wc)

templfile = sys.argv[1]  # typically A.TEMPL

flow_exponents = (1, 2, 3.3)
elbowbends     = (0.1, 0.4, 0.8, -1)  # watercut where elbow has it bend. values < 0 are special
q              = linspace(0,20,11)    # m3/d/segment
dp_maxes       = 140*linspace(0.25, 4, 9) # this is the 'strength' of the valve (ala VPJ)
watercuts      = linspace(0,1,11)

i = 0
dr = 1./(len(flow_exponents)*len(elbowbends)) # just for coloring
for x in flow_exponents:
   i += 1
   j = 0
   for wc2 in elbowbends:
      j += 1
      if wc2 < 0:
         _wcf = lambda wc: 1.         # no elbow, just a constant 1.
      else:
         _wcf = lambda wc: wcf(wc, wc2)
      k = 0
      for dp_max in dp_maxes: 
         k += 1
         dp, txt = ECL.create_vfp_table2(q, lambda q: q**x, watercuts, _wcf, dp_max, 2355)
         repl = (('_VFPPROD_', txt), )
         datafile = '%s%i%i%i.DATA' % (UT.basename(templfile), i,j,k)
         print 'creating', datafile
         UT.replace_in_file(repl, templfile, datafile)
         if k == 3: # dont plot everything
            figure()
            for m in range(len(dp)):
               r = m*dr
               plot(dp[m], color=(r,0,1-r), label='wc=%.2f'%watercuts[m])
            title('%i,%i'%(i,j))
            xlabel('flow rate')
            ylabel('dp [bar]')
            grid(True)
            legend(loc='best')

show()
