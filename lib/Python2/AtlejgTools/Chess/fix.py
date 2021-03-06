def _create_pgnfile(gm, matchno):
   fnm = '%s_%02i.pgn' % (UT.basename(pgnfile), matchno)
   f = open(fnm, 'w')
   for gm_ in gm:
      f.write(gm_)
   f.close()
   print fnm, 'has been created'

pgnfile = 'carlsen_vs_anand_2013_wch.pgn'
lines = open(pgnfile).readlines()
gm = []
matchno = 0
for line in l:
   if 'Event ' in line:
      if gm:
         _create_pgnfile(gm, matchno)
      gm = [line]
      matchno += 1
      continue
   gm.append(line)
if gm:
   _create_pgnfile(gm, matchno)
