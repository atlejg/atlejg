#!/usr/bin/env python
'''
atle j. gyllensten
agy@statoil.com

Reads a dbg-file which has been created using the WELDEBUG <wellname> 3000 keyword in Eclipse.
Creates a plain ascii file with a header describing the content.

Note! It assumes only one well is being 'debuged'.

Usage   : <this file> <case name> segment
Example : <this file> RCPCASE 34
Make sure the segment is where a inflow-valve is found
'''

import sys

casenm = sys.argv[1]
segm_no = int(sys.argv[2])

potential_solving = False

# create ascii
fname = '%s_segm=%i.txt' % (casenm, segm_no)
print 'creating file', fname
fout = open(fname, 'w')
fout.write('# time density viscosity flowrate dp\n')

# read dbg-file
fin = open(casenm+'.DBG')
for line in fin :
   if ' 0--TIME STEP STARTING' in line:
      rec = line.split()
      time = float(rec[7])
      potential_solving = False
   if 'CALCULATING WELL POTENTIAL' in line:
      potential_solving = True
   if potential_solving : continue
   # array filled redundantly (for each iteration), only last one is interesting
   if 'WCFRSG' in line:
      rec = line.split()
      if not int(rec[1])-segm_no ==0 : continue
      fout.write('%s %s %s %s %s\n' % (time, rec[4], rec[3], rec[2], rec[7]))
fin.close()
fout.close()
