# atle. j. gyllensten
# R1: Frank M. White
# R2: fluent doc: http://www.fluentusers.com/fluent/doc/ori/v121/fluent/fluent12.1/help/html/ug/node233.htm

from math import pi

class Struct: pass

def inertial_coeff_for_pipeflow(diam, rho, mu, q_ref):
   '''
   all input in si-units.
   q_ref: typical flowrate for the given problem. if it varies along the pipe,
   result will be incorrect.
   '''
   r = diam / 2.
   A = pi*r**2
   v = q_ref / A
   re = rho*v*diam/mu; # reynolds number
   s = Struct()
   if re < 2300.:
      s.turbulent = False
      s.C1 = 32. / diam**2 # see R1 p. 342 and R2
      s.C2 = 0
   else:
      f = 0.316*re**-0.25 # blasius. R1 p. 345 and R2
      s.turbulent = True
      s.C1 = 0
      s.C2 = f
   return s
