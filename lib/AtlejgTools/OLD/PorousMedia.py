import numpy
import AtlejgTools.RCP_project.Valve_characteristics as VC

_tiny = numpy.finfo(float).eps
JOINT_LENGTH = 12.                                               # length of joint [m]

def darcy_dp(flowrate, perm, visc, area, length): 
   ''' assuming all units are si '''
   return flowrate*visc*length / (perm*area)

def relperm(sw, wirr, nwcs=None):
   ''' sw : saturation of wetting phase
       wirr: irreducible wetting phase satuation ("connate water")
       nwcs: nonwetting critical saturation
       if nwcs is given, it is assumed that relperm for nonwetting phase is asked for,
       else relperm for wetting phase is assumed.
       simple linear function.
   '''
   if nwcs == None:
      # this is wetting phase
      if sw <= wirr: return _tiny
      else        : return (sw - wirr) / (1. - wirr)
   else:
      # non-wetting phase
      if   sw <= wirr: return 1.
      elif sw >= nwcs: return _tiny
      else           : return 1. - (sw - wirr) / (nwcs - wirr)

class Fluid(object):
   def __init__(self, mu, rho, nm):
      self.mu  = mu
      self.rho = rho
      self.nm  = nm

class ReservoirZone(object):
   '''
   NB! only saturations of 0 and 1 is handled. *very* simple relperm-effect included.
   flow_resistance is sort of a scaled permeability.
   '''
   fluid1 = None                                            # must be set
   fluid2 = None                                            # must be set
   def __init__(self, s1, q, vpj, vt):
      self.s1        = s1                                   # saturation of fluid1
      self.q         = q
      self.vpj       = vpj
      self.vt        = vt                                   # valve type
      self.PI_scaler = 1.                                   # > 1 if hi-perm zone
      self.kr_end    = 0.35                                 # relperm @ s1 == 0
   def relperm(self): # not very intelligent
      return 1.0 if self.s1 > 0 else self.kr_end
   def viscosity(self):
      return ReservoirZone.fluid2.mu*(1-self.s1) + ReservoirZone.fluid1.mu*self.s1
   def density(self):
      return ReservoirZone.fluid2.rho*(1-self.s1) + ReservoirZone.fluid1.rho*self.s1
   def PI(self):
      return 1./self.flow_resistance/self.viscosity() * self.relperm() * self.PI_scaler
   def dp_reservoir(self):
      return self.l * self.q/self.PI()
   def nvalves(self):
      return self.l/JOINT_LENGTH * self.vpj
   def dp_valve(self):
      '''
      for a given flow, calculates pressure drop for a number of given valves.
      [q] = m3/d (will be converted to l/h for RCP-valve)
      '''
      if self.vt[:2] in ('ar', 'tr', 'bh', 'ha'):
         return VC.rcp_dp2(self.vt, self.density(), self.viscosity(), self.q*1000/24./self.nvalves())
      elif 'mm' in self.vt: # icd nozzle with fixed diameter - like 3.0mm
         diam = float(self.vt.replace('mm',''))
         return VC.nozzle_dp(self.q, self.density(), diam, self.nvalves())
      else: return -1.  # unphysical...
   def dp(self):
      return self.dp_reservoir() + self.dp_valve()
   def flow_pr_valve(self):
      qv = self.q / self.nvalves() # m3/d
      return qv*1000/24.           # l/h
   def __str__(self):
      return '''
 length          = %i m
 flowrate        = %i m3/d
 saturation      = %.2f (%s)
 flow_resistance = %.2e
 dp              = %.2f bar
 dp_res          = %.2f bar
 dp_valve        = %.2f bar
 valve type      = %s
 vpj             = %.2f
 flow pr valve   = %.1f l/h
 ''' % (round(self.l), self.q, self.s1, Fluid1.nm, self.flow_resistance, self.dp(), self.dp_reservoir(), self.dp_valve(), self.vt, self.vpj, self.flow_pr_valve())



if __name__ == '__main__':
   import AtlejgTools.SimulationTools.UnitConversion as u
   A = 0.01
   L = 0.1
   v = 0.001
   q = A*v
   wirr = 0.2
   nwcs = 0.85
   sw = 1. - r_[0.1 : 1 : 0.1]
   perm = 4 * u.DARCY
   visc = 1 * u.CP
   #
   rp = []
   for sw_ in sw:
      rp.append(relperm(sw_, wirr))
   dp = darcy_dp(q, perm, visc, A, L) / array(rp)
