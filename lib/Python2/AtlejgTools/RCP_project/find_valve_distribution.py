'''
Atle J. Gyllensten
agy@statoil.com
Feb 2016

Tries to find number of valves in gas-zone & oil-zone so that the total flow-rate of gas and oil is honored.
Typically used for a Troll-well.
(This is an expansion of the method proposed by Sigurd (and Martin Halvorsen) to estimate how many
valves see gas and how many see oil).
It will look through all combinations of valves and pick the best solution.
If input value 'mode' is set to fixed_dp, it will try to match the given input-dp, if not, it will try to match
the dp in gas-zone to the dp in the oil-zone

(I tried to solve it very generally, but ran into problems (matrix inversion , division by zero etc), so
i just do it by brute force now...)

 Zone 1 is the gas-zone
 Zone 2 is the oil-zone
 *_g : gas-property
 *_o : oil-property

 Input is an input-file that looks like this (see e.g. Data/I13_AYH_Y1_drawdown=15.inp):

# in-situ flowrates [m3/d]
float qo    = 1550.  # I13 AYH Y1 (oyvind m.) (really q-liquid, but assume its oil...)
float qg    = 1850.  # I13 AYH Y1 (oyvind m.)
# the assumed dp's
string mode = fixed_dp  # fixed_dp or find_dp
float dp1   = 15.       # only used for fixed_dp mode
float dp2   = 8.       # only used for fixed_dp mode
# fluid props
float rho_g = 116.
float rho_o = 805.
float mu_g  = 0.018
float mu_o  = 1.8
# inflow valves
int   nv    = 231            # total number of valves
string valve_type = tr7_troll
# try these gvf's
eval  gvf1  = [1, 0.8, 0.6]  # gas fraction in gas-zone
eval  gvf2  = [0, 0.2, 0.4]  # gas fraction in oil-zone
#eval  gvf1  = [1]
#eval  gvf2  = [0]
# extra info?
bool debug = False

'''

from scipy.optimize import fsolve
import sys

iv = UT.InputValues(sys.argv[1])

def _transition_func(gvf,b,v):
   '''
   juan's transition function - to enhance WSEGAICD function
   '''
   return (1-exp(-b*(1-gvf)))**v * (1+exp(-b))**v / ((1+exp(-b*(1-gvf)))**v * (1-exp(-b))**v)

def dpv2(gvf, q, b=0.5, v=0.2):
   '''
   valve pressure drop.
   using juan's transition function
   '''
   rho_cal, mu_cal, cnst, x, y = VC.rcp_params(iv.valve_type)
   f = _transition_func   # for convinience
   denom = f(gvf,b,v)*(rho_cal/(cnst*iv.rho_o**2)*(iv.mu_o/mu_cal)**y)**(1./x) + (1-f(gvf,b,v))*(rho_cal/(cnst*iv.rho_g**2)*(iv.mu_g/mu_cal)**y)**(1./x) 
   return denom**(-x) * q**x 

def dpv(gvf, q):
   '''
   valve pressure drop.
   NOT using juan's transition function
   '''
   # mix-props
   rho = gvf*iv.rho_g + (1-gvf)*iv.rho_o
   mu  = gvf*iv.mu_g  + (1-gvf)*iv.mu_o
   q *= 1000/24.  # m3/d -> l/h
   dp = VC.rcp_dp2(iv.valve_type, rho, mu, q)
   return dp

def calc_dp(gvf1, gvf2, nv1):
   # gvf1 & gvf2 cant both be 0 or 1.
   if gvf1+gvf2 in (0,2): return -99999, 99999
   # gas fraction in gas-zone *must* be higher than in the oil zone
   if gvf1 <= gvf2:       return -99999, 99999
   #
   if gvf1 != 0:
      q2 = (iv.qo-iv.qg*(1-gvf1)/gvf1) / ((1-gvf2) - gvf2/gvf1*(1-gvf1))
      q1 = (iv.qg-gvf2*q2) / gvf1
   else:
      q2 = iv.qg/gvf2
      q1 = iv.qo -q2*(1-gvf2)
   # we could get bad answers
   if q1 <= 0 or q2 <= 0: return -99999, 99999
   # calculate dp-valve for each zone
   nv2 = iv.nv - nv1
   qv1 = q1/float(nv1)  # flow per valve
   qv2 = q2/float(nv2)  # flow per valve
   return dpv2(gvf1,qv1), dpv2(gvf2,qv2)

results = []
print 'gvf1 gvf2 nv1 nv2 dp1  dp2    (mode= %s)' % iv.mode
for gvf1 in iv.gvf1:
   for gvf2 in iv.gvf2:
      r = UT.Struct()  # result
      r.err = Inf
      for nv1 in range(1,iv.nv):
         dp = calc_dp(gvf1, gvf2, nv1)
         if iv.mode == 'fixed_dp':
            err = norm([dp[0]-iv.dp1, dp[1]-iv.dp2])  # try to match given dp1, dp2
         else:
            err = abs(dp[0]-dp[1])                    # just try to make them equal
         if err < r.err:
            r.err = err
            r.nv1 = nv1
            r.nv2 = iv.nv - r.nv1
            r.dp1, r.dp2 = dp
         if iv.debug: print 'DEBUG %.2f %.2f %03i %03i %.2f %.2f' % (gvf1, gvf2, r.nv1, r.nv2, dp[0], dp[1])
      print '%.2f %.2f %03i %03i %.2f %.2f' % (gvf1, gvf2, r.nv1, r.nv2, r.dp1, r.dp2)
      r.gvf1 = gvf1
      r.gvf2 = gvf2
      results.append(r)

