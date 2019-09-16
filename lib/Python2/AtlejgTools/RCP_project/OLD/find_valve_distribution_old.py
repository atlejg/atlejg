'''
Tries to find number of valves in gas-zone & oil-zone so that the total flow-rate of gas and oil is honored.
Typically used for a Troll-well.
It will look through all combinations of valves and pick the one that gives the smallest difference in dp.

 Zone 1 is the gas-zone
 Zone 2 is the oil-zone
 _g : gas
 _o : oil
 Input is an input-file that looks like this (see e.g. I13_AYH_Y1.inp):

# in-situ flowrates
float qo    = 1200.  # I13 AYH Y1 (oyvind m.) (really q-liquid, but assume its oil...)
float qg    = 1650.  # I13 AYH Y1 (oyvind m.)
# fluid props
float rho_g = 116.
float rho_o = 805.
float mu_g  = 0.018
float mu_o  = 1.8
# number of valves in total
int   nv    = 231
# try these gvf's
eval  gvf1  = [1, 0.8, 0.6]  # gas fraction in gas-zone
eval  gvf2  = [0, 0.2, 0.4]  # gas fraction in oil-zone

'''

from scipy.optimize import fsolve
import sys

iv = UT.InputValues(sys.argv[1])

'''
old stuff. tried to solve it very generally, but ran into problems (matrix inversion , division by zero etc)
def calc_dp_old(gvf1, gvf2, nv1):
   ## solving eq-system
   ##A = matrix(([gvf1,gvf2], [1-gvf1, 1-gvf2]))
   ##b = [qg, qo]
   ##q1, q2 = solve(A,b)      # flow rate in gas-zone + oil-zone, respectively
   q1 = q2 = 1
   if (0 < gvf1 < 1) and (0 < gvf2 < 1):
      q2 = (qo-qg*(1-gvf1)/gvf2)  /  ((1-gvf2) - gvf2/gvf1*(1-gvf1))
      q1 = (qg-gvf2*q2) / gvf1
   elif gvf1 == 0 and gvf2 > 0:
      q2 = qg/gvf2
      q1 = qo -q2*(1-gvf2)
   elif gvf1 == 1 and gvf2 < 1:
      q2 = qo/(1-gvf2)
   elif gvf2 == 0 and gvf1 > 0:
      q1 = qg/gvf1
      q2 = qo - q1*(1-gvf1)
   elif gvf2 == 1 and gvf1 < 1:
      q1 = qo/(q-gvf1)
      q2 = qg - gvf1*q1
   else:
      print "not covered"
   #print gvf1, gvf2, q1, q2
   # mix-props
   rho1 = gvf1*rho_g + (1-gvf1)*rho_o
   rho2 = gvf2*rho_g + (1-gvf2)*rho_o
   mu1  = gvf1*mu_g  + (1-gvf1)*mu_o
   mu2  = gvf2*mu_g  + (1-gvf2)*mu_o
   # try to match dp-valve
   nv2 = nv - nv1
   #print nv, nv1, nv2
   if nv1*nv2 == 0:
      print 'nv1 or nv2 is 0', nv1, nv2
      return 0, 0
   qv1 = q1/float(nv1)  # flow per valve
   qv2 = q2/float(nv2)  # flow per valve
   #print gvf1, nv1, nv2, qv1, qv2, mu1, mu2, rho1, rho2
   return dpv(rho1,mu1,qv1), dpv(rho2,mu2,qv2)

def func(x):
   gvf2 = 0.                # gas-fraction for "oil-valves". hardcoded
   gvf1, nv1 = x            # gas-fraction for "gas-valves" + number of valves in gas-zone
   dp = calc_dp(gvf1, gvf2, nv1)
   err1 = (dp1 - dp[0])**2
   err2 = (dp2 - dp[1])**2
   return (err1, err2)

x0 = (1., int(nv/2))
f = fsolve(func, x0, full_output=True)
print f[-1]
gvf1, nv2 = f[0]
'''

def dpv(rho,mu,q):
   q *= 1000/24.  # m3/d -> l/h
   dp = VC.rcp_dp2('tr7_troll', rho, mu, q)
   return dp

def calc_dp(gvf1, gvf2, nv1):
   # gvf1 & gvf2 cant both be 0 or 1.
   if gvf1+gvf2 in (0,2): return -99999, 99999
   # gas fraction in gas-zone *must* be higher than in the oil zone
   if gvf1 <= gvf2:       return -99999, 99999
   #
   if gvf1 != 0:
      q2 = (iv.qo-iv.qg*(1-gvf1)/gvf1) / ((1-gvf2) - gvf2/gvf1*(1-gvf1))
      q1 = (qg-gvf2*q2) / gvf1
   else:
      q2 = iv.qg/gvf2
      q1 = iv.qo -q2*(1-gvf2)
   # we could get bad answers
   if q1 <= 0 or q2 <= 0: return -99999, 99999
   # mix-props
   rho1 = gvf1*iv.rho_g + (1-gvf1)*iv.rho_o
   rho2 = gvf2*iv.rho_g + (1-gvf2)*iv.rho_o
   mu1  = gvf1*iv.mu_g  + (1-gvf1)*iv.mu_o
   mu2  = gvf2*iv.mu_g  + (1-gvf2)*iv.mu_o
   # calculate dp-valve for each zone
   nv2 = iv.nv - nv1
   qv1 = q1/float(nv1)  # flow per valve
   qv2 = q2/float(nv2)  # flow per valve
   return dpv(rho1,mu1,qv1), dpv(rho2,mu2,qv2)

for gvf1 in iv.gvf1:
   for gvf2 in iv.gvf2:
      dpmin = Inf
      indx  = -1
      for indx in range(2,iv.nv-1):
         dp = calc_dp(gvf1, gvf2, indx)
         ddp = abs(dp[0]-dp[1])
         if ddp < dpmin:
            dpmin = ddp
            nv1 = indx
            nv2 = iv.nv - nv1
      print gvf1, gvf2, nv1, nv2, '%.2f'%dpmin, '%.2f'%calc_dp(gvf1, gvf2, nv1)[0]

