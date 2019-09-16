def dp_dupuit_thiem(q, mu, k, Rw, R):
   return q*mu*log(R/float(Rw)) / (2*pi*k) * U.ATM / U.BAR

# input
l_well = 1600            # length of injection well. m
Q0     = 1600.           # typical injection rate. m3/d for one well 
k      = 4.              # reservoir permeablity. Darcy
mu     = 10.             # fluid viscosity. cP
Rw     = 8.5*U.INCH / 2. # wellbore radius. m

# want to plot dp as function of Q
Q = Q0*linspace(0.1, 2, 2)
q = Q/U.CM**3 / U.DAY / (l_well/U.CM)  # in 'darcy units', cm2/s

figure()
for R in [1,10,100,1000]:
   dp = dp_dupuit_thiem(q, mu, k, Rw, R)
   plot(Q, dp, label='R=%im' % R)

plot([Q0,Q0], [0,dp[-1]], 'k--')

legend(loc='best')
xlabel('injection rate [m3/d]')
ylabel('injection pressure [bar]')
grid(True)
title('%iD permeability, %icP water, %im long well' % (k, mu, l_well))
show()

