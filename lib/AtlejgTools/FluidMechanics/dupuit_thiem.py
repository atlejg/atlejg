import AtlejgTools.SimulationTools.UnitConversion as U
import AtlejgTools.FluidMechanics.PorousMedia as PM

# input
l_well = 1000            # length of injection well. m
Q0     = 1000.           # typical injection rate. m3/d for one well 
k      = 5.              # reservoir permeablity. Darcy
mu     = 0.5             # fluid viscosity. cP
Rw     = 8.5*U.INCH / 2. # wellbore radius. m

# want to plot dp as function of Q
q = Q0*linspace(0.1, 2, 2)

figure()
for R in [1,10,100,1000]: # varying reservoir radius
   dp = PM.dp_radialflow(q, mu, k, Rw, l_well, R)
   plot(q, dp, label='R=%im' % R)

plot([Q0,Q0], [0,dp[-1]], 'k--')

legend(loc='best')
xlabel('injection rate [m3/d]')
ylabel('injection pressure [bar]')
grid(True)
title('%iD permeability, %icP water, %im long well' % (k, mu, l_well))
show()

