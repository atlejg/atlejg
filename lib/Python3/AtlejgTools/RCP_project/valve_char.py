q = linspace(0,2.5, 100)
p=VC.rcp_params('ar3_grane')
rhow = 1060.
rhoo = 865.
muw = 0.6
muo = 12.
wc = linspace(0,1, 6)
rho = rhow*wc + rhoo*(1-wc)
mu = muw*wc + muo*(1-wc)
figure()
r = linspace(0,1,len(wc))
#[plot(q, VC.rcp_dp3(p, rho[i], mu[i], 1000*q), color=(r[i],0,1-r[i]), label='wc=%.1f'%wc[i]) for i in range(len(wc))]
[plot(q, VC.rcp_dp3(p, rho[i], mu[i], 1000*q), label='wc=%.1f'%wc[i]) for i in range(len(wc))]
ylim(0,40)
xlabel('dp [bar]')
ylabel('dp [bar]')
xlabel('flow rate [m3/h]')
grid(True)
legend(loc='upper left')
show()
