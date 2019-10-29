# north komsomolskoja for ali.
rho_o = 920.                 # kg/m3
rho_g = 87.5                 # kg/m3
mu_o  = 115.                 # cp
mu_g  = 0.012                # cp
vp = VC.rcp_params('tr7_peregrino')

a = float(sys.argv[1])
b = float(sys.argv[2])
shift = float(sys.argv[3])
c2    = float(sys.argv[4])

COLOUR = 'gbkmcr'

rho_cal, mu_cal, cnst, x, y = vp

# plot valve characteristics for varying gvf's
gvf = r_[0, 0.25, 0.5, 0.75, 0.9]
gvf = linspace(0.,1, 350)
q = linspace(1, 2000, 100)
rm = VC.rho_mix(rho_o,rho_g, gvf)
vm = VC.mu_mix(mu_o,mu_g, gvf)
f1 =  (mu_cal/vm)**y * cnst
f2 =  (a-vm/mu_cal) * b
f3 =  (mu_cal/(vm+shift))**y * c2
fig = figure()
plot(gvf, f1, label='std')
plot(gvf, f2, label='linear')
plot(gvf, f3, label='std, shifted')
legend(loc='best')
xlabel('gvf')
show()
