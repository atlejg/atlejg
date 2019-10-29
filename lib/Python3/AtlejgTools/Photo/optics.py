f0 = sqrt(2)
# circle of confusion is dependent on camera (or at least the size of the bildebrikke)
# see f.ex. http://photo.tutsplus.com/tutorials/shooting/using-hyperfocal-distance-to-ensure-maximum-depth-of-field-in-landscape-photography/
_coc = 0.02e-03 # circle of confusion for eos500d

N = arange(0,10)
fstops = array([f0**n for n in N])

def hyperfocal_dist(f, fstop, coc):
    return f**2/coc/fstop

def calc_dof(f, fstop, s):
    '''
    f    : focal length [m]
    fstop: aperture, fstop in [1, 1.4, 2, 2.8 ...]
    s    : distance to object in m. could be vector
    '''
    H = hyperfocal_dist(f, fstop, _coc)
    print('Hyperfocal distance H = %.2fm' % H)
    Dn  = H*s / (H+s)    # near limit
    Df  = H*s / (H-s)    # far limit
    ii = (s > H).nonzero()[0]
    if len(ii) > 0: Df[ii[0]:] = Inf
    DoF = Df - Dn        # depth of field
    return (Dn, Df, DoF)

indx = 1
figure(1); clf()
figure(2); clf()
f = 24e-03          # focal length [m]
s = arange(1, 20)   # distance to object [m]
i = -1
DoFs = []
Dns  = []
Dfs  = []
for fstop in fstops:
    i += 1
    red = i/float(len(fstops))
    Dn, Df, DoF = calc_dof(f, fstop, s)
    figure(1)
    plot(s, Dn, '*-',   color=(red,0,1-red), label='f/%.1f'%fstop)
    plot(s, Df, '*--',  color=(red,0,1-red))
    figure(2)
    plot(s, DoF, '*-', color=(red,0,1-red), label='f/%.1f'%fstop)
    DoFs.append(DoF[indx])
    Dns.append(Dn[indx])
    Dfs.append(Df[indx])
figure(1)
xlabel('subject distance [m]')
ylabel('distance [m]')
grid(True)
title('whole line: near limit. dashed line: far limit')
legend(loc='best')
ylim(0,40)
figure(2)
xlabel('subject distance [m]')
ylabel('distance [m]')
grid(True)
title('DoF')
legend(loc='best')
ylim(0,40)
figure(3); clf()
plot(fstops, DoFs, label='DoF')
plot(fstops, Dns,  label='near lim')
plot(fstops, Dfs,  label='far lim')
xlabel('aperture')
ylabel('distance [m]')
grid(True)
title('DoF, distance = %im' % s[indx])
#ylim(0,40)
legend(loc='best')
show()
