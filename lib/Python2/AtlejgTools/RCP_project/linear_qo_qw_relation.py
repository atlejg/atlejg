'''
trying to check implications of sigurd's statement
that there is a linear relation between oil-rate and water-rate
when running with varying wct at a constant dP.
i.e. qw = qw0 - k*qo
'''
# user input
qw0 = 500.      # flow rate water @ dp_max
k = 1/3.        # water/oil ratio
dp_max = 30.    # dp @ qw0

# calculatins
qos = linspace(0,qw0/k, 101)
wcts = linspace(0,1, 11)

def qw_func(qo):  return qw0-k*qo
def wc_func(qo):  return (qw0-k*qo)/(qw0+(1-k)*qo)
def qo_func(wct): return qw0*(1-wct**2) / (wct+k*(1-wct**2))

def dp_func(q, y=2.5):
   dp = q**y
   return dp/max(dp)*dp_max

figure()
dp = dp_func(qos)
plot(k*qos, dp, 'b-', label='wat')
plot(qos, dp, 'r-', label='oil')
for i,wct in enumerate(wcts):
   r = i/float(len(wcts))
   qo = qo_func(wct)
   qw = qw_func(qo)
   plot(qo+qw, dp_max, 'o', color=(1-r,0,r))
grid(1)
xlabel('flow rate')
ylabel('dP')
legend(loc='best')
ylim(0, dp_max*1.05)
show()
