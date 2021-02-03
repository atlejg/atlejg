import pandas as pd
import yfinance as yf

tickr = sys.argv[1]         # typically EQNR.OL
yy = float(sys.argv[2])     # percent yield per year. typically 20 [%]
period1 = int(sys.argv[3])  # [days]. length of streak
period2 = int(sys.argv[4])  # [days]. length of post-streak streak

oneday = datetime.timedelta(1)
period1 = datetime.timedelta(period1)  # days
period2 = datetime.timedelta(period2)  # days

growth = yy/100 * period1.days/365

#s = yf.Ticker(tickr).history(period1='max')
#s = s.resample('1d').mean().interpolate()

def _is_streak(c, period, growth, std_lim=0.1):
    c0 = c[0]
    c1 = c[-1]
    grad = (c1-c0)/(period.days-1)
    if grad/c0 < growth: return False
    c_lin = c0 + grad*arange(period.days)
    if std(c-c_lin)/(c1-c0) < std_lim: return True
    return False

def _evaluate(c, period, growth):
    c_lin = c[0] * (1 + growth/period.days*arange(period.days))
    figure()
    plot(c)
    plot(c_lin)
    return all(greater_equal(c, c_lin)), c, c_lin

streaks = []
res = []
t0 = s.index[0]
while True:
    if t0 + period1 > s.index[-1]: break
    c = s.loc[t0:t0+period1-oneday].Close.values
    if _is_streak(c, period1, growth):
        print('streak', t0)
        streaks.append(t0)
        t0 = t0 + period1
        res.append(_evaluate(s.loc[t0:t0+period2-oneday].Close.values, period2, growth))
    t0 = t0 + oneday


#for i in range(len(res)):

show()
