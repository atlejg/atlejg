import pandas as pd
import yfinance as yf

#tickr = sys.argv[1]         # typically EQNR.OL
yy = float(sys.argv[1])     # percent yield per year. typically 20 [%]
period1 = int(sys.argv[2])  # [days]. length of streak
period2 = int(sys.argv[3])  # [days]. length of post-streak streak

DEBUG = False

tickrs = unique(tickers[:])

gpd = yy/100 /365.       #  growth per day

class struct:
    pass



def _is_streak(ps, period, gpd, std_lim=0.1):
    if len(ps) < period: return False, ps, []
    if any(ps<0)       : return False, ps, []
    p0 = ps[0]
    p1 = ps[-1]
    grad = (p1-p0)/p0/(period-1)
    p_lin = p0*(1 + grad*arange(period))
    if grad < gpd: return False, ps, p_lin
    if DEBUG: print(grad, gpd)
    return (std(ps-p_lin)/(p1-p0) < std_lim), ps, p_lin

stocks = []

for tickr in tickrs:
    stock = struct()
    tickr_ = tickr if tickr.endswith('.OL') else tickr+'.OL'     # make sure it's OSL
    s = yf.Ticker(tickr_).history(period='max')
    if len(s) < period1: continue
    s = s.resample('1d').mean().interpolate()
    #
    ps = s.Close.values
    dates = s.index.values
    streaks1 = []
    streaks2 = []
    i0 = 0
    while i0 + period1 < len(ps):
        if _is_streak(ps[i0:i0+period1], period1, gpd)[0]:
            if DEBUG: print('streak', i0)
            streaks1.append(i0)
            i0 = i0 + period1
            r = _is_streak(ps[i0:i0+period2], period2, gpd)
            if r[0]: streaks2.append([i0, r[1], r[2]])
        i0 = i0 + 1
    #
    success = len(streaks2) / len(streaks1) * 100 if len(streaks2) else 0
    print(tickr, 'n-streaks1', len(streaks1))
    print(tickr, 'success = %.1f %%' % success)
    #
    if DEBUG and success > 10:
        figure()
        for i0, p, p_lin in streaks2:
            plot_date(dates[i0:i0+period2], p, 'k-')
            plot_date(dates[i0:i0+period2], p_lin, 'k--o')
        title(tickr)
    stock.tickr    = tickr
    stock.s        = s
    stock.ps       = ps
    stock.dates    = dates
    stock.success  = success
    stock.streaks1 = streaks1
    stock.streaks2 = streaks2
    stock.period1  = period1
    stock.period2  = period2
    stock.gpd      = gpd
    stock.yy       = yy
    stocks.append(stock)


show()
