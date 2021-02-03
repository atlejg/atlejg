import pandas as pd
import yfinance as yf

#tickr = sys.argv[1]         # typically EQNR.OL
yy = float(sys.argv[1])     # percent yield per year. typically 20 [%]
period1 = int(sys.argv[2])  # [days]. length of streak
period2 = int(sys.argv[3])  # [days]. length of post-streak streak

DEBUG = False

tickrs = unique(tickers[:])

gpd = yy/100 /365.       #  growth per day

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

class Stock(object):
#
    def __init__(self, tickr):
        self.tickr = tickr
        self.invalid = True
        self.success = 0
#
    def get_data(self, period='max'):
        tickr_ = tickr if tickr.endswith('.OL') else tickr+'.OL'     # make sure it's OSL
        s = yf.Ticker(tickr_).history(period=period)
        print('Getting data for', tickr_)
        if len(s) == 0:
            return
        self.invalid = False
        s = s.resample('1d').mean().interpolate()
        #
        self.ps    = s.Close.values
        self.dates = s.index.values
        self.s     = s
#
    def analyze(self, period1, period2, gpd):
        if self.invalid: return
        self.period1  = period1
        self.period2  = period2
        self.gpd      = gpd
        self.streaks1 = []
        self.streaks2 = []
        #
        i0 = 0
        while i0 + period1 < len(self.ps):
            if _is_streak(self.ps[i0:i0+period1], period1, gpd)[0]:
                if DEBUG: print('streak', i0)
                self.streaks1.append(i0)
                i0 = i0 + period1
                r = _is_streak(self.ps[i0:i0+period2], period2, gpd)
                if r[0]: self.streaks2.append([i0, r[1], r[2]])
            i0 = i0 + 1
        #
        self.success = len(self.streaks2)/len(self.streaks1)*100 \
                       if len(self.streaks1) else 0
#
    def plot_streaks2(self):
        if self.invalid or not self.success > 0:
            print(self.tickr, 'No streaks2')
            return
        figure()
        for i0, p, p_lin in self.streaks2:
            plot_date(dates[i0:i0+self.period2], p, 'k-')
            plot_date(dates[i0:i0+self.period2], p_lin, 'k--o')
        title(tickr)
        show()

stocks = []
for tickr in tickrs:
    stock = Stock(tickr)
    stock.get_data()
    stock.analyze(period1, period2, gpd)
    stocks.append(stock)
    print('%-10s - %.1f%%' %(stock.tickr, stock.success))
