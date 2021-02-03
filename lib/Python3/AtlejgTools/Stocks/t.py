import pandas
import pandas as pd
import yfinance as yf
yf
r = yf.Ticker('EQNR.OL')
r
r.info
r = yf.Ticker('EQNR.OL')
r = yf.Ticker?
r
r.info?
r.history(proxy='http://www-proxy.statoil.no:80')
r.history()
h = r.history()
len(h)
h
len(yf.Ticker('EQNR.OL').history())
len(yf.Ticker('MSFT').history())
len(yf.Ticker('MSFT').history(period='max'))
len(yf.Ticker('EQNR.OL').history(period='max'))
h = yf.Ticker('EQNR.OL').history(period='max')
h
h.plot()
h.plot(['Close'])
h.plot('Date',['Close'])
h
h.plot?
h.plot(x='Date',y=['Close'])
h
h.Date
h.keys
h.keys()
h.plot(y=['Close'])
c = h.Close
c = h.Close.values
c
max(c)
figure()
plot(c)
c = h.Low.values
l = h.Low.values
c = h.Close.values
plot(l)
figure()
plot(l)
h
h.max
h.max()
'atle'.split('r')
'atle'.split('r')[1]
pwd
ls
run -i tickers.py stocklist.osl
run -i tickers.py stocklist.osl
tickers
run -i tickers.py stocklist.osl
history
len(yf.Ticker('EQNR.OL').history(period='max'))
[len(yf.Ticker(ticker).history(period='max')) for ticker in tickers[:10]]
hs = [yf.Ticker(ticker).history(period='max') for ticker in tickers]
len(hs)
hs
h
h.columns
h.index
run -i tickers.py stocklist.osl
len(hs)
len(tickers)
run -i tickers.py stocklist.osl
h
c
c
len(c)
c = h.Close.values
type(c)
c = h.Close
type(c)
c.is_monotonic_increasing?
c.is_monotonic?
c.is_monotonic??
ygr = 20 # percent yield per year
yy = 20 # percent yield per year
period = 30  # days
growth = yy/100 * period/365
growth
h
h.index
t0 = h.index[-100]
t0
t0 = h.index[-300]
t0
h[t0]
h
h.at[t0]
h
h.iat[t0]
t0
h.at?
h.loc(t0)
h.loc?
h.loc[t0]
h.loc[t0+period]
t0+period
datetime.timedelta(period)
datetime.timedelta(period) + t0
h.loc[datetime.timedelta(period) + t0]
h.loc[datetime.timedelta(period) + t0].Close
h.loc[datetime.timedelta(0) + t0].Close
h.interpolate?
h.resample?
h
h.resample('1d').mean()
h.resample('1d').mean?
h.resample('1d').mean
h.resample('1d').ffill()
h.resample('1d')
h.resample('1d').mean()
hh = h.resample('1d').mean()
hh.interpolate()
h.resample('1d').mean().interpolate()
h = h.resample('1d').mean().interpolate()
len(h)
history -f t.py
c = h.loc[t0]
c = h.loc[t0].Close
c0 = h.loc[t0].Close
c0 = h.loc[t0].Close
period = datetime.timedelta(30)  # days
c1 = h.loc[t0+period].Close
(c1-c0)/c0
growth
c1
c0
h
h[1:20]
h.loc[t0:t0+period]
h.loc[t0:t0+period].Close
h.loc[t0:t0+period].Close.values
plot(h.loc[t0:t0+period].Close.values)
figure()
figure()
plot(h.loc[t0:t0+period].Close.values)
yy = c0 + (c1-c0)/period*arange(period)
period.days
yy = c0 + (c1-c0)/period*arange(period.days)
yy = c0 + (c1-c0)/period.days*arange(period.days)
yy
plot(yy)
h.loc[t0:t0+period].Close.values
close('all')
plot(h.loc[t0:t0+period].Close.values)
plot(yy)
len(h.loc[t0:t0+period])
len(yy)
yy = c0 + (c1-c0)/period.days*arange(period.days)
period.days
len(h.loc[t0:t0+period])
len(h.loc[t0:t0+period-1])
oneday = datetime.timedelta(1)
len(h.loc[t0:t0+period-oneday])
clf()
plot(h.loc[t0:t0+period-oneday].Close.values)
plot(yy)
yy
yy = c0 + (c1-c0)/period.days*arange(period.days)
yy[0]
c0
c1
yy[-1]
c1 = h.loc[t0+period].Close
h.loc[t0+period].Close
c1
h.loc[t0:t0+period-oneday].Close.values
h.loc[t0+period].Close
h.loc[t0:t0+period].Close.values
c1 = h.loc[t0+period-oneday].Close
c1
yy = c0 + (c1-c0)/period.days*arange(period.days)
yy[-1]
yy
c0,c1
yy
len(yy)
yy = c0 + (c1-c0)/(period.days-1)*arange(period.days)
yy
clf()
plot(h.loc[t0:t0+period-oneday].Close.values)
plot(yy)
dc = c1 - c0
std?
var?
std?
history
c = h.loc[t0:t0+period-oneday].Close.values
cc = c0 + (c1-c0)/(period.days-1)*arange(period.days)
clf()
plot(c)
plot(cc)
std(c -cc)
c1 - c0
history -f t.py
