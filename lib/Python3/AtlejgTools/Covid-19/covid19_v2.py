#!/usr/bin/env python
# coding: utf-8

'''
 pip install covid-data-api
 https://pypi.org/project/covid-data-api/
'''
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

from covid.api import CovId19Data

testing_url = "https://ourworldindata.org/d08107e0-17d5-4c5c-bc09-5e8af43ee10b"

class CovidCountry(object):
#
    def __init__(self, nm, days, dates, sick, dead, tst=None, hosp=None, icu=None, rec=None):
        self.nm     = nm
        self.cnm    = nm.capitalize()   # capitalized name
        self.days   = days
        self.dates  = dates
        self.tst    = tst
        self.sick   = sick
        self.dead   = dead
        self.hosp   = hosp
        self.icu    = icu
        self.rec    = rec
        self.nd     = len(days)
        #
        if self.nd:
            # relative values (in percent)
            self.sick_r = self.sick/max(self.sick)*100
            self.dead_r = self.dead/max(self.dead)*100
            #
            self.platou = self._platou()                # boolean
            self.delay  = self._delay()                 # lag between sick and dead
            #
            self.dr = self.dead[-1] / self.sick[-1]     # death-rate
        else:
            self.sick_r = r_[]
            self.dead_r = r_[]
            self.platou = False                         # boolean
            self.delay  = 0                             # lag between sick and dead
#
    def _platou(self, ndays=20, lvl=95):
        return len((self.sick_r > lvl).nonzero()[0]) > ndays
#
    def _delay(self):
        y1 = self.sick_r
        y2 = self.dead_r
        n = len(y1)
        v = []
        for i in range(1,int(n/2)):
            v.append(norm(y1[:-i]-y2[i:])/(n-i)) 
        return argmin(v) + 1
#
    def plot_delayed_deads(self):
        figure()
        plot(self.days, self.sick_r)
        plot(self.days-self.delay, self.dead_r)
        legend(['sick', 'dead'])
        title('%s - delay = %i days' % (self.cnm, self.delay))
        grid(1)
        ylabel('%')
        xlabel('Days since start of Covid-19 in this country')
        show()

# confusing. using both country names and 'labels' for lookups.
#_countries = list(cv.get_all_records_by_country().keys())
#_map    = {}
#for c,l in zip(_countries, _labels): _map[c] = l


CONFIRM_LIMIT = 80    # indicates day-0 of covid
CONFIRM_YMIN  = 80
FATAL_LIMIT   = 2    # indicates day-0 of covid
FATAL_YMIN    = 2.5

force_read = int(sys.argv[1])
varnm        = sys.argv[2]

# Read country data sets
cov19 = CovId19Data(force=force_read)

_countries = [x['label'].lower() for x in list(cov19.get_all_records_by_country().values())]
# note: countries *must* be in _countries
countries = ['norway', 'sweden', 'denmark', 'korea, south', 'italy', 'spain', 'us']
countries = ['norway', 'sweden', 'denmark', 'italy', 'spain', 'us', 'china']
countries = _countries


if varnm == 'confirmed':
    y_lim = CONFIRM_LIMIT
    y_min = CONFIRM_YMIN
elif varnm == 'deaths':
    y_lim = FATAL_LIMIT
    y_min = FATAL_YMIN
else:
    raise Exception('No such variable: %s' % varnm)

def get_history(country):
    data = list(cov19.get_history_by_country(country).values())[0]['history']
    df = pd.DataFrame.from_dict(data, orient='index')
    df.index = pd.to_datetime(df.index)
    df.replace(r'na', np.nan,inplace=True)
    df.change_confirmed = pd.to_numeric(df.change_confirmed, downcast='float')
    df.change_deaths = pd.to_numeric(df.change_deaths, downcast='float')
    df['time_to_double'] = np.log(2)/np.log(1 + df.change_confirmed)
    return df

dfs = {}
cs  = []
ymax = -1
dmax = -1
for country in countries:
    try:
        df = get_history(country)
        df = df[df[varnm] >= y_lim]
        if len(df) > 0:
            df['days'] = [d.days for d in (df.index - df.index[0])]
            ymax = max(ymax, nanmax(df[varnm]))
            dmax = max(dmax, nanmax(df['days']))
        else:
            df['days'] = []
        dfs[country] = df
        c = CovidCountry(country, df['days'].values, df.index.values, df['confirmed'].values, df['deaths'].values)
        cs.append(c)
    except:
        print('could not handle', country)

def doubling(days, dbl, y_lim):
    return y_lim*2**(days/dbl)

figure()

# plot real data
for country, df in dfs.items():
    #print(country)
    semilogy(df['days'], df[varnm], basey=2, label=country)

# plot doubling-lines
days = r_[0, dmax]
for dbl in [2, 7]:
    semilogy(days, doubling(days, dbl, y_lim),  '--', color='k', lw=2,  basey=2, label='doubling every %i days'%dbl)

#nmax = int(log(ymax/y_min)/log(2)) + 1
#yt = y_min*array([2**n for n in range(nmax)])
yt_base = r_[1, 2.5, 5]
expon = 0
yt = []
i  = 0
while True:
    y = yt_base[i]*10**expon
    if y >= y_min: yt.append(y)
    if y > ymax: break
    i += 1
    if i == len(yt_base):
        i = 0
        expon += 1
yticks(yt, ['%d'%y for y in yt])

ylim(0.9*y_lim, 1.05*ymax)
legend(loc='best')
ylabel(varnm)

show()
