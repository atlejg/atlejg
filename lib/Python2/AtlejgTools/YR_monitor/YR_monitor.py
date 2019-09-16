#!/usr/bin/env python
# -*- coding: utf-8 -*-

import urllib
import sys, os, glob
import datetime, time
import AtlejgTools.Utils as UT
from pylab import *

# constants
NPERIODS_ADAY   = 4
PERIODE_TIME    = 24 / NPERIODS_ADAY
N_FORECAST_DAYS = 9   # number of days for the forecast
YR_URL_PREFIX   = 'http://www.yr.no/sted'

# globals
_rootdir = ''
_d       = ''

def _usage():
   return '''
   usage: %s -[adf] inifile
   Use -a for anlayzing and -d for downloading once, -f for downloading forever
   where inifile is like:
[params]
location  = Norge/Telemark/Kragerø/Jomfruland
rootdir   = .
   ''' % sys.argv[0]


def _read_xml_forecast(fname):
   '''
   makes sure only the part inside the tabular-tag is included
   '''
   fp = open(fname)
   lines = []
   skip = True
   for line in fp:
      if 'tabular' in line: skip = False
      if skip: continue
      lines.append(line)
      if '/tabular' in line: break
   fp.close()
   return lines

def _get_value(line, split_by):
   return float(line.split(' %s="' % split_by)[1].split('"')[0])

def _periode(hour):
   return int(hour) / PERIODE_TIME

def download(location):
   now = datetime.datetime.now()
   path = _rootdir + '/' + location
   if not os.path.exists(path): os.makedirs(path)
   fname = '%s/%i-%02i-%02i_%i.xml' % (path, now.year, now.month, now.day, _periode(now.hour))
   url = "%s/%s/varsel.xml" % (YR_URL_PREFIX, location)
   urllib.urlretrieve(url, fname)
   print 'saving file', fname
   return fname

def download_forever(location, interval=None):
   '''
   keeps downloading until killed.
   interval is in hours. default is PERIODE_TIME
   '''
   if interval==None: interval = PERIODE_TIME
   while True:
      download(location)
      time.sleep(interval*3600)

def _time(y, m, d, h):
   return date2num(datetime.datetime(y, m, d, h))

def _time_from_tstamp(tstamp):
   date, period = tstamp.split('_')
   y, m, d = date.split('-')
   return _time(int(y), int(m), int(d), int(period)*PERIODE_TIME)

def _time_from_fname(fname):
   tstamp = os.path.split(fname)[1].split('.')[0]
   return _time_from_tstamp(tstamp)

def _nperiods(fname1, fname2):
   dt = _time_from_fname(fname1) - _time_from_fname(fname2)
   return int(dt / (1./NPERIODS_ADAY)) + 1

def _nan(nx, ny):
   d    = zeros((nx, ny))
   d[:] = NaN
   return d

def _nperiods_forecast(fname):
   tstamp = os.path.split(fname)[1].split('.')[0]
   period = int(tstamp.split('_')[1])
   np = N_FORECAST_DAYS * NPERIODS_ADAY
   if period % 2 == 1: np = np - 1
   return np

def analyze(location, pattern, ndays=999999, plotit=False):
   '''
   plots contour-maps showing how the forecast is changing
   '''
   fnames = glob.glob(location+'/%s*.xml' % pattern)
   nfiles = NPERIODS_ADAY * ndays
   if nfiles < len(fnames): fnames = fnames[-nfiles:]
   fnames.sort()
   # add 2 because it seems to vary from file to file :-(
   nperiods = _nperiods(fnames[-1], fnames[0]) + _nperiods_forecast(fnames[-1]) + 2
   print nperiods
   temp  = _nan(nperiods, N_FORECAST_DAYS*NPERIODS_ADAY + 2)
   prec  = _nan(nperiods, N_FORECAST_DAYS*NPERIODS_ADAY + 2)
   windD = _nan(nperiods, N_FORECAST_DAYS*NPERIODS_ADAY + 2)
   windS = _nan(nperiods, N_FORECAST_DAYS*NPERIODS_ADAY + 2)
   pres  = _nan(nperiods, N_FORECAST_DAYS*NPERIODS_ADAY + 2)
   for fname in fnames:
      print fname
      offset = _nperiods(fname, fnames[0])
      lines = _read_xml_forecast(fname)
      np = -1
      for line in lines:
         if 'time from'     in line: np += 1
         if 'temperature'   in line: temp[offset + np, np]  = _get_value(line, 'value')
         if 'precipitation' in line: prec[offset + np, np]  = _get_value(line, 'value')
         if 'windDirection' in line: windD[offset + np, np] = _get_value(line, 'deg')
         if 'windSpeed'     in line: windS[offset + np, np] = _get_value(line, 'mps')
         if 'pressure'      in line: pres[offset + np, np] = _get_value(line, 'value')
   #
   if plotit:
      figure()
      contourf(temp.T)
      title('temperature [Celsius]')
      xlabel('period')
      ylabel('delta period')
      colorbar()
      #
      figure()
      contourf(prec.T)
      title('precipitation [mm/period]')
      xlabel('period')
      ylabel('delta period')
      colorbar()
      #
      figure()
      contourf(windS.T)
      title('wind speed [m/s]')
      xlabel('period')
      ylabel('delta period')
      colorbar()
      #
      figure()
      contourf(windD.T)
      title('wind direction [deg]')
      xlabel('period')
      ylabel('delta period')
      colorbar()
      #
      show()
   s = UT.Struct()
   s.temp = temp
   s.prec  = pres
   s.windD = windD
   s.windS = windS
   s.pres  = pres
   return s

def analyze2(location, pattern='*', plotit=True):
   '''
   just get the first value (more or less current time) from each file.
   '''
   fnames = glob.glob(location+'/%s*.xml' % pattern)
   nfiles = len(fnames)
   fnames.sort()
   temp  = []
   prec  = []
   windD = []
   windS = []
   pres  = []
   t     = []
   for fname in fnames:
      print fname
      t.append(_time_from_fname(fname))
      lines = _read_xml_forecast(fname)
      for line in lines:
         if 'temperature'   in line: temp.append(_get_value(line, 'value'))
         if 'precipitation' in line: prec.append(_get_value(line, 'value'))
         if 'windDirection' in line: windD.append(_get_value(line, 'deg'))
         if 'windSpeed'     in line: windS.append(_get_value(line, 'mps'))
         if 'pressure'      in line: pres.append(_get_value(line, 'value'))
         if 'pressure'      in line: break
   #
   if plotit:
      if temp:
         fig = figure()
         plot_date(t, temp)
         fig.autofmt_xdate()
         grid(1)
         ylabel('temperature [Celsius]')
         xlabel('time')
      #
      if prec:
         fig = figure()
         plot_date(t, prec)
         fig.autofmt_xdate()
         grid(1)
         ylabel('precipitation [mm/period]')
         xlabel('time')
      #
      if windS:
         fig = figure()
         plot_date(t, windS)
         fig.autofmt_xdate()
         grid(1)
         ylabel('wind speed [m/s]')
         xlabel('time')
      #
      if windD:
         fig = figure()
         plot_date(t, windD)
         fig.autofmt_xdate()
         grid(1)
         ylabel('wind direction [deg]')
         xlabel('time')
      #
      if pres:
         fig = figure()
         plot_date(t, pres)
         fig.autofmt_xdate()
         grid(1)
         ylabel('Pressure [hPa]')
         xlabel('time')
      show()
   s = UT.Struct()
   s.temp = temp
   s.prec  = pres
   s.windD = windD
   s.windS = windS
   s.pres  = pres
   return s
   s = UT.Struct()
   s.temp = array(temp)
   s.prec  = array(pres)
   s.windD = array(windD)
   s.windS = array(windS)
   s.pres  = array(pres)
   return s

def _read_stats_html(location): pass

def _evaluate(location):
   '''
   not in use...
   '''
   global _d
   data = _read_stats_html(location)
   url = '%s/%s/detaljert_statistikk.html' % (YR_URL_PREFIX, location)
   website = urllib.urlopen(url)
   lines = website.readlines()
   website.close()
   _d = lines
   started = True
   count = 0
   data = []
   for line in lines:
      if '/thead' in line:
         started = True
         continue
      if started and '/tbody' in line: break
      if started and '<tr>' in line: d.append([])
      if started and '<td>' in line: pass

if __name__ == '__main__':

   if len(sys.argv) < 3:
      print(_usage())
      sys.exit(1)
   ini = UT.read_inifile(sys.argv[2])
   _rootdir = ini.rootdir
   if   sys.argv[1] == '-a':
      if len(sys.argv) == 4: pattern = sys.argv[3]    # like 2012-2
      else                 : pattern = '*'
      s = analyze(ini.location, pattern, int(ini.ndays))
   elif sys.argv[1] == '-d': download(ini.location)
   elif sys.argv[1] == '-f': download_forever(ini.location)
   #elif sys.argv[1] == '-e': _evaluate(ini.location)
