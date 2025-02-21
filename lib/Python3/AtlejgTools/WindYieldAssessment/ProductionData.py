#!/usr/bin/env python

'''
agy@equinor.com
2022-02-25

reads production data from Bazefield.
based on daily_wing_gen.ipynb by Claire.
simplified by getting all data in one go.
this is too simple if we need PBA (production availability) per time-step.

to browse the production data 'manually', and where to find the github library, 
see links from Claire in Teams-meeting-chat 25/2-2022

input is a yaml like this:

    farm_nm: DOW
    timezone: Europe/London
    date1: 2017-01-01
    date2: 2021-01-01             # special entry: today
    resolution: 1h                # 10m, 30m, 1h, 6h, 1d, 7d ...

notes
 * i do not use ProducedMWh.AD.SUM (just ProducedMWh) since averaging out
   the data, gives constant value for each 24h.
 * tagnames for various assets/turbines can be found here
   https://statoilsrm.sharepoint.com/sites/DPA/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2Fsites%2FDPA%2FShared%20Documents%2FData%2FBazefield%2FTaglists%2FTaglist%20%2D%20turbine
'''

import os, sys, time
import pandas as pd
import datetime
import  AtlejgTools.Utils as UT

libdir = r'D:\OneDrive - Equinor\ACTIVE\Wind\BazeStats\ReportingTools'
assert os.path.exists(libdir), f'necessary bazestat-library {libdir} does not exist!'
if not libdir in sys.path: sys.path.append(libdir)

from bazestats import BazeStats

ALLOC_TYPE = 7                                     # 7 = Equinor Availability (not in use here..)
KEYS       = ['export', 'ws']                      # use these column names
UNITS      = dict(zip(KEYS, ['MWh', 'm/s']))

# authorization. Claire uses an azure key-vault.
# for now, i just use system environment (ENV)
API_KEY = os.environ.get('BF_API_KEY')             # got this from Claire. should have my own!

def bazefield_tags(keys, farm_nm, allowed_keys=KEYS):
    farmSiteTag = "SITE-PI" if farm_nm == 'SHS' else "PI"
    #
    tags = []
    for key in keys:
        assert key in allowed_keys
        if key == 'export':
            tags.append(f'{farm_nm}-{farmSiteTag}-ProducedMWh')
        elif key == 'ws':
            tags.append(f'{farm_nm}-Site-WindSpeed-WTGs')
    #
    return tags

def get_bazefield_data(yaml_file=None, bz=None, resolution=None, date1=None, date2=None, timezone=None,
                api_key=API_KEY, keys=KEYS, alloc_type=ALLOC_TYPE):
    '''
    get production data from Bazefield.
    - inputs
      * yaml_file   input file specifying the data to be retrieved
                    set to None if you use bz
                    yaml_file *must* look like:
                        farm_nm: DOW
                        timezone: Europe/London
                        date1: 2019-01-01
                        date2: today
      * bz:         a BazeStats object.
                    set to None if you use yaml_file
                    it's not really that much to gain by re-using the
                    bz-object, so might as well just stick to using yaml_file option
      * resolution: 10m, 30m, 1h, 6h, 1d, 7d ...
      * date1:      from-date [YYYY-MM-DD].
                    if None, use from-date from bz-object or the yaml
      * date2:      to-date [YYYY-MM-DD] or special key 'today'
                    if None, use to-date from bz-object or the yaml
      * timezone:   typically 'Europe/London' (TODO: check!)
      * api_key:    *must* be provided if yaml_file is provided
      * keys:       list of keys. each key must be an allowed key
                    (see bazefield_tags)
    - returns
      * data, as DataFrame
      * BazeStats object (if yaml_file was provided)
    - notes
      * either bz OR yaml_file *must* be provided
      * will set variable _timezone for bz-object since this useful when re-using the bz
    '''
    if yaml_file:
        inp = UT.get_yaml(yaml_file)
        if not timezone: timezone = inp.timezone
        if not resolution: resolution = inp.resolution
        if not date1: date1 = inp.date1
        if not date2: date2 = inp.date2
        bz = BazeStats(api_key, inp.farm_nm)
        bz._timezone = timezone
        #
    else:
        assert not resolution is None, 'resolution *must* be provided'
        if not timezone: timezone = bz._timezone
        if not date1: date1 = bz.startTime
        if not date2: date2 = bz.endTime
    #
    if type(date2) is str and date2.lower() == 'today':
        date2 = datetime.datetime.now()
    date1 = pd.to_datetime(f'{date1} 00:00:00').tz_localize(timezone)
    date2 = pd.to_datetime(f'{date2} 04:00:00').tz_localize(timezone)
    #
    # update BazeStats-object with new dates
    bz(date1, date2, alloc_type)
    #
    # get data
    tags = bazefield_tags(keys, bz.farmName)
    aggrs = ['TIMEAVERAGE']*len(tags)
    data = bz.MultiTagTimeSeries(tags, aggrs, resolution)
    data.columns = keys
    #
    if yaml_file:
        return data, bz
    else:
        return data

if __name__ == '__main__':

    tic = time.time()

    yaml_file = sys.argv[1]                       # typically dow.yaml
    if len(sys.argv) == 3:
        resolution = sys.argv[2]              # typically 10m, 30m, 1h, 6h, 1d, 7d ...

    data1, bz = get_bazefield_data(yaml_file, resolution=resolution)

    #wtgs = bz.farm                 # name, location etc.
    #pe = bz.GetAllocationStats()   # production efficency. note: SLOW!

    toc = time.time()
    print(f'runtime: {(toc-tic):.1f} secs. n-entries: {len(data1)}. resolution: {resolution}')
