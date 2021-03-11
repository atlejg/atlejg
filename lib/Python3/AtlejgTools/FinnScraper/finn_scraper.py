#! /prog/sdpsoft/python3.6.4/bin/python3

'''
Created on Jan 18, 2021
@author: Atle J. Gyllensten

This script will 'scrape' data from finn.
It was first and foremost made for extracting all real-estate data.
We need to fetch data 'kommune' by 'kommune' (called locations) since
finn will only show up to 50 'pages' for each search. So to make sure
we will have less than 50 pages for each search, we search 'kommune' by 'kommune'

It expects a YAML-input file like this:

    url_templ:
        https://www.finn.no/realestate/homes/search.html?%s&sort=PUBLISHED_DESC&is_new_property=false

    location_file:             # if location_file is empty, it will use all locations

    max_pages:
        !!int   51

    sleep_tm:                  # number of secs to sleep between each location (dont wanna stress web-server)
        !!float 2

    prefix:                    # prefix of the csv-files. if empty, will use basename of *this* file

    save_to_file:
        !!bool true

NOTES:
    n1: I'm not sure it works if is_new_property=true (in the url_templ)

'''

from bs4 import BeautifulSoup
import urllib.request
import sys  
import ssl
import time
import pandas as pd
import datetime
import yaml
import pathlib

MAX_PAGES = 51

MAPPER = {
 'location'                 :(1,0, str),
 'price_suggestion'         :(2,1, float),
 'price_total'              :(2,1, float),
 'price_shared_cost'        :(2,1, float),
 'size_from'                :(1,1, float),
 'size_to'                  :(1,1, float),
 'number_of_bedrooms'       :(1,1, int),
 'organisation_name'        :(1,0, str),
 'local_area_name'          :(1,0, str),
 'owner_type_description'   :(1,0, str),
 'property_type_description':(1,0, str),
 'heading'                  :(1,0, str),
}

# fix for CERTIFICATE_VERIFY_FAILED,
# see https://gankrin.org/how-to-fix-ssl-certificate_verify_failed-error-error-in-python/
ssl._create_default_https_context = ssl._create_unverified_context   


def map_data(rec, loc):
    data = [loc]
    for keyno, key in enumerate(MAPPER.keys()):
        for i, item in enumerate(rec):
            if item == key:
                j,k, typ = MAPPER[key]
                val = rec[i+j][k:]
                if typ in (float ,int):
                    try:    val = typ(val.replace(',',''))
                    except: val = -1
                data.append(val)
                break
        #
        # set default values for missing data
        if not len(data) == keyno+2: 
            data.append(None)
    return data

def get_html_txt(url):
    #req = urllib.request.Request(url, headers={'User-Agent':'Mozilla/5.0'})
    req = urllib.request.Request(url)
    ir = urllib.request.urlopen(req)
    ih = ir.read()
    bs = BeautifulSoup(ih,"html.parser")
    return bs.find_all('script')[11].getText()

def get_all_locations(url):
    txt = get_html_txt(url)
    lines = txt.split('display_name')[1:-1]
    locs = []
    for line in lines:
        if 'location' in line and 'filter_items' in line and 'value' in line:
            line = line.replace('"','')
            recs = line.split(':')
            if len(recs) < 3: continue
            loc =  recs[3].split(',')[0]
            if loc.startswith('0'): continue
            name = recs[1].split(',')[0]
            locs.append((loc, name))
    return locs

def get_locations(fnm):
    '''
    fnm: file of locations - just plain list of 'kommune-id' & 'kommune-navn'
    '''
    return [x.strip().split() for x in open(fnm).readlines()]

def get_finn_data(url_templ, locs, max_pages=MAX_PAGES, sleep_tm=0., verbose=True):
    '''
    sleep_tm: number of secs to sleep between each location (dont wanna stress web-server)
    '''
    nap_tm = sleep_tm / 10.
    res = []
    for loc in locs:
        page = 1
        while True:
            if verbose: print(loc[1], page)
            txt = get_html_txt(url_templ%('location=%s&page=%i'%(loc[0],page)))
            # first get finn-codes (not the same place as the rest of the data)
            recs = txt.split('finnkode=')[1:-1]
            finncodes = [rec.split('"')[0] for rec in recs]
            #
            # now, get the rest of the data
            recs = [x.split('ad_link')[0] for x in txt.split('SEARCH_ID_REALESTATE')[1:-1]]
            for finncode, rec in zip(finncodes,recs):
                data = [x for x in rec.split('"') if not x in (":", ",", "'", "},", ":{")]
                rec = map_data(data, loc[1])
                res.append([finncode]+rec)
            if not len(recs): break
            page += 1
            if page > max_pages: break
            if nap_tm > 0: time.sleep(nap_tm)
        if sleep_tm > 0: time.sleep(sleep_tm)
    #
    return pd.DataFrame(res, columns=['finncode','area']+list(MAPPER.keys()))

def save_cvs(df, prefix, verbose=True):
    now = datetime.datetime.now()
    fnm = '%s_%i_%02i_%02i.csv' % (prefix, now.year, now.month, now.day)
    df.to_csv(fnm, sep=';')
    if verbose: print(fnm, 'was created')




if __name__ == '__main__':

    inp_file = sys.argv[1]
    inp = yaml.load(open(inp_file))

    if inp['location_file']:
        locs = get_locations(inp['location_file'])
    else: 
        locs = get_all_locations(inp['url_templ']%'')

    max_pages = inp['max_pages'] if inp['max_pages'] else MAX_PAGES

    finn = get_finn_data(inp['url_templ'], locs, max_pages=max_pages, sleep_tm=inp['sleep_tm'])

    if inp['save_to_file']:
        prefix = inp['prefix'] if inp['prefix'] else pathlib.Path(inp_file).stem
        save_cvs(finn, prefix)
