import re

l1 = open('covid19_testing_20200404.html').readlines()

l2 = []

rec = []
for l in l1:
    m = re.search('coronavirus pandemic in (.*)</a>', l)
    if m:
        if rec: l2.append(rec)
        rec = [m.groups()[0]]
        continue
    m = re.search('/span><span>(\\d.*)$', l)
    if m:
        rec.append(m.groups()[0])
        continue
l2.append(rec)


data = {}
for rec in l2:
    s = UT.Struct()
    s.nm     = rec[0].replace('the ','')
    if s.nm in data: continue
    if len(rec) < 6:
        print('cannot handle', s.nm)
        continue
    s.tests  = int(rec[1].replace(',','').replace(' ','').replace('*',''))
    s.pos    = int(rec[2].replace(',','').replace(' ','').replace('*',''))
    s.test_r = float(rec[4].replace(',','').replace(' ','').replace('*',''))
    s.pos_r  = float(rec[5].replace(',','').replace(' ','').replace('*',''))
    data[s.nm] = s
