fnm = sys.argv[1]

tickers = []
for line in open(fnm).readlines():
    if not 'ose|' in line: continue
    ticker = line.split('ose|')[1][:-3]
    print(line, ticker)
    tickers.append(ticker)
