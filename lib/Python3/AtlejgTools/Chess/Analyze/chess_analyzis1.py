import chess.engine as ce
import chess
import stockfish
import pickle, sys, os
import chess.pgn as pgn
import numpy as np
import pandas as pd

pgn_file = sys.argv[1]          # typically chess_com_games_2019-09.pgn

if not 'IS_INITIALIZED' in dir():
    IS_INITIALIZED = False

MULTIPV           = 10                   # number of 'best moves' to be returned
TIME_LIM          = 1.0                  # analysis time for chess-engine [s]
STOCKFISH_EXE     = 'c:/MineTing/Tools/stockfish_14_x64_avx2.exe'
PCK_FILE          = '%s_t=%is.pck' % (pgn_file.split('.')[0], TIME_LIM)
READ_HEADERS_ONLY = False
INF               = 999999.


def header_id(h):
    if h['White'] == 'atlejg':
        colour = 'w'
        player = h['Black']
    else:
        colour = 'b'
        player = h['White']
    return '%s:%s:%s' % (h['Date'], colour, player)

def game_id(game):
    return header_id(game.headers)

def get_score(r):
    if r['score'].is_mate():
        if r['score'].white().mate() > 0: return  INF
        else:                             return -INF
    return r['score'].white().cp/100.

def fix(x):
    if len(x) == 4: return
    val = 0 if abs(x[-1]) < np.Inf else x[-1]
    x += [val]*(4-len(x))


def color(ply):
    return 'b' if ply % 2 else 'w'

if not IS_INITIALIZED:
    ENGINE = chess.engine.SimpleEngine.popen_uci(STOCKFISH_EXE)
    IS_INITIALIZED = True

def read_games(pgn_file):
    f = open(pgn_file)
    games   = []
    headers = []
    while True:
        if READ_HEADERS_ONLY:
            h = pgn.read_headers(f)
            if not h: break
            h0 = h.__dict__['_tag_roster']
            headers.append(h0)
        else:
            game = pgn.read_game(f)      # **much** slower to read games, rather than headers
            if not game: break
            games.append(game)
    f.close()
    return games


if not os.path.exists(PCK_FILE):
    res = {}
    games = read_games(pgn_file)
    for game in games[:]:
        gid = game_id(game)
        print('analyzing game', gid)
        score = []
        res[gid] = score
        brd = game.board()
        moves = [x for x in game.mainline_moves()]
        for ply, move in enumerate(moves):
            moveno = f'{(ply//2)+1:d}{color(ply)}'
            print(f'  analyzing move... {moveno} {move}')
            brd.push(move)
            rs = ENGINE.analyse(brd, limit=chess.engine.Limit(time=TIME_LIM), multipv=MULTIPV)
            evaluation = [get_score(r) for r in rs]
            reverse = True if color(ply) == 'w' else False
            evaluation = sorted([get_score(r) for r in rs], reverse=reverse)
            score.append([moveno, move]+evaluation) 
            stop
    for key, r in res.items():
        for x in r: fix(x)
        res[key] = pd.DataFrame(r)
    f = open(PCK_FILE, 'wb')
    pickle.dump(res, f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()
else:
    f = open(PCK_FILE, 'rb')
    res = pickle.load(f)
    f.close()

ms = []
for r in res.values():
    print(len(r))
    if not len(r): continue
    ms.append(np.vstack(r))

#figure()
#for m in ms: plot(m[:,0], 'o')

#print('done')
#print([game_id(game) for game in games])

rr = r.loc[:,2:]

def get_index(col, v, vals):
    if col == 'w' and v >= max(vals): return 0
    if col == 'b' and v <= min(vals): return 0
    if col == 'w' and v <= min(vals): return len(vals)-1
    if col == 'b' and v >= max(vals): return len(vals)-1
    func =  interp1d(vals, range(len(vals)), fill_value='extrapolate')
    return round(float(func(v)))

m = rr.values
actual = []
for ply in range(len(rr)-1):
    v = m[ply+1, 0]
    vals = m[ply,:]
    actual.append(get_index(color(ply), v, vals))

if True:
    figure(figsize=(3, 7))
    imshow(rr.clip(lower=-3, upper=3), origin='lower', cmap='coolwarm', aspect=0.5)
    title(game_id(games[0]))
    xlabel('Stockfish #')
    ylabel('Trekk #')
    plot(actual, range(len(actual)), 'ks--', ms=4)
    show()
