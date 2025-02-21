import chess.engine as ce
import chess
import stockfish
import pickle, sys, os
import chess.pgn as pgn
import numpy as np
import pandas as pd

pgn_file  = sys.argv[1]          # typically chess_com_games_2019-09.pgn
do_pickle = int(sys.argv[2])

if not 'IS_INITIALIZED' in dir():
    IS_INITIALIZED = False

N_BEST_MOVES      = 10                   # number of 'best moves' to be returned
TIME_LIM          = 1.0                  # analysis time for chess-engine [s]
STOCKFISH_EXE     = 'C:\Appl\stockfish-windows-x86-64-avx2.exe'
PCK_FILE          = '%s_t=%is.pck' % (pgn_file.split('.')[0], TIME_LIM)
INF               = 999999.

#pm = {'Minimum Thinking Time':2000}
sf = stockfish.Stockfish(path=STOCKFISH_EXE, depth=15, parameters={})

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

def fix(x):
    if len(x) == 4: return
    val = 0 if abs(x[-1]) < np.Inf else x[-1]
    x += [val]*(4-len(x))


def color(ply):
    return 'b' if ply % 2 else 'w'

def read_games(pgn_file):
    f = open(pgn_file)
    games   = []
    headers = []
    while True:
        game = pgn.read_game(f)
        if not game: break
        games.append(game)
    f.close()
    return games

def evaluate(top, move, col, n_vals=N_BEST_MOVES):
    #
    act = -1
    vals = []
    for i, m in enumerate(top):
        if m['Move'] == move:
            act = i
        if m['Mate']:
            cp = INF if col=='w' else -INF
        else:
            cp = m['Centipawn'] / 100
        vals.append(cp)
    #
    if len(vals) < n_vals:
        vals = vals + [NaN]*(n_vals-len(vals))
    #
    if act < 0:
        act = i+1
    #
    return vals, act


if not os.path.exists(PCK_FILE):
    res = {}
    games = read_games(pgn_file)
    for game in games:
        gid = game_id(game)
        print('analyzing game', gid)
        res[gid] = {'eval':[], 'actual':[]}
        moves = [str(x) for x in game.mainline_moves()]
        for ply, move in enumerate(moves):
            moveno = f'{(ply//2)+1:d}{color(ply)}'
            print(f'  analyzing move... {moveno} {move}')
            top = sf.get_top_moves(N_BEST_MOVES)
            vals, act = evaluate(top, move, color(ply))
            res[gid]['eval'].append([moveno, move]+vals) 
            res[gid]['actual'].append(act)
            sf.make_moves_from_current_position([move])
    for key in res.keys():
        res[key]['eval'] = pd.DataFrame(res[key]['eval'])
    #
    if do_pickle:
        f = open(PCK_FILE, 'wb')
        pickle.dump(res, f)
        f.close()
else:
    f = open(PCK_FILE, 'rb')
    res = pickle.load(f)
    f.close()

r = list(res.values())[0]
rr = r['eval'].loc[:,2:]
actual = r['actual']
act_w = actual[::2]
act_b = actual[1::2]
n_moves_w = list(range(0, len(actual), 2))
n_moves_b = list(range(1, len(actual), 2))


if 1:
    n = len(rr)
    clf()
    #figure(figsize=(3, 7))
    aspect = min(15 * N_BEST_MOVES / n, 0.4)
    ms = min(3, 3*(100/n))
    subplot(121)
    imshow(rr.clip(lower=-3, upper=3), origin='lower', cmap='coolwarm', aspect=aspect)
    title(pgn_file.split('.')[0])
    xlabel('Stockfish #')
    ylabel('Trekk #')
    xticks(range(N_BEST_MOVES), ['' for x in range(N_BEST_MOVES)])
    yticks(range(0,n,20), [f'{x/2:.0f}' for x in range(0,n,20)])
    #colorbar(location='bottom', fraction=0.1, shrink=0.3)
    #subplot(122, xscale=0.4)
    subplot(122)
    plot(act_w, n_moves_w, 'ro', ms=ms)
    plot(act_b, n_moves_b, 'bs', ms=ms)
    ylim(0, len(actual))
    xticks(range(N_BEST_MOVES), ['' for x in range(N_BEST_MOVES)])
    yticks([])
    show(block=False)
