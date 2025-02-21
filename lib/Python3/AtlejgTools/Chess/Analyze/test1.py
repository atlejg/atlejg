import chess.engine as ce
import chess
import stockfish
history
pwd
cd OneDrive/PASSIVE/Chess/PGNs/Mine/Chess.com/
ls
b = chess.Board()
import chess.pgn as pgn
pgn
f = open('chess_com_games_2019-09-01.pgn')
g1 = pgn.read_game(f)
g1
g1.headers
b = g1.board()
import stockfish

for move in g1.mainline_moves():
    b.push(move)
board
b
for move in g1.mainline_moves():
    b.push(move)
    print(b)
b = g1.board()
for move in g1.mainline_moves():
    b.push(move)
    print(b)
plot(b)
print(b)
b
b = g1.board()
for move in g1.mainline_moves():
    b.push(move)
    print(b + '\n')
for move in g1.mainline_moves():
    b.push(move)
    print(b)
    print('\n')
b = g1.board()
for move in g1.mainline_moves():
    b.push(move)
    print(b)
    print('\n')
pwd
history -f test1.py
e = chess.engine.SimpleEngine.popen_uci('E:\Tools\stockfish\stockfish-10-win\Windows\stockfish_10_x64.exe')
e.analyse?
e.analyse(b)
e.analyse(b, limit=1)
e.analyse(b, limit=chess.engine.Limit(time=0.100))
b = g1.board()
b
e.analyse(b, limit=chess.engine.Limit(time=0.100))
e.analyse(b, limit=chess.engine.Limit(time=0.500))
for move in g1.mainline_moves([:10]):
    b.push(move)
for move in g1.mainline_moves()[:10]:
    b.push(move)
b = g1.board()
moves = [x for x in g1.mainline_moves()]
moves
b
for move in moves[:10]:
    b.push(move)
b
e.analyse(b, limit=chess.engine.Limit(time=0.500))
e.analyse(b, limit=chess.engine.Limit(time=0.100))
e.analyse(b, limit=chess.engine.Limit(time=2))
e.analyse(b, limit=chess.engine.Limit(time=4))
e
chess.engine
e.analyse?
e.analyse??
e.analyse(b, limit=chess.engine.Limit(time=2), multipv=2)
e.analyse(b, limit=chess.engine.Limit(time=2), multipv=4)
rs = e.analyse(b, limit=chess.engine.Limit(time=2), multipv=4)
rs
rs[0]
rs[0]['pv']
rs[0]['pv'][0]
rs[0]['score']
sc = rs[0]['score']
sc
sc + 1
sc
sc.relative
sc.pov
sc.pov()
sc.relative
sc.relative
sc.relative.cp
sc.relative.cp + 1
history -f test1.py
b = g1.board()
e.analyse(b, limit=chess.engine.Limit(time=4))
f.close()
UT
f = open('chess_com_games_2019-09-01.pgn')
history
g1 = pgn.read_headers?
g1 = pgn.read_headers
headers = pgn.read_headers()
headers = pgn.read_headers(g1)
g1
g1 = pgn.read_game(f)
headers = pgn.read_headers(g1)
pgn
pgn.read_headers
pgn.read_headers?
headers = pgn.read_headers(f)
headers
headers = pgn.read_headers(f)
headers
import pickle
pickle?
p = pickle(headers)
str(headers)
s = str(headers)
eval(s)
headers
header.__dict__
headers.__dict__
h = pgn.read_headers(f).__dict__
h
str(h)
s = str(h)
s
h2 = eval(s)
h2
s
headers
history -f test1.py
