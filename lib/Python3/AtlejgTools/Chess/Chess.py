'''
abbreviations:
   pc   : piece
   ptype: type of piece (PAWN, QUEEN, ...)
   instr: pgn instruction (ala a3, Rxg4)
   sq   : square on chess board
'''

import pgn
import pylab as pl
import sys, pdb
import AtlejgTools.Utils as UT
from Tools.sammon import sammon_plot
#pdb.set_trace()

WHITE  = 'white'
BLACK  = 'black'
PAWN   = 'pawn'
BISHOP = 'bishop'
KNIGHT = 'nait'    # a bit stupid...
ROOK   = 'rook'
QUEEN  = 'queen'
KING   = 'king'

FILES  = 'abcdefgh'
RANKS  = list(range(1, 9))

def _game_over(instr):
    if '0-1'     in instr: return True
    if '1-0'     in instr: return True
    if '1/2-1/2' in instr: return True
    return False

def _rank(pos):
    return int(pos[1])

def _prepare_instruction(instr):
    for ch in 'x+#':
        instr = instr.replace(ch, '')
#   for ch in 'bnrqko':
#      instr = instr.replace(ch, ch.upper())
    return instr

def _value(ptype):
    if ptype == PAWN:    return 1
    if ptype == BISHOP:  return 3
    if ptype == KNIGHT : return 3
    if ptype == ROOK  :  return 5
    if ptype == QUEEN :  return 9
    if ptype == KING  :  return 2 # could say its Inf, but as a playing piece, 2 is ok

def _ptype(ch):
    if ch == 'N': return KNIGHT
    if ch == 'B': return BISHOP
    if ch == 'R': return ROOK
    if ch == 'Q': return QUEEN
    if ch == 'K': return KING

class Piece(object):
    def __init__(self, color, name, ptype, pos):
        self.color = color
        self.ptype = ptype
        self.name  = name
        self.pos   = pos
        self.is_alive = True
        self.path  = Path()
        self._value = _value(ptype)
    def value(self):
        if not self.is_alive: return 0
        else                : return self._value
    def file(self):
        return self.pos[0]
    def rank(self):
        return _rank(self.pos)

class Pieces(object):
    def __init__(self, color):
        self.pawns    = []
        self.officers = []
        if color == WHITE:
            row1 = '1'
            row2 = '2'
        else:
            row1 = '8'
            row2 = '7'
        for x in FILES:
            self.pawns.append(Piece(color, '%s_%s'%(PAWN, x), PAWN, x+row2))
            if x in ('a', 'h'):
                self.officers.append(Piece(color, '%s_%s'%(ROOK, x), ROOK, x+row1))
            elif x in ('b', 'g'):
                self.officers.append(Piece(color, '%s_%s'%(KNIGHT, x), KNIGHT, x+row1))
            elif x in ('c', 'f'):
                self.officers.append(Piece(color, '%s_%s'%(BISHOP, x), BISHOP, x+row1))
            elif x == 'd':
                self.officers.append(Piece(color, QUEEN, QUEEN, x+row1))
            else:
                self.officers.append(Piece(color, KING, KING, x+row1))
        self.all = self.pawns + self.officers
    def get(self, pos):
        for pc in self.all:
            if pc.pos == pos and pc.is_alive: return pc
        return None # not found

def _ij(pos):
    # i is column (file), j is row (rank)
    i = ord(pos[0]) - ord('a')
    j = int(pos[1]) - 1
    return (i,j)

def _pos(i,j):
    return FILES[i]+str(RANKS[j])

def _dist(p1, p2):
    # p1 is (i1,j1), p2 is (i2,j2)
    return max(abs(p1[0]-p2[0]), abs(p1[0]-p2[0]))

def _color(pos):
    i, j = _ij(pos)
    if pl.mod(i+j, 2): return WHITE
    else          : return BLACK

class Square(object):
    def __init__(self, pos):
        self.pos = pos
        self.color = _color(pos)
        self.pc    = None
    def set(self, pc):
        self.pc = pc
    def get(self):
        return self.pc

class Board(object):
#
    def __init__(self):
        self.wpieces = Pieces(WHITE)
        self.bpieces = Pieces(BLACK)
        self.player  = WHITE
        self.drawno  = 0
        g = []
        for x in FILES:
            g.append([])
            for i in RANKS:
                pos = '%s%i' % (x, i)
                sq = Square(pos)
                pc = self.wpieces.get(pos)
                if not pc: pc = self.bpieces.get(pos)
                if pc: sq.pc = pc
                g[-1].append(sq)
        self.grid = g
#
    def get(self, *args):
        if len(args) == 1: i,j = _ij(args[0])  # got postion (a1 etc)
        else             : i,j = args          # got i,j
        return self.grid[i][j]
#
    def get_piece(self, *args):
        if len(args) == 1: i,j = _ij(args[0])  # got postion (a1 etc)
        else             : i,j = args          # got i,j
        return self.get(i,j).pc
#
    def move_piece(self, to_pos, pc):
        #pdb.set_trace()
        self.get(pc.pos).pc = None
        pc0 = self.get_piece(to_pos)
        if pc0: pc0.is_alive = False
        self.get(to_pos).pc = pc
        pc.pos = to_pos
        return pc # convinent
#
    def _pieces(self):
        return self.wpieces if self.player == WHITE else self.bpieces
#
    def _active_pieces(self, ptype):
        return [pc for pc in self._pieces().all if pc.ptype == ptype and pc.is_alive]
#
    def _move(self, ptype, to_pos, from_pos=None, file=None, rank=None):
        if self.drawno == STOP_AT: pdb.set_trace()
        if from_pos:
            pc = self.get_piece(_ij(from_pos))
        elif ptype == PAWN:
            # *** wrong ***
            if not file: file = to_pos[0]
            to_rank = _rank(to_pos)
            pc = None
            scaler = 1 if self.player == WHITE else -1
            for pawn in self._pieces().pawns:
                if not pawn.is_alive:     continue
                if not pawn.file() == file: continue
                rank = _rank(pawn.pos)
                if scaler*(to_rank - rank) == 1: # found it?
                    pc = pawn
                    break
                if scaler*(to_rank - rank) == 2: # pawn moved to steps?
                    pc2 = pawn
            if not pc: pc = pc2
        else:
            pcs = self._active_pieces(ptype)
            if len(pcs) == 1:
                pc = pcs[0]
            else:
                if file or rank:
                    # we have been given a clue
                    if file and file in pcs[0].pos: pc = pcs[0]
                    if file and file in pcs[1].pos: pc = pcs[1]
                    if rank and rank in pcs[0].pos: pc = pcs[0]
                    if rank and rank in pcs[1].pos: pc = pcs[1]
                else:
                    # must search among legal moves
                    done = False
                    for pc in pcs:
                        for pos in self.legal_moves(pc):
                            if pos == to_pos:
                                done = True
                                break
                        if done: break
        self.move_piece(to_pos, pc)
#
    def get_state(self):
        m = pl.zeros((8,8))
        for j in range(8,0,-1):
            for i in RANKS:
                pc = self.get_piece(i-1,j-1)
                if pc: m[i-1,j-1] = pc.value()
        return m
#
    def __str__(self):
        s = '\n   -------------------------------\n'
        for j in range(8,0,-1):
            s += '%i | ' % j
            for i in RANKS:
                pc = self.get_piece(i-1,j-1)
                if pc and pc.is_alive:
                    p = pc.ptype[0]
                    if pc.color == WHITE: p = p.upper()
                    s += p
                else:
                    s+= ' '
                s += ' | '
            s += '\n   -------------------------------\n'
        s += '    a   b   c   d   e   f   g   h \n'
        return s



#
    def is_free(self, i,j):
        return True if not self.get_piece(i,j) else False
#
    def legal_moves_old(self, pc, ij=False):
        '''
        pc: piece
        ij: return i,j-pairs (boolean)
        '''
        if not pc.is_alive: return []
        lm = [] # legal moves
        i,j = _ij(pc.pos)
        if pc.ptype == PAWN:
            # note - dont include moves where pawn takes
            if pc.color == WHITE:
                jump = 1
                init_j = 1
            else:
                jump = -1
                init_j = 6
            j2 = j + jump
            if self.is_free(i,j2): lm.append((i,j2))
            if j == init_j: # first move
                j2 += jump
                if self.is_free(i,j2): lm.append((i,j2))
        if pc.ptype == BISHOP or pc.ptype == QUEEN or pc.ptype == KING:
            # diagonal move
            stop = False
            i2, j2 = (i+1, j+1)
            while i2 < 8 and j2 < 8:
                if pc.ptype == KING and _dist((i,j), (i2,j2)) > 1: break
                pc0 = self.get_piece(i2,j2)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i2,j2))
                if stop: break
                i2 += 1; j2 += 1
            stop = False
            i2, j2 = (i+1, j-1)
            while i2 < 8 and j2 >= 0:
                if pc.ptype == KING and _dist((i,j), (i2,j2)) > 1: break
                pc0 = self.get_piece(i2,j2)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i2,j2))
                if stop: break
                i2 += 1; j2 -= 1
            stop = False
            i2, j2 = (i-1, j+1)
            while i2 >= 0 and j2 < 8:
                if pc.ptype == KING and _dist((i,j), (i2,j2)) > 1: break
                pc0 = self.get_piece(i2,j2)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i2,j2))
                if stop: break
                i2 -= 1; j2 += 1
            stop = False
            i2, j2 = (i-1, j-1)
            while i2 >= 0 and j2 >= 0:
                if pc.ptype == KING and _dist((i,j), (i2,j2)) > 1: break
                pc0 = self.get_piece(i2,j2)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i2,j2))
                if stop: break
                i2 -= 1; j2 -= 1
        if pc.ptype == ROOK or pc.ptype == QUEEN or pc.ptype == KING:
            # straigth move
            stop = False
            i2 = i+1
            while i2 < 8:
                if pc.ptype == KING and _dist((i,j), (i2,j)) > 1:  break
                pc0 = self.get_piece(i2,j)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i2,j))
                if stop: break
                i2 += 1
            stop = False
            i2 = i-1
            while i2 >= 0:
                if pc.ptype == KING and _dist((i,j), (i2,j)) > 1:  break
                pc0 = self.get_piece(i2,j)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i2,j))
                if stop: break
                i2 -= 1
            stop = False
            j2 = j+1
            while j2 < 8:
                if pc.ptype == KING and _dist((i,j), (i,j2)) > 1:  break
                pc0 = self.get_piece(i,j2)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i,j2))
                if stop: break
                j2 += 1
            stop = False
            j2 = j-1
            while j2 >= 0:
                if pc.ptype == KING and _dist((i,j), (i,j2)) > 1:  break
                pc0 = self.get_piece(i,j2)
                if pc0:
                    if   pc0.color == self.player: break
                    else                         : stop = True
                lm.append((i,j2))
                if stop: break
                j2 -= 1
        if pc.ptype == KNIGHT:
            jumps = [(1,2), (-1,2), (-2,1), (-2,-1), (-1,-2), (1,-2), (2,-1), (2,1)]
            for jump in jumps:
                i2, j2 = (i+jump[0], j+jump[1])
                if not (0 <= i2 < 8) or not (0 <= j2 < 8): continue
                pc0 = self.get_piece(i2,j2)
                if pc0 and pc0.color == self.player:       continue
                lm.append((i2,j2))
        if ij: return lm
        else : return [_pos(x[0],x[1]) for x in lm]
    def _legal_moves(self, pc, vec):
        i,j = _ij(pc.pos)
        lm = []
        for scaler in (-1, 1):
            for v in vec:
                n = 1
                stop = False
                while True:
                    i2 = i + n*scaler*v[0]
                    if not 0 <= i2 < 8:                                break
                    j2 = j + n*scaler*v[1]
                    if not 0 <= j2 < 8:                                break
                    if pc.ptype == KING and _dist((i,j), (i1,j2)) > 1: break
                    pc0 = self.get_piece(i2,j2)
                    if pc0:
                        if pc0.color == self.player: break  # cannot take own piece
                        else                       : stop = True
                    lm.append((i2,j2))
                    if stop: break                         # cannot move behind foe's piece
                    n += 1
        return lm
    def legal_moves(self, pc, ij=False):
        '''
        pc: piece
        ij: return i,j-pairs (boolean)
        '''
        if not pc.is_alive: return []
        lm = [] # legal moves
        i,j = _ij(pc.pos)
        if pc.ptype == PAWN:
            # note - dont include moves where pawn takes
            if pc.color == WHITE:
                jump = 1
                init_j = 1
            else:
                jump = -1
                init_j = 6
            j2 = j + jump
            if self.is_free(i,j2): lm.append((i,j2))
            if j == init_j: # first move
                j2 += jump
                if self.is_free(i,j2): lm.append((i,j2))
        else:
            if pc.ptype == BISHOP or pc.ptype == QUEEN or pc.ptype == KING:
                # diagonal move
                vec = ((1,1), (1,-1))
                lm.extend(self._legal_moves(pc, vec))
            if pc.ptype == ROOK or pc.ptype == QUEEN or pc.ptype == KING:
                # straigth move
                vec = ((0,1), (1,0))
                lm.extend(self._legal_moves(pc, vec))
            if pc.ptype == KNIGHT:
                jumps = [(1,2), (-1,2), (-2,1), (-2,-1), (-1,-2), (1,-2), (2,-1), (2,1)]
                for jump in jumps:
                    i2, j2 = (i+jump[0], j+jump[1])
                    if not (0 <= i2 < 8) or not (0 <= j2 < 8): continue
                    pc0 = self.get_piece(i2,j2)
                    if pc0 and pc0.color == self.player:       continue
                    lm.append((i2,j2))
        if ij: return lm
        else : return [_pos(x[0],x[1]) for x in lm]
#
    def _find(self, ptype, to_pos, i=None, j=None):
        pcs = self._active_pieces(ptype)
        if len(pcs) == 1: return pcs[0]
        for pc in pcs:
            i1, i1 = _ij(pc.pos)
            if i and (i == i1) or j and (j==j1): return pc
            if to_pos in self.legal_moves(pc)         : return pc
#
    def castle(self,instr):
        king = self._pieces().officers[4]
        rank = 1 if self.player == WHITE else 8
        if instr == 'O-O':
            # short castling
            rook = self._pieces().officers[7]
            self.move_piece('g%i'%rank, king)
            self.move_piece('f%i'%rank, rook)
        else:
            # long castling
            rook = self._pieces().officers[0]
            self.move_piece('c%i'%rank, king)
            self.move_piece('d%i'%rank, rook)
#
    def make_move(self, instr):
        '''
        instr: pgn short hand instruction (like Rac8, e3 ...)
        '''
        # make sure we shift player
        self.drawno += 1
        if pl.mod(self.drawno, 2): self.player = WHITE
        else                     : self.player = BLACK
        #
        # cleanup
        instr = _prepare_instruction(instr)
        #
        # some special cases
        if 'O' in instr:
            self.castle(instr)
            return
        if _game_over(instr):
            return
        #
        # analyze instruction
        file = None
        to_pos = instr[-2:]
        if len(instr) > 2:
            x = instr[0]
            if x in FILES:
                ptype = PAWN
                file  = x
            else:
                ptype = _ptype(x)
        else: ptype = PAWN
        #
        # now it's time to move
        if len(instr) <= 3: self._move(ptype, to_pos, file=file)
        else :
            x = instr[1]
            if x in FILES:     self._move(ptype, to_pos, file=x)
            else         :     self._move(ptype, to_pos, rank=x)

class Path(object):
    def __init__(self):
        self.hist = []   # list of Square-s
    def add(sq):
        self.hist.append(sq)

class Game(object):
    def __init__(self):
        self.board   = Board()

if __name__ == '__main__':
    STOP_AT = 999*2 - 0 # for debugging
    DEBUG   = False
    figsize = pl.rcParams['figure.figsize']
    pl.rcParams['figure.figsize'] = (8,8)
    pgnfiles = UT.glob(sys.argv[1])
    pgnfiles.sort()
    for pgnfile in pgnfiles:
        print(pgnfile)
        pgn_instructions = pgn.read_pgn(pgnfile)
        g = Game()
        i = 1
        m = []
        for instructions in pgn_instructions:
            if DEBUG: print(g.board)
            print(i, instructions)
            g.board.make_move(instructions[0])
            if len(instructions) >= 2:
                g.board.make_move(instructions[1])
            m.append(g.board.get_state().flatten())
            i += 1
        print(g.board)
        m = pl.array(m)
        sammon_plot(m, list(range(len(m))))
        roundno = UT.grep_column(pgnfile, 'Round ', 2, False)[0][1:-2]
        result  = UT.grep_column(pgnfile, 'Result ', 2, False)[0][1:-2]
        white   = UT.grep_column(pgnfile, 'White ', 2, False)[0][1:-1]
        black   = UT.grep_column(pgnfile, 'Black ', 2, False)[0][1:-1]
        titl = 'Round %s. %s vs. %s. Result: %s' % (roundno, white,black,result)
        pl.title(titl)
        pl.axis('equal')
    pl.show()
    pl.rcParams['figure.figsize'] = figsize
