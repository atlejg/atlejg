import string
import math

class SimulationCounter:
   '''
   something ala ParameterStudy
   usage:
   sc = SimulationCounter()
   sc.load('dp', [1, 3, 5], [0])
   sc.load('diam', [10, 30, 50], [0])
   sc.load('rho', [700, 800, 900, 1000], [0, 1])
   sc.load('visc', [1, 2, 3, 4], [0])
   sc.load('length', arange(0,1000,10))
   n2 = []
   while sc.has_next('dp'):
      dp = sc.next('dp')
      sc.reset('diam')
      while sc.has_next('diam'):
         diam = sc.next('diam')
         sc.reset('rho')
         while sc.has_next('rho'):
            rho = sc.next('rho')
            sc.reset('visc')
            while sc.has_next('visc'):
               visc = sc.next('visc')
               sc.reset('length')
               while sc.has_next('length'):
                  length = sc.next('length')
                  n2.append(sc.str())
   '''
   def __init__(self):
      self.vars     = {}
      self.selected = {}
      self.varnames = [] # need correct order of variables
      self.index    = {}

   def load(self, varnm, vals, selected=None):
      self.vars[varnm] = vals
      self.varnames.append(varnm)
      self.reset(varnm)
      if selected: self.selected[varnm] = selected
      else       : self.selected[varnm] = range(len(vals)) # all
   
   def has_next(self, varnm):
      return self.index[varnm] < len(self.selected[varnm]) - 1

   def get(self, varnm):
      if self.index[varnm] < 0 or self.index[varnm] > len(self.vars[varnm])-1:
         raise Exception('index out of range')
      return self.vars[varnm][self.selected[varnm][self.index[varnm]]]

   def next(self, varnm):
      self.index[varnm] += 1
      return self.get(varnm)

   def reset(self, varnm):
      self.index[varnm] = -1

   def str(self):
      s = ''
      for varnm in self.varnames:
         fmt = '%%0%ii' % int(math.log10(len(self.vars[varnm])) + 1) # padding with zeros
         s += fmt % self.selected[varnm][self.index[varnm]]
      return s

   def num(self):
      return int(self.str())


class LoopCounter:
   '''
   Useful when creating simulation files based on parameter index.
   Example:
   lc = LoopCounter(3)                       # 3 is number of nested loops
   for i in range(4):
      lc.inc(1)                              # first loop
      for j in range(5):
         lc.inc(2)                           # second loop
         for k in range(6):
            lc.inc(3)                        # third loop
            fname = lc.str() + '.DATA'       # using the number as string
   '''
   def __init__(self, nloops):
      self.nloops = nloops
      self.cnt = []
      for i in range(nloops): self.cnt.append(-1)

   def inc(self, loopno):
      self.cnt[loopno-1] += 1
      if loopno < self.nloops: self.cnt[loopno] = -1

   def str(self):
      return string.join([str(int(i)) for i in self.cnt], '')

   def num(self):
      return int(self.str())

if __name__ == '__main__':
   lc = LoopCounter(3)
   n1 = []
   for i in range(4):
      lc.inc(1)
      for j in range(5):
         lc.inc(2)
         for k in range(6):
            lc.inc(3)
            n1.append(lc.num())

   sc = SimulationCounter()
   sc.load('dp', [1, 3, 5], [0])
   sc.load('diam', [10, 30, 50], [0])
   sc.load('rho', [700, 800, 900, 1000], [0, 1])
   sc.load('visc', [1, 2, 3, 4], [0])
   sc.load('length', arange(0,1000,10))
   n2 = []
   while sc.has_next('dp'):
      dp = sc.next('dp')
      sc.reset('diam')
      while sc.has_next('diam'):
         diam = sc.next('diam')
         sc.reset('rho')
         while sc.has_next('rho'):
            rho = sc.next('rho')
            sc.reset('visc')
            while sc.has_next('visc'):
               visc = sc.next('visc')
               sc.reset('length')
               while sc.has_next('length'):
                  length = sc.next('length')
                  n2.append(sc.str())

