import AtlejgTools.LogTools as LogTools
import sys

_log = LogTools.Log(LogTools.LOG_LEVEL_WARN)

_MODE_TO_SI   = 1
_MODE_FROM_SI = 2

ONE     = 1.0
BAR     = 100000.0
ATM     = 101315.0
PSI     = 6.8948e3
HOUR    = 3600.0
MINUTE  = 60.0
DAY     = HOUR*24
WEEK    = DAY*7
MONTH   = DAY*30
YEAR    = DAY*365.24
LITRE   = 1e-03
ML      = LITRE/1000.
CM      = 1e-02
MM      = 1e-03
INCH    = 2.54e-02
M       = ONE
CP      = 1e-03
DARCY   = 9.869233e-13
DEG_ZERO= 273.15
FEET    = 0.3048
BARREL  = 158.987295/1000. # m3

# assuming conversion is done via y[si] = ax[non-si] + b
_cnv = {#          a          b          si-unit
      'bar'    : [BAR,        0.,         'pa'],
      'atm'    : [ATM,        0.,         'pa'],
      'psi'    : [PSI,        0.,         'pa'],
      'hour'   : [HOUR,       0.,         's'],
      'day'    : [DAY,        0.,         's'],
      'week'   : [WEEK,       0.,         's'],
      'year'   : [YEAR,       0.,         's'],
      'litre'  : [LITRE,      0.,         'm3'],
      'ml'     : [ML,         0.,         'm3'],
      'l/h'    : [LITRE/HOUR, 0.,         'm3/s'],
      'cp'     : [CP,         0.,         'kg/m-s'],
      'mm'     : [MM,         0.,         'm'],
      'cm'     : [CM,         0.,         'm'],
      'darcy'  : [DARCY,      0.,         'm2'],
      'inch'   : [INCH,       0.,         'm'],
      'feet'   : [FEET,       0.,         'm'],
      'barrel' : [BARREL,     0.,         'm3'],
      'degc'   : [ONE,        DEG_ZERO,   'degK'],
      'degf'   : [5./9.,      459.67*5/9.,'degK'],
      # si units follows
      'pa'     : [ONE,        0.,         'pa' ],
      's'      : [ONE,        0.,         's'  ],
      'm3'     : [ONE,        0.,         'm3' ],
      'm3/s'   : [ONE,        0.,         'm3/s'],
      'm'      : [ONE,        0.,         'm'  ],
      'm2'     : [ONE,        0.,         'm2' ],
      '-'      : [ONE,        0.,         '-'  ],
      'kg/m3'  : [ONE,        0.,         'kg/m3']
      }

def _convert(mode, unit, val):
   unit = unit.lower()  # case insenstitive
   if not _cnv.has_key(unit): raise Exception('no such unit: %s' % unit)
   if _cnv[unit][2] == unit: return val  # nothing to do
   if mode == _MODE_TO_SI:
      return _cnv[unit][0]*val + _cnv[unit][1]
   else:
      return (val - _cnv[unit][1]) / _cnv[unit][0]

def to_si(unit,val=ONE):
   return _convert(_MODE_TO_SI,unit,val)

def from_si(unit,val=ONE):
   return _convert(_MODE_FROM_SI,unit,val)

def si_unit(unit):
   unit = unit.lower()  # case insenstitive
   if not _cnv.has_key(unit): raise Exception('no such unit: %s' % unit)
   return _cnv[unit][2]

def oil_api(val, to_si = True):
   # ref: http://www.engineeringtoolbox.com/api-gravity-d_1212.html
   #
   rho_w = 983.2     # water at 60degC
   if to_si:
      sg = 141.5 / (val+131.5)
      return sg*rho_w
   else:
      sg = val/rho_w # specific gravity, relative to water 
      return (141.5/sg) - 131.5

# aliases. yes, it's the best way of naming them..
to = from_si
si = to_si

# test code
if __name__ == "__main__":

   print oil_api(float(sys.argv[1]))
