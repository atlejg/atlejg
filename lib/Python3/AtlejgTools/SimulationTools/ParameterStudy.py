'''
Atle J. Gyllensten
agy@statoil.com

General info:
This module is for creating/running/analazing simulations (e.g. Fluent simulations) with multiple varying parameters.

Technical info:
What could be a bit confusing is that 'selected' indexes are positive numbers ala Matlab (i.e. first element in an
array is counted as being #1 etc), while the internal index starts at 0. this means we have a lot of +/- 1 in the code:-(

Here is a typical example of a ini-file to be used for varying simulation parameters (this is a 'dynsim' example (troll box model))
   templ       = A.templ
   prefix      = A
   overwrite   = YES
   suffix      = .csv

   [branches]
   linkto      = None
   values      = ['mainbore', 'lateral_1', 'lateral_2']
   unit        = -
   selected    = ALL
   mnemonic    = compl
   format      = %s
   description = []
   replace     = $BRANCH$
   [owcfile]
   linkto      = branches
   values      = ['owc_m.txt', 'owc_1.txt', 'owc_2.txt']
   unit        = -
   selected    = ALL
   mnemonic    = file
   format      = %s
   description = []
   replace     = $OWC_FILE$

   [completion]
   linkto      = None
   values      = ['Netool', '1:1', '1:2', '1:3', '1:4']
   unit        = -
   selected    = ALL
   mnemonic    = compl
   format      = %s
   description = []
   replace     = $COMPL$
   [honour_netool]
   linkto      = completion
   values      = ['True', 'False', 'False', 'False', 'False']
   unit        = -
   selected    = ALL
   mnemonic    = honour
   format      = %s
   description = []
   replace     = $HONOUR$

   [oil_column]
   linkto      = None
   values      = [1, 3, 5]
   unit        = -
   selected    = ALL
   mnemonic    = pay
   format      = %.1f
   description = []
   replace     = $OILCOLUMN$




Nov 09 : Version 0.1

Nov 10 : found a bug. when using linkto, this must be the last entries of the ini-file.
         dont know why (this means only one set of linked values may be used:-( )
         --> FIXED!
'''

import re, sys, os, math, configparser, copy
import AtlejgTools.LogTools as LogTools
import AtlejgTools.Utils as Utils
from . import UnitConversion
import pdb
import pylab # useful to have it available when reading parameters with eval

_log = LogTools.Log()
_init_value = -1

def set_loglevel(level):
    _log.set_level(level)

class Struct : pass

class Counter(object):
    '''
    Need a counter that can be shared  between linked Param's
    '''

    class NoSuchIndex(Exception): pass

    def __init__(self,start_val=_init_value):
        self._n = start_val

    def prev(self):
        self._n -= 1
        if self._n == _init_value:
            raise NoSuchIndex()

    def __next__(self):
        self._n += 1

    def reset(self):
        self._n = _init_value

    def get(self):
        return self._n

    def set(self,n):
        self._n = n

def _indent(level):
    return '   ' * level

def _formatted(val,fmt):
    formatstr = "'%s' %% %s" % (fmt, str(val))
    return eval(formatstr)


class Param(object):

    ALL = 'ALL';

    def __init__(self,name,vals,unit,sel,mnemonic,format,linkto,descr, replace) :
        self.vals = vals
        self._nm   = name
        self._unit= unit
        self._mnem = mnemonic
        self._fmt  = format
        self._cnt  = Counter()
        self.selected(sel)
        if linkto == 'None': self._linkto = None
        else               : self._linkto = linkto
        self._descr= descr
        self._linkfrom = []
        self._replace = replace
        # number of digits (based on max possible values)
        self._ndigits = int(math.log(len(self.vals))/math.log(10.0)) + 1

    def next(self,use_si=False):
        next(self._cnt)
        _log.debug('next: self._cnt= ' + str(self._cnt.get()))
        return self.value(use_si)

    def prev(self,use_si=False):
        self._cnt.prev()
        _log.debug('prev: self._cnt= ' + str(self._cnt.get()))
        return self.value(use_si)

    def unit(self,use_si=False):
        if use_si: return UnitConversion.si_unit(self._unit)
        else:      return self._unit

    def index(self):
        '''
        index is in [0,N)
        '''
        return self._sel[self._cnt.get()] - 1

    def counter(self,n=None):
        '''
        get or set method
        '''
        if n != None: self._cnt.set(n)
        return self._cnt.get()

    def value(self,use_si=False, index=-1):
        if index < 0: index = self.index()
        v = self.vals[index]
        if use_si: return UnitConversion.to_si(self._unit,v)
        else:      return v

    def get_selected_values(self,use_si=False):
        return [self.value(index=i-1,use_si=use_si) for i in self._sel]

    def description(self,id):
        if self._descr: descr = self._descr[id-1]
        else          : descr = self._mnem + "=" + _formatted(self.value(index=id-1),self._fmt) + self.unit()
        return descr

    def has_next(self):
        _log.debug('has_next: param: %s self._cnt= %i' % (self._nm, self._cnt.get()))
        if self._cnt.get() < len(self._sel) - 1:
            return True
        else:
            return False

    def has_prev(self):
        if self._cnt.get() > _init_value:
            return True
        else:
            return False

    def reset(self):
        self._cnt.reset()

    def is_reset(self):
        return self._cnt.get() == _init_value

    def selected(self,sel=None):
        '''
        get or set method
        '''
        if sel != None:
            if sel == Param.ALL:
                self._sel = list(range(1,len(self.vals)+1))
            else:
                self._sel  = sel
            self.reset()
        return self._sel

    def __str__(self):
        v_str = ''
        for v in self.vals:
            v_str += (str(v) + ' ')
        return self._mnem + ' [' + v_str[0:-1] + '] ' + self._unit

    def id(self):
        return self._sel[self._cnt.get()]

    def link(self,pm):
        if isinstance(pm,ParamCollection): pm = pm.get(self._linkto)
        # Param's will share conuter and selection
        self._sel = pm._sel
        self._cnt = pm._cnt
        self._linkto = pm # more convinient to have the reference directly
        pm._linkfrom.append(self)

    def max(self, use_si=False) :
        return max(self.get_selected_values(use_si))

    def min(self, use_si=False) :
        return min(self.get_selected_values(use_si))

# end of class Param

class ParamRealization(object):
    '''
    this class defines an object that holds one value for each paramter,
     i.e. parameter values for one case
     '''

    class NoSuchParamNm(Exception):
        def __init__(self): pass

    def __init__(self,pm_coll):
        self._pms     = {}
        self._indices = {} # index directly into Params.vals.
        self._pm_coll = pm_coll

    def copy(self):
        pr = ParamRealization(self._pm_coll)
        pr._indices  = self._indices.copy()
        pr._pms = self._pms.copy()
        for key in list(self.__dict__.keys()):
            if (key != '_indices' and key != '_pms'):
                pr.__dict__[key] = copy.copy(self.__dict__[key])
        return pr

    def set(self,pm,index):
        self._pms[pm._nm]     = pm
        self._indices[pm._nm] = index

    def id(self, numeric=True, skip_these=[], single=None):
        s = ''
        for pm in self._pm_coll._params:
            if pm._nm in skip_these       : continue
            if single and pm._nm != single: continue
            if numeric:
                if pm._linkto: continue # numeric id represented by the linkto Param
                s += _formatted(self._indices[pm._nm]+1, '%%0%ii' % pm._ndigits)
            else:
                if pm._descr:
                    s += pm._descr[self._indices[pm._nm]] + '_'
                elif pm._linkto and pm._linkto._descr: continue  # the linked param handles it
                else:
                    formatstr = "'%s=%s%s_' %% pm.value(False,index=self._indices[pm._nm])" % (pm._mnem,pm._fmt ,pm.unit(False))
                    _log.debug('ParamRealization.id: formatstr = ' + formatstr)
                    s = s + eval(formatstr)
        if not numeric:
            s = s[:-1] # remove last delimiter
        return s

    def __str__(self):
        return self.id(numeric=False)

    def get_parameter(self,nm):
        '''
        finds parameter based on name or mnemonic
        '''
        for pm in list(self._pms.values()):
            if pm._nm == nm or pm._mnem == nm:
                return pm
        raise NoSuchParamNm()

    def get(self,nm,use_si=False,as_string=False):
        '''
        finds value based on name or mnemonic
        '''
        pm = self.get_parameter(nm)
        val = pm.value(use_si,index=self._indices[pm._nm])
        if type(val) is str: return val  # dont do anything with strings
        if as_string:
            if not use_si : return _formatted(val,pm._fmt)
            else          : return _formatted(val, '%e') # avoid problems ...
        else:
            return val

    def get_index(self,nm):
        '''
        finds index based on name or mnemonic
        note: this is the same index as in the 'selected' array, i.e. starts with 1 (not 0)
        '''
        id = None
        for pm in list(self._pms.values()):
            if pm._nm == nm or pm._mnem == nm:
                id = self._indices[pm._nm] + 1
                break
        if id == None: raise NoSuchParamNm()
        return id

    def reset(self, pm) :
        '''
        useful if were using one ParamRealization to loop parameter values.
        should be used together with next_index() since it sets index to one lower
        than the first selected index
        '''
        if isinstance(pm, str) : pm = self.get_parameter(pm) # are we given a name?
        self._indices[pm._nm] = pm._sel[0] - 2 # index prior to first selected

    def next_index(self, pm, shift=1) :
        '''
        useful if were using one ParamRealization to loop parameter values.
        shift could be negative
        '''
        if isinstance(pm, str) : pm = self.get_parameter(pm) # are we given a name?
        # mixing indices and selected values - shaky stuff...
        i = 0
        for sel in pm._sel :
            if self._indices[pm._nm] + shift + 1 == sel : break
            i += 1
        #
        if 0 <= i < len(pm._sel) :
            next_index = pm._sel[i] - 1
        else :
            next_index = -1
        _log.debug('next_index= %i' % next_index)
        #
        return next_index

    def set_index(self, pm, index) :
        '''
        useful if were using one ParamRealization to loop parameter values.
        '''
        if isinstance(pm, str) : pm = self.get_parameter(pm) # are we given a name?
        self._indices[pm._nm] = index

    def shift_index(self, pm, shift=1) :
        '''
        useful if were using one ParamRealization to loop parameter values.
        shift could be negative
        '''
        if isinstance(pm, str) : pm = self.get_parameter(pm) # are we given a name?
        self._indices[pm._nm] += shift

    def get_replace_array(self, use_si=True):
        '''
        replacements put into array which can be fed to Utils.replace_in_file
        '''
        replace = []
        for pm in list(self._pms.values()) :
            replace.append([pm._replace, self.get(pm._nm, use_si, True)])
        return replace

# end of class ParamRealization(object)


class ParamCollection(object):

    def __init__(self):
        self._params = [] # Params
        self._prs    = [] # ParamRealizations

    def append(self,pm):
        # only append Params with selected values
        if pm._sel:
            self._params.append(pm)

    def get(self,nm):
        for pm in self._params:
            if pm._nm == nm or pm._mnem == nm:
                return pm

    def loop_all(self,cnst,opts,init_func,main_func,exit_func):

        NL         = '\n'

        # we're gonna build dynamic code
        codestr = ''
        if init_func != None: codestr = 'init_func(cnst,opts,self._params)' + NL

        pr = ParamRealization(self)

        n1 = -1 # for indenting
        n2 = -1 # for indexing parameters
        for pm in self._params:
            _log.debug('ParamCollection::loop_all param= %s' % pm._nm)
            if pm._linkto: continue
            n1 += 1
            n2 += 1
            codestr += _indent(n1+0) + 'self._params[' + str(n2) + '].reset()' + NL
            codestr += _indent(n1+0) + 'while self._params[' + str(n2) + '].has_next():' + NL
            codestr += _indent(n1+1) + 'self._params[' + str(n2) + '].next()' + NL
            codestr += _indent(n1+1) + 'pr.set(self._params[' + str(n2) + '], self._params[' + str(n2) + '].index())' + NL
            for i in range(0, len(pm._linkfrom)):
                codestr += _indent(n1+1) + 'link = self._params[' + str(n2) + ']._linkfrom[' + str(i) + ']' + NL
                codestr += _indent(n1+1) + 'pr.set(link, link.index())' + NL
            n2 += len(pm._linkfrom)

        codestr += _indent(n1+1) + 'main_func(cnst,opts,pr)' + NL
        codestr += _indent(n1+1) + 'self._prs.append(pr.copy())' + NL
        if exit_func != None: codestr += 'exit_func(cnst,opts,self)' + NL

        _log.debug('ParamCollection::loop_all codestr=' + NL + codestr)

        # done building code, now execute
        exec(codestr)

# end of class ParamCollection(object)

def default_fluent_func(cnst, opts, pr) :
    """
    useful when all you want to do is the standard stuff, i.e.
    replace __XXX__ expressions in bc and runfiles and start fluent.
    assumes you have __RUNID__ in runfile and bcfile.
    the cnst object must have the following attributes:
    * fluent_cmd
    * geom_id
    * run_in_bg
    the opts object must have the following attributes:
    * run_id
    """
    #
    # filenames
    #
    casedir = cnst.geom_id + pr.id()
    _log.debug('casedir = ' + casedir)
    templ_bc  = os.path.join(cnst.geom_id,'%i.bc'  % (opts.run_id))
    templ_jou = os.path.join(cnst.geom_id,'%i.jou' % (opts.run_id))
    bcfile  = os.path.join(casedir, '%i.bc'  % opts.run_id)
    runfile = os.path.join(casedir, '%i.jou' % opts.run_id)
    #
    # run new-command
    #
    cmd = "new -c %s.cas case %s" % (cnst.geom_id, casedir)
    _log.info('default_fluent_func: cmd = ' + cmd)
    Utils.tcsh(cmd)
    #
    # create bc-file and runfile
    #
    replace = []
    replace.append(('__RUNID__', str(opts.run_id)))
    replace.extend(pr.get_replace_array())
    # copy template files, replacing __*__ expressions
    Utils.replace_in_file(replace, templ_bc,  bcfile)
    Utils.replace_in_file(replace, templ_jou, runfile)
    #
    # run fluent
    #
    cmd = 'cd %s; %s -i %s >&! %i.tra' % (casedir, cnst.fluent_cmd, os.path.basename(runfile), opts.run_id)
    if int(cnst.run_in_bg) : cmd += ' &'
    _log.info('default_fluent_func: cmd = ' + cmd)
    Utils.tcsh(cmd)

def get_section_names(inifile):
    # cannot use RawConfigParser.sections() (=cf.sections()) because it does not preserve order

    file = open(inifile)

    secs = []
    for line in file:
        match = re.search('^\[\s*(\w+)\s*\]', line)
        if match:
            sec = match.group(1)
            _log.debug('sec = ' + sec)
            if not sec == 'constants': secs.append(sec)

    file.close()
    _log.debug(secs)
    return secs

def new_param(cf, section_nm):
    '''
    input : cf = ConfigParser
            section_nm = name of section in ini-file
    creates Param from a section in an ini-file
    '''
    linkto = cf.get(section_nm,'linkto')
    vals = eval(cf.get(section_nm,'values'))
    sel  = cf.get(section_nm,'selected')
    unit = cf.get(section_nm,'unit')
    mnemonic = cf.get(section_nm,'mnemonic')
    format = cf.get(section_nm,'format')
    descr  = eval(cf.get(section_nm,'description'))
    replace = cf.get(section_nm,'replace')

    if sel == 'ALL':
        sel = Param.ALL
    else:
        sel = eval(sel)

    return Param(section_nm,vals,unit,sel,mnemonic,format,linkto,descr,replace)

def create_inifile(filenm) :
    inifile = open(filenm, 'w')
    inifile.write("""\
 [constants]
 c1          = 1
 [Parameter_1]
 linkto      = None
 values      = []
 unit        = -
 selected    = []
 mnemonic    = p1
 format      = -
 description = []
 replace     = ____
    """)
    inifile.close()
    print("create_inifile : created inifile " + filenm)

def get_constants(inifile):
    '''
    nb! in the inifile, the section names and the options must be all lowercase!!!
    nb2! : all constants are read as strings and must be handled as such
    '''
    cf = configparser.ConfigParser()
    cf.read(inifile)

    cnst = Struct()

    for opt in cf.options('constants'):
        str = "cnst.%s = cf.get('constants','%s')" % (opt, opt)
        exec(str)

    _log.debug(cnst)
    return cnst

def get_ParamCollection(inifile,skip_these=[]):
    '''
    nb! in the inifile, the options must be all lowercase!!!
    '''

    cf = configparser.ConfigParser()
    cf.read(inifile)
    secs = get_section_names(inifile)

    _log.debug(cf.__dict__)

    pm_coll = ParamCollection()

    for sec in secs:
        if sec in skip_these: continue
        pm = new_param(cf,sec)
        if pm._linkto: pm.link(pm_coll)
        pm_coll.append(pm)

    return pm_coll

def run(inifile,opts,init_func,main_func,exit_func,skip_these=[]):
    cnst   = get_constants(inifile)
    pm_coll = get_ParamCollection(inifile,skip_these)
    pm_coll.loop_all(cnst,opts,init_func,main_func,exit_func)
    return pm_coll # could be useful

def create_files(buildfile, verbose=True, use_si=False):
    '''
    useful function for quickly creating files based on a template file.
    buildfile is an ini-file that should look like this:
 [constants]
 templ       = a.TEMPL
 prefix      = A
 overwrite   = NO
 suffix      = .DATA
 [licd]
 linkto      = None
 values      = [1,2]
 unit        = -
 selected    = ALL
 mnemonic    = l
 format      = %i
 description = []
 replace     = _LICD_
    the template file (a.TEMPL) then needs to have the token _LICD_ inside (which
    will be replaced).
    '''
    if not os.path.exists(buildfile):
        print('%s does not exist' % buildfile)
        return
    files = []
    def main_func(cnst, opts, pr):
        fname = cnst.prefix + pr.id() + cnst.suffix
        overwrite = ('Y' in cnst.overwrite.upper())
        if overwrite or (not os.path.exists(fname)):
            if verbose: print('creating file', fname)
            Utils.replace_in_file(pr.get_replace_array(use_si), cnst.templ, fname)
            files.append(fname)
    pm_coll = run(buildfile, None, None, main_func, None)
    return files, pm_coll

# TESTING CODE
if __name__ == "__main__":

    def main_func(cnst,opts,pr)      : print(pr.id(), pr.id(numeric=False,skip_these=['SlotSide','SlotStart']))
    def exit_func(cnst,opts,pm_coll) : pass
    def init_func(cnst,opts,params)  : pass

    inifile = sys.argv[1]

    if len(sys.argv) == 3:
        set_loglevel(int(sys.argv[2]))
    else:
        set_loglevel(LogTools.LOG_LEVEL_WARN)

    run(inifile,None,init_func,main_func,exit_func)

    #print _log.get_level()
