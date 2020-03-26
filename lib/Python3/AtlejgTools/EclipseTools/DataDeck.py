import os, inspect
import pylab as pl
import pdb
import AtlejgTools.Utils              as UT

# end-of-record for Eclipse input
DELIM       = '/'

# the main sections of the DATA-file
SECTIONS    = ['RUNSPEC', 'GRID', 'EDIT', 'PROPS', 'REGIONS', 'SOLUTION', 'SUMMARY', 'SCHEDULE', 'OPTIMIZE']

# some keywords (like EQUALS) are special because they are followed by lines with other keywords (like PERMX)
# as their first entry. this is problematic and must be handled specially
SPECIAL_KWS = ['EQUALS', 'ADD', 'COPY', 'MULTIPLY']

MONTH_MAP    = {'JAN':1, 'FEB':2, 'MAR':3, 'APR':4, 'MAY':5, 'JUN':6, 'JUL':7, 'AUG':8, 'SEP':9, 'OCT':10, 'NOV':11, 'DEC':12} # useful when converting dates

class DataDeck:
    '''
    atle j. gyllensten
    july 2015
    #
    this is a somewhat advanced and generic Eclipse data-file reader.
    it uses the (available) Eclipse keywords for splitting the data-file into chunks of data.
    then each chunk is splitted into (possibly) multiple lines by the DELIM separator.
    (it is important to realize that Eclipse input does not really care too much about the \\n (newline),
    even if humans will typically organize the data file according to this)
    it also handles default values by giving you a 1* (if raw data is asked for) or a NaN (if numeric value is asked for)
    #
    note1: it reads the existing keywords from a text-file. this file must be found in *this* directory and has the format eclipse_keywords.<version> where
      version is something like 2014.2. this file can be made by doing a pdf2asc (or pdf2ps + ps2ascii) on the Eclipse ref-manual. in this ascii-file you will find an index
      of all keywords. delete everything in front and after this index, and then do some grep'ing:
      > grep -v Example ecl.asc | grep -v ECLIPSE | g -v Notes | g -v Legal | g '^[A-Z][A-Z]' | awk '{print $1}' > eclipse_keywords.2014.2
      do a manual inspection to make sure this file looks ok.
      if there is a keyword in the DATA-file that is not found in the list of existing keywords, this code will probably fail big time.
      you are allowed to add any keyword by using the xtra_kws parameter found in the constructor.
    note2 : it does not follow INCLUDE-files, so the keyword you are looking for must be found in the provided file.
    note3: one technicality which makes it more difficult to handle the keywords generically, is that some of them, like COMPDAT,
      can have multiple wells inside one instance of the keyword, while others, like WELSEGS only allows one well per keyword (which
      implies multiple keywords). this is handled internally by adding a counter to each keyword, indicating how many times it has
      been found (like COMPDAT-1, WELSEGS-1, WELSEGS-2 ...)
    note4: another problem is that you are allowed to put several records in one line (separated by DELIM) and also that you are
      allowed to write comments after the DELIM. use remove_after_delim if you have one line per record and comments after DELIM.
    note5: for the SPECIAL_KWS (like EQUALS) it assumes the format is like this (with the ending '/' on a separate line):
       EQUALS
        PERMX 1000 1 11 1 19 1 1 /
        PORO  0.3  /
       /
    '''
    ALL_KEYWORDS = [] # all existing Eclipse keywords
#
    def __init__(self, fname, remove_after_delim=False, version='2016.1', xtra_kws=[]):
        '''
        fname   : the DATA or INCLUDE-file to read
        version : version of Eclipse keywords to use
        xtra_kws: list of new keywords (not found in the file with existing keywords)
        '''
        self.set_keywords(version, xtra_kws)
        self._chunks     = {}                   # one chunk for each keyword - more or less
        self._chunks0    = []                   # just a list of chunks in the original order. potentially useful...
        self._kw_counter = {}                   # must handle multiple keyword-entries (like WELSEGS)
        self.fname       = fname
        key = ''                                # current hash key
        #pdb.set_trace()
        for line in open(fname):
            line = line.strip()                  # 2016-10-26: changed from 'rstrip'
            if not line             : continue   # skip empty lines
            if line.startswith('--'): continue   # skip comments
            # need to treat specials separately
            if self._dekey(key) in SPECIAL_KWS:
                if line.startswith(DELIM): key = ''
                else                     : self._chunks[key].append(line)
                continue
            w = line[:8].rstrip()                # look for keywords
            # usually we get here
            if w in self.ALL_KEYWORDS:
                key = self._make_key(w)
                self._chunks[key] = []
                self._chunks0.append((key,self._chunks[key]))
            else:
                self._chunks[key].append(line)
        self._clean(remove_after_delim)
        self.get = self.get_values              # just an alias
#
    def set_keywords(self, version, xtra_kws, force_update=False):
        '''
        reads all existing Eclipse keywords into ALL_KEYWORDS, from file.
        can be forced to re-read, but usually it is done only once.
        '''
        if self.ALL_KEYWORDS and not force_update: return  # dont redo
        self.ALL_KEYWORDS = SECTIONS  # section keywords must be included
        self.ALL_KEYWORDS += xtra_kws
        dirname = os.path.dirname(inspect.getfile(self.__init__))
        fname = '%s/eclipse_keywords.%s' % (dirname, version)
        if os.path.exists(fname):
            self.ALL_KEYWORDS += [x.rstrip() for x in open(fname).readlines()]
            #print 'read keywords from file', fname
        else:
            print('cannot find file', fname)
            print('here is what i found:', UT.glob('%s/eclipse_keywords.*'%dirname))
#
    def _make_key(self, kw):
        # make keyword on the format COMPDAT-1 etc
        if kw in self._kw_counter: self._kw_counter[kw] += 1
        else                           : self._kw_counter[kw] = 1
        return '%s-%i' % (kw, self._kw_counter[kw])
#
    def _dekey(self, kw):
        # get back the original keyword from the format obtained from _make_key
        # - eg COMPDAT-1 => COMPDAT
        return kw.split('-')[0]
#
    def _clean(self, remove_after_delim):
        for key in list(self._chunks.keys()):
            chunk = self._chunks[key]
            if not chunk: continue                # some keywords take no parameters
            s = ''
            for line in chunk:
                line = line.strip()
                if line == DELIM: continue
                # remove inline comments
                pos = line.find('--')
                if pos > 0: line = line[:pos]
                # remove comments after DELIM - if asked for
                if remove_after_delim:
                    line = line.split(DELIM)[0] + ' ' + DELIM  # put DELIM back at the end
                # replace special characters with empty strings
                line = line.replace("'", "").replace("\\","").replace('\t',' ')
                # keep this
                s += ' ' + line
            if s and s[-1] == DELIM: s = s[:-1]     # it usually ends with a DELIM
            self._chunks[key] = [x.replace(DELIM,'').strip() for x in s.split(DELIM)]
#
    def _chunk(self, kw, indx):
        return self._chunks['%s-%s'%(kw,indx)]
#
    def _raw_data(self, chunk):
        raw = []           # to be returned
        for line in chunk:
            rec = line.split()
            s = ''
            for v in rec:
                if not v[0].isdigit():
                    s += v+' '
                    continue
                if not '*' in v:
                    s += v+' '
                    continue
                # need to handle things like 3* or 3*250
                mult, val = v.split('*')
                for n in range(int(mult)):
                    if val: s += val+' '
                    else  : s += '1* '
            if s: raw.append(s.split())
        return raw
#
    def has_kw(self, kw):
        return kw in self._kw_counter
#
    def get_raw_data(self, kw, identif='', get_all=False):
        '''
        returns the parameters for the given keyword in text form as a
        list of lists. even if this method is called 'get_raw_data', data
        is not raw-raw - they have been cleaned (comments removed, default values
        expanded etc.)
        if mulitple entries exists, like WELSEGS, you need to provide a
        unique identifier too search for or use get_all to actually get all
        #
        kw     : keyword (like WELSEGS)
        identif: unique identifier (like well name)
        get_all: if multiple entries exists, get everything
                 (but then you cannot use raw2values so be prepared to do some parsing of the result..)
        '''
        # - handle all special cases first
        if not self.has_kw(kw):
            # -no chunks
            print('no such keyword: %s' % kw)
            return []
        if self._kw_counter[kw] > 1 and get_all:
            # -get all chunks
            raws = []
            for i in range(1, self._kw_counter[kw]+1):
                raws.extend(self._raw_data(self._chunk(kw, i)))
            return raws
        if self._kw_counter[kw] > 1 and not identif:
            # -more than 1 chunk, need identifier
            print('multiple data exists - need an identifier')
            return []
        if self._kw_counter[kw] > 1:
            # -more than 1 chunk. has identifier
            foundit = False
            for i in range(1, self._kw_counter[kw]+1):
                chunk = self._chunk(kw, i)
                for line in chunk:
                    if identif in line: foundit = True
                    break
                if foundit: break
            if foundit:
                return self._raw_data(chunk)
            else:
                print('multiple data exists - did not find identifier "%s"' % identif)
                return []
        # - one unique chunk. most often this is the case
        chunk = self._chunk(kw, 1)
        return self._raw_data(self._chunk(kw, 1))
#
    def raw2values(self, raw, col1, col2, skip=0, identif='', match_col=0, matrix=False):
        '''
        extracts numeric values from the data provided by get_raw_data.
        sometimes you want to process the raw data yourself, but often this
        method will save you some programming.
        raw       : output from get_raw_data
        col1      : start-column
        col2      : end-column
        skip      : number of headerlines to skip. typically 0 or 1 (could be higher)
        identif   : some string to look for, typically well name. note: '' matches anything
        match_col : which column to search for 'identif'
        matrix    : force result to be matrix
        note1: column-indexing uses the natural numbering (not python-index)
        note2: gives you a NaN for defaulted values
        note3: returns array if col1==col2 or if only one line is found, else it gives a matrix.
               this can be overridden using the 'matrix' parameter
        '''
        ncols = col2 - col1 + 1
        m = []  # to be returned
        for rec in raw[skip:]:
            if not identif in rec[match_col-1]: continue
            row = pl.NaN*pl.ones(ncols)
            for j in range(col1-1,col2):
                indx = j-col1+1  # starting at 0
                if len(rec) <= j: row[indx] = pl.NaN # values 'far right' are defaulted
                else:
                    val = rec[j]
                    if val == '1*': row[indx] = pl.NaN
                    else          : row[indx] = float(val)
            m.append(row)
        if matrix:       return pl.array(m) # matrix
        if len(m) == 1 : return m[0]        # array
        m = pl.array(m)
        if col1 == col2: return m[:,0]      # array
        else           : return m           # matrix
#
    def get_values(self, kw, col1, col2, skip=0, identif='', matrix=False):
        '''
        an easy way of getting data (often sufficient) - more simple than dd.raw2values(dd.get_raw_data(...)),
        but not as versatile
        note1: EQUALS is handled specially since it is often used. must have identif in this case (e.g. 'PORO' or 'PERMX')
        note2: if it fails to give you want you need, try using get_raw_data instead - maybe with get_all=True
        '''
        if kw == 'EQUALS':
            ll = self.get_raw_data('EQUALS', get_all=True)
            for l in ll[0]:
                if l[0] == identif:
                    return [float(x) for x in l[col1:col2+1]]
            print('Did not find identifier %s in EQUALS' % identif)
            return []
        return self.raw2values(self.get_raw_data(kw, identif), col1, col2, skip=skip, matrix=matrix)
#
    def get_dates(self):
        '''
        useful for getting dates from SCHEDULE section.
        for now, only works for DATES keyword
        '''
        if not self.has_kw('DATES'): return []
        dates = [pl.datetime.datetime(int(x[2]), MONTH_MAP[x[1]], int(x[0])) for x in self.get_raw_data('DATES', get_all=1)]
        return pl.date2num(dates)
