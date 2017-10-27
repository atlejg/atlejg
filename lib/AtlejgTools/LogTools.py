LOG_LEVEL_DEBUG = 999
LOG_LEVEL_INFO  = 3
LOG_LEVEL_WARN  = 2
LOG_LEVEL_ERROR = 1
LOG_LEVEL_NONE  = 0

class Log(object):

   def __init__(self,level=LOG_LEVEL_INFO):
      self.level = level

   def _output(self,prefix,level,msg):
      if level > self.level : return
      print prefix + " : ",msg

   def debug(self,msg):
      self._output('debug',LOG_LEVEL_DEBUG,msg)

   def info(self,msg):
      self._output('info',LOG_LEVEL_INFO,msg)

   def warn(self,msg):
      self._output('warn',LOG_LEVEL_WARN,msg)

   def error(self,msg):
      self._output('error',LOG_LEVEL_ERROR,msg)

   def write(self,msg,level):
      self._output('write',level,msg)

   def set_level(self,level):
      self.level = level

   def get_level(self):
      return self.level
