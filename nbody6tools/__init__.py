from os.path import expanduser
import sys
try:
    from configparser import ConfigParser
except ImportError :
    from ConfigParser import ConfigParser
sys.path.append( expanduser("~/.nbody6tools") )
options = ConfigParser()
options.read( expanduser("~/.nbody6tools/nb6.conf" )  )
