from os.path import expanduser
import sys
from configparser import ConfigParser
sys.path.append( expanduser("~/.nbody6tools") )
options = ConfigParser()
options.read( expanduser("~/.nbody6tools/nb6.conf" )  )
