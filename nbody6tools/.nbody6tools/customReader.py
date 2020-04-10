## This is a configuration file for nbody6tools
## Any nesessary customization in order to read snapshots
## from Nbody6 or a custom versio should be here

## A copy of the default Snapshot object is here for reference

from nbody6tools.Datamodel import Snapshot
from scipy.io import FortranFile
import numpy

#check Datamodel.__init__.py for the original file
class Snapshot(Snapshot):
    def __read_record(self):
      NTOT,MODEL,NDRUN,NK = tuple(self.__recordfile.read_ints(dtype=numpy.int32))
      self._ntot = NTOT
      record = self.__recordfile.read_record(
              numpy.dtype( (numpy.float32,NK)   ) #AS
             ,numpy.dtype( (numpy.float32,NTOT) ) #bodys
             ,numpy.dtype( (numpy.float32,NTOT) ) #rhos
             ,numpy.dtype( (numpy.float32,NTOT) ) #xns
             ,numpy.dtype( (numpy.float32,(NTOT*3) ) )#x
             ,numpy.dtype( (numpy.float32,(NTOT*3) ) )#v
             ,numpy.dtype( (numpy.float32,NTOT  ) )#phi
             ,numpy.dtype( (numpy.int32  ,NTOT) ) #name
             ,numpy.dtype( (numpy.float32,NTOT) ) #time
           )
      return record
