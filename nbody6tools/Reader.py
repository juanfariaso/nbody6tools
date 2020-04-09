import nbody6tools
import subprocess
from nbody6tools import Datamodel
import numpy
from scipy.io import FortranFile
import os
import glob

try:
    from customReader import Snapshot
    print("Custom reader imported OK")
except (ImportError):
    from nbody6tools.Datamodel import Snapshot

Options = nbody6tools.options
singlefile = Options.getboolean("SNAPSHOTS","singlefile")

def parse_inputfile(inputfilename,**kw):
    """%s
    """%Datamodel.parse_inputfile.__doc__
    return Datamodel.parse_inputfile(inputfilename,**kw)

def get_number_of_snapshots(folder):
    l=[x.replace("%sconf.3_"%folder,"") for x in glob.glob(folder+"conf.3*") ]
    l.sort(key=float)
    if len(l) == 0:
        raise ValueError("No snapshots in this folder.")
    return len(l)

def get_times(folder,dtype=float):
    l=[x.replace("%sconf.3_"%folder,"") for x in glob.glob(folder+"conf.3*") ]
    l.sort(key=float)
    l = numpy.array(l,dtype=dtype)
    return l

def read_snapshot(folder,snapshot=0,inputfilename="input",singlefile=singlefile,
        snapshotfile = "conf.3"):

    inputfile ="%s/%s"%(folder,inputfilename)

    if not singlefile:
        if folder[-1] != "/" : folder +="/"
        opt = Datamodel.parse_inputfile(folder+"/"+inputfilename)
        kz = opt["KZ"]
       
        l=[x.replace("%sconf.3_"%folder,"") for x in glob.glob(folder+"%s*"%snapshotfile ) ]
        l.sort(key=float)
        if len(l) == 0:
            raise ValueError("No snapshots in this folder.")

        snapshotfile = "conf.3_%s"%l[snapshot]
    else:
        if not os.path.isfile( "%s/%s"%(folder,snapshotfile)):
            raise FileNotFoundError("Snapshot file %s/%s not found"%(folder,snapshotfile))

    return Snapshot("%s/%s"%(folder,snapshotfile),inputfile,singlefile,snapshot=snapshot)

def read_binaries(folder,snapshot=0,inputfilename="input"):
    """
    Returns two dictionaries with binary properties
    regularized binaries , wide binaries
    """
    #TODO raise error when needed kz option is not set
    times = get_times(folder,dtype=str)
    widefile = "%s/bwdat.19_%s"%(folder,times[snapshot])
    hardfile = "%s/bdat.9_%s"%(folder,times[snapshot] )
    
    return Datamodel.get_binaries_from_files(hardfile,widefile)
