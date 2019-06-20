import subprocess
from nbody6tools import Datamodel
import numpy
from scipy.io import FortranFile
import os
import glob

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

def get_times(folder):
    l=[x.replace("%sconf.3_"%folder,"") for x in glob.glob(folder+"conf.3*") ]
    l.sort(key=float)
    l = numpy.array(l,dtype=float)
    return l

def read_snapshot(folder,snapshot=0,inputfilename="input"):
    if folder[-1] != "/" : folder +="/"
    opt = Datamodel.parse_inputfile(folder+"/"+inputfilename)
    kz = opt["KZ"]
   
    l=[x.replace("%sconf.3_"%folder,"") for x in glob.glob(folder+"conf.3*") ]
    l.sort(key=float)
    if len(l) == 0:
        raise ValueError("No snapshots in this folder.")

    snapshotfile ="%s/conf.3_%s"%(folder,l[snapshot])
    inputfile ="%s/%s"%(folder,inputfilename)
    return Datamodel.Snapshot(snapshotfile,inputfile)
