import nbody6tools
from nbody6tools import Datamodel
from nbody6tools.Datamodel import Interpolators
import numpy
import os
import glob

try:
    from customReader import Snapshot
    print("Custom reader imported OK")
except (ImportError):
    from nbody6tools.Datamodel import Snapshot

Options = nbody6tools.options
singleFile = Options.getboolean("CONFIG","singlefile",fallback=False)
snapshotFile = Options.get("CONFIG","snapshotFile",fallback="conf.3_" )
inputFile = Options.get("CONFIG","inputFile",fallback="input" )


def get_globals(folder):
    "find last block of globals that is good"
    x = open("%s/global.30"%folder,"r").readlines()
    it = [ i for i in range(len(x)) if "TIME" in x[i] ]  #find headers if more than one
    #header = x[it[-1]]
    header = "TIME[NB] TIME[Myr] TCR[Myr] DE BE(3) RSCALE[PC] RTIDE[PC] RDENS[PC] RC[PC]  RHOD[M*] RHOM[M*] MC[M*] CMAX <Cn> Ir/R RCM VCM AZ EB/E EM/E VRMS N NS NPAIRS NUPKS NPKS NMERGE MULT <NB> NC NESC NSTEPI NSTEPB NSTEPR NSTEPU NSTEPT NSTEPQ NSTEPC NBLOCK NBLCKR NNPRED NIRRF NBCORR NBFLUX NBFULL NBVOID NICONV NLSMIN NBSMIN NBDIS NBDIS2 NCMDER NFAST NBFAST NKSTRY NKSREG NKSHYP NKSPER NKSMOD NTTRY NTRIP NQUAD NCHAIN NMERG NEWHI"
    header = header.split()
    data = numpy.loadtxt(x[it[-1]+1:] )
    result = dict()
    for i,h in enumerate(header):
        result[h] = data[:,i]
    return result

def get_events(folder):
    "find last block of globals that is good"
    x = open("%s/event.35"%folder,"r").readlines()
    it = [ i for i in range(len(x)) if "TIME" in x[i] ]  #find headers if more than one
    #header = x[it[-1]]
    header = "TIME[Myr] NDISS  NTIDE  NSYNC NCOLL NCOAL NDD NCIRC NROCHE NRO NCE NHYP NHYPC NKICK EBIN EMERGE ECOLL EMDOT ECDOT EKICK ESESC EBESC EMESC DEGRAV EBIND MMAX NMDOT NRG NHE NRS NNH NWD NSN NBH NBS ZMRG ZMHE ZMRS ZMNH ZMWD ZMSN ZMDOT"
    header = header.split()
    data = numpy.loadtxt(x[it[-1]+1:] )
    result = dict()
    for i,h in enumerate(header):
        result[h] = data[:,i]
    result["NTYPE"]= [0]
    for i in range(42,57):
        result["NTYPE"].append( data[:,i] )
    return result

def parse_inputfile(inputfilename=inputFile,**kw):
    """%s
    """%Datamodel.parse_inputfile.__doc__
    return Datamodel.parse_inputfile(inputfilename,**kw)

def get_number_of_snapshots(folder,inputfilename=inputFile,**kw):
    return len(get_times(folder,inputfilename))

def get_times(folder,nbody=False,inputfilename=inputFile,
        snapshotfile = snapshotFile,singlefile=singleFile):
    glfile = get_globals(folder)
    inputfile = parse_inputfile("%s/%s"%(folder,inputfilename) )
    key = "TIME[NB]" if nbody else "TIME[Myr]"

    gltimes = glfile[key][0::inputfile["NFIX"] ]

    if not singlefile:
        #sometimes the last snapshot is not saved. Check that it exist.
        lsnaps=[x.split("_")[-1] for x in glob.glob(folder+"%s*"%snapshotfile ) ]

    return gltimes[:len(lsnaps)]

def read_snapshot(folder,snapshot=0,time=None,inputfilename=inputFile,singlefile=singleFile,
        snapshotfile = snapshotFile):
    if time is not None:
        times = numpy.array(get_times(folder))
        snapshot = numpy.argmin( abs( times - time) )
        #snapshot = numpy.where( times >= time )[0][0]

    inputfile ="%s/%s"%(folder,inputfilename)

    if not singlefile:
        if folder[-1] != "/" : folder +="/"
        #opt = Datamodel.parse_inputfile(folder+"/"+inputfilename)
        #kz = opt["KZ"]

        #l=[x.replace("%s%s"%(folder,snapshotfile),"") for x in glob.glob(folder+"%s*"%snapshotfile ) ]
        l=[x.split("_")[-1] for x in glob.glob(folder+"%s*"%snapshotfile ) ]
        l.sort(key=float)
        if len(l) == 0:
            raise ValueError("No snapshots in this folder.")

        snapshotfile = "%s%s"%(snapshotfile,l[snapshot])
    else:
        if not os.path.isfile( "%s/%s"%(folder,snapshotfile)):
            raise FileNotFoundError("Snapshot file %s/%s not found"%(folder,snapshotfile))

    return Snapshot("%s/%s"%(folder,snapshotfile),inputfile,singlefile,snapshot=snapshot)

def read_binaries(folder,snapshot=0,inputfilename=inputFile):
    """
    Returns two dictionaries with binary properties
    regularized binaries , wide binaries
    """
    #TODO raise error when needed kz option is not set and check how to deal
    # with single snapshotFilecase
    l=[x.replace("%s%s"%(folder,snapshotFile),"")
            for x in glob.glob(folder+"%s*"%snapshotFile) ]
    l.sort(key=float)
    times = numpy.array(l,dtype=str)

    #times = get_times(folder,nbody=True)
    widefile = "%s/bwdat.19_%s"%(folder,times[snapshot])
    hardfile = "%s/bdat.9_%s"%(folder,times[snapshot] )

    return Datamodel.get_binaries_from_files(hardfile,widefile)

def get_orbit_interpolator(folder,outputfile = None):
    return Interpolators.ClusterOrbitInterpolator(folder = folder, 
                                              outputfile = outputfile)
