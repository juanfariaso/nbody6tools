from nbody6tools.Reader import read_snapshot,get_number_of_snapshots
from nbody6tools.ext import qparameter
from nbody6tools import Utilities

import numpy

local_variables = locals()

def compute(folder,function,args="",output=None,overwrite=False,fmt = "%10.7g",doc=False,**kw):
    """ 
    Evaluate function at every snapshot and saves result
    in a file
    folder                 : Simulation folder
    funcion(snapshot,**kw) : the function to evaluate, should accept
                        a snapshot object as input and return
                        a tuple of the resulting values (if more than one)
    args                   : extra arguments (list of strings like: ['arg1=value']
                             'value' can only be floats.
    output                 : output file with results. Default 'function_name'.dat is given.
    overwrite              : if False raise Error if output file exists. Otherwise overwrite. Default: True
    fmt                    : format for numbers in columns. Default: '%10.7g' 
    """

    function = local_variables[function]

    if doc:
        print(function.__doc__)
        return

    #parse the extra arguments if any
    kwargs = dict()
    for s in args:
        k,v=s.split("=")
        try :
            kwargs[k] = float(v)
        except ValueError:
            raise("Can only use float arguments for now")

    ns = get_number_of_snapshots(folder)
    if output is None:
        output =  function.__name__+".dat"

    mode = "w" if overwrite else "wx"
    resultfile =  open(output,mode,1)   #fail if already exists for safety

    resultfile.write("# function %s with extra arguments: %s \n"%(
        function.__name__,str(args) ) )
    resultfile.write("# time   time[Myr]  results  \n")
    for i in range(ns):
        print("Snapshot: %i/%i    "%(i,ns))
        sn = read_snapshot(folder,i,inputfilename="input")
        t = sn.parameters["time"]
        tscale = sn.parameters["tscale"]
        resultfile.write( (fmt+"  "+fmt)%( t,t*tscale ) )
        fout = function(sn,**kwargs)
        print(str(fout)+"\r")
        if hasattr(fout,"__len__"):
           nout = len(fout) if hasattr(fout,"__len__") else 1
           for j in range(len(fout)):
               resultfile.write( " "+fmt%(fout[j] ) )
        resultfile.write("\n")
    resultfile.close()

def Qpar(snapshot,average=1,zeroaxis=1,rmax=0.9,**args):
    """ 
    Calculate the Q parameter (##add ref).
    extra arguments:
    average : if 1 average the three projections. Default: 1
    zeroaxis: if average=0, will set the axis which position will be taken as zero. keys are: 1-> x ; 2-> y, 3-> z
    rmax    : Mass fraction radius to use as maximum. Default : 0.9
    
    returns "mean Qparameter","standard error of mean"
    """
    #snapshot.to_physical() #unnesessary
    rmax = Utilities.get_mass_radius(snapshot.stars,rmax)
    x = snapshot.stars["x"]
    y = snapshot.stars["y"]
    z = snapshot.stars["z"]
    r = numpy.sqrt(x**2+y**2+z**2)
    mask =  r < rmax
    return qparameter(x[mask],y[mask],z[mask],int(average),int(zeroaxis),rmax)