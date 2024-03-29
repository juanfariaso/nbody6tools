from nbody6tools.Reader import read_snapshot,get_number_of_snapshots
#from nbody6tools.ext import qparameter
from nbody6tools import Utilities
from nbody6tools.Datamodel import get_binaries_from_files
import sys

import numpy
from multiprocessing import Pool

try :
    from tqdm import tqdm
    fancy_bar = True
except ImportError:
    fancy_bar = False
    

local_variables = locals()

def compute(folders,function,args=None,output=None,overwrite=False,
            fmt ="%10.7g",doc=False,stdout=None,badTol=10,bar_pos=None,
            leave_bar=True,strict =False,**kw):
    """ 
    Evaluate function at every snapshot and saves result
    in a file
    folder                 : Simulation folder
    funcion(snapshot,**kw) : the function to evaluate, should accept
                        a snapshot object as input and return
                        a tuple of the resulting values (if more than one).
                        If function has a docstring, a line will be searched
                        with the format:
            headers: header1  header2 header3 ..
                        everithing after "headers:" will be used as header for the
                        outupt file. header1 should identify the first output
                        value of function, etc.
    args                   : extra arguments (list of strings like: ['arg1=value']
                             'value' can only floats. comma sepparated arguments will be used as lists
                             e.g. 'arg1=1.0,2,3.0,4' is interpreted as 'arg1 = [1.0,2.0,3.0,4.0]'
    output                 : output file with results. Default 'function_name'.dat is given.
    overwrite              : if False raise Error if output file exists. Otherwise overwrite. Default: True
    fmt                    : format for numbers in columns. Default: '%10.7g' 
    badTol                 : Number snapshots that fails to evalueate function
                             that are ignored. Above this calculation is
                             aborted for this folder.
    """

    #check if functuion is a customized function, otherwise look in this module
    if not callable(function):
        function = local_variables[function]


    if doc:
        print(function.__doc__)
        return

    headers = ""
    if function.__doc__ is not None:
        for line in function.__doc__.split("\n"):
            if "headers:" in line: 
                headers = line.split(":")[-1].replace("\n","")
        
    

    if strict:
        iException = AssertionError
    else:
        iException = Exception

    scriptmode = False
    if stdout is not None:
        sys.stdout = open(stdout,"a")
        sys.stderr = open(stdout+"err","a")
        scriptmode=True

    kwargs = dict()
    if args is not None:
        if type(args) == str:
            for s in args:
                k,v=s.split("=")
                try :
                    if "," in v:
                        #list case
                        kwargs[k] = [float(val) for val in v.split(",")]
                    else:
                        #single float case
                        kwargs[k] = float(v)
                except ValueError:
                    raise ValueError("arg can only be a float or "
                                     "comma separated floats")
        else :
            kwargs = args

    for folder in folders:
        nbad = 0 
        #parse the extra arguments if any
        if len(folders) > 1 and not scriptmode and not type(folders)==str:
            print("Working on folder %s"%folder)

        ns = get_number_of_snapshots(folder)
        if output is None:
            outputfile  =  folder +"/"+ function.__name__+".dat"
        else:
            outputfile = output

        #print(outputfile)
        mode = "w" if overwrite else "x"
        try:
            resultfile =  open(outputfile,mode,1)   #fail if already exists for safety
        except FileExistsError:
            print("File %s exists, skipping. Try --overwrite option"%outputfile)
            sys.stdout.flush()
            sys.stderr.flush()
            if len(folders) == 1:
                raise FileExistsError()
            continue

        resultfile.write("# function %s with extra arguments: %s \n"%(
            function.__name__,str(args) ) )
        resultfile.write("# time   time[Myr]  %s  \n"%headers )
        if not scriptmode and fancy_bar:
            pbar = tqdm(total=ns, position=bar_pos,
                    desc="[%s]%s"%(bar_pos,folder),leave=leave_bar)
        nresult = 1
        for i in range(ns):
            if not scriptmode:
                if fancy_bar :
                    pbar.update(1)
                else:
                    print("Snapshot: %i/%i    "%(i,ns),end="\r")
            try:
                sn = read_snapshot(folder,i,inputfilename="input")
            except :
                break
            t = sn.parameters["time"]
            tscale = sn.parameters["tscale"]
            resultfile.write( (fmt+"  "+fmt)%( t,t*tscale ) )
            try:
                fout = function(sn,**kwargs)
            except iException:
                nbad+=1
                if nbad > badTol :
                    print("Warning: Too many bad evaluations of function %s"%
                          function.__name__)
                    print(" skipping folder %s"%folder)
                    break
                else:
                    for _ in range(nresult) :
                        resultfile.write("  nan ")
                    resultfile.write("\n")
                    continue
            if hasattr(fout,"__len__") :
                nresult = len(fout)
            #print(str(fout)+"\r")
            nout = len(fout) if hasattr(fout,"__len__") else 1
            if nout > 1:
                for j in range(nout):
                    resultfile.write( " "+fmt%(fout[j] ) )
            else:
                resultfile.write( " "+fmt%(fout) )
            resultfile.write("\n")
        if not scriptmode:
            if fancy_bar:
                pbar.close()
            else:
                print("\n")
        else:
            print("%s Done"%outputfile)
            sys.stdout.flush()
            sys.stderr.flush()
        resultfile.close()
    return 0

#def Qpar(snapshot,average=1,zeroaxis=1,rmax=0.7,**args):
#    """ 
#    Calculate the Q parameter (##add ref).
#    extra arguments:
#    average : if 1 average the three projections. Default: 1
#    zeroaxis: if average=0, will set the axis which position will be taken as zero. keys are: 1-> x ; 2-> y, 3-> z
#    rmax    : Mass fraction radius to use as maximum. Default : 0.9
#    
#    returns "mean Qparameter","standard error of mean", current rmax
#    """
#    #snapshot.to_physical() #unnesessary
#    stars = snapshot.unresolved_stars
#    rmax = Utilities.get_mass_radius(stars,rmax)
#    x = stars["x"]
#    y = stars["y"]
#    z = stars["z"]
#    r = numpy.sqrt(x**2+y**2+z**2)
#    mask =  r < rmax
#    result = qparameter(x[mask],y[mask],z[mask],int(average),int(zeroaxis),rmax)
#    return result[0],result[1],rmax

def lrad(snapshot,mfrac = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.0],**args) :
    """
    Custom Lagrangian radii calculation.
    extra arguments:
    mfrac : list of lagrangian radii to calculate
            default: 0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1.0

    """
    snapshot.to_physical()
    result = []
    for fr in mfrac:
        result.append(Utilities.get_mass_radius(snapshot.stars,fraction=fr) )
    return result

def lradb(snapshot,mfrac = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.0],**args) :
    """
    Custom Lagrangian radii calculation for bound set.
    extra arguments:
    mfrac : list of lagrangian radii to calculate
            default: 0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1.0

    """
    snapshot.to_physical()
    result = []
    for fr in mfrac:
        stars = snapshot.bound_stars_unresolved
        result.append(Utilities.get_mass_radius(stars,fraction=fr) )
    return result

def mr_track(snapshot,mass_fraction=0.5,**args):
    snapshot.to_physical()
    radius = Utilities.get_mass_radius(snapshot.stars,fraction=mass_fraction) 
    mass = snapshot.stars.mass.sum()
    return radius,mass

def bound_fraction(snapshot,**args):
    """
    Compute bound mass fraction and basic statistics of bound stars
      return:
      bound_fraction, sigma_bound ,rh_bound
    """
    bound_set = snapshot.bound_stars_unresolved
    sigma = bound_set.velocity_dispersion()
    rh = bound_set.half_mass_radius()
    return bound_set.mass.sum() / snapshot.stars.mass.sum(), sigma,rh ,len(bound_set)

def virial_ratio(snapshot,bound = False):
    """ Calculate the virial ratio of snapshot:

    bound : if True only consider bound stars
    smoothing_length : Softhening length for avoid potential get too high.
    """
    if bound :
        stars = snapshot.bound_stars_unresolved
    else:
        stars = snapshot.unresolved_stars
    return stars.kinetic_energy() / snapshot.virial_energy_of_set(stars)

def binary_fraction(snapshot):
    widefile = snapshot._snapshotfile.replace("conf.3_","bwdat.19_")  
    hardfile = snapshot._snapshotfile.replace("conf.3_","bdat.9_") 
    hard,wide = get_binaries_from_files(hardfile,widefile)
    nb = len(hard["primary"]) + len(wide["primary"] ) 
    return float(nb)/ (snapshot.n - nb ) 

def binary_energy(snapshot):
    widefile = snapshot._snapshotfile.replace("conf.3_","bwdat.19_")  
    hardfile = snapshot._snapshotfile.replace("conf.3_","bdat.9_") 
    hard,wide = get_binaries_from_files(hardfile,widefile)
    ebin_hard = hard["ebin"].sum()
    ebin_wide = hard["ebin"].sum()

    return ebin_hard,ebin_wide,ebin_hard+ebin_wide

def velocity_dispersion(snapshot,bound=False) :
    snapshot.to_physical()
    if bound:
        stars = snapshot.bound_stars_unresolved
    else:
        stars = snapshot.unresolved_stars
    return stars.velocity_dispersion()

def number_density(snapshot,mass_fraction_radii = 0.5,bound=False) :
    """
    Compute the number density on a radius given by mass_fraction_radius.
    Default : 0.5, i.e. Half mass radius.
    returns : number density at given radii in stars / pc**3
    """
    snapshot.to_physical()
    #print(bound, "b is false?",bound == False)
    #if bound :
    #    stars = snapshot.bound_stars_unresolved #unresolved is faster
    #else:
    #    stars = snapshot.stars[snapshot.stars.r < 1e6]
    stars = snapshot.stars
    r = stars.mass_radius(fraction=mass_fraction_radii) 
    stars.to_center()
    stars = stars[ stars.r < r ]

    n = (stars.name <= snapshot.n).sum() + 2*(stars.name > snapshot.n ).sum()
    return n * 3./4./numpy.pi/r**3 

def core_radius(snapshot):
    rcore = snapshot.parameters["rc"] * snapshot.parameters["rbar"]
    snapshot.to_physical()
    stars = snapshot.unresolved_stars
    #rcenter = stars.mass_radius(fraction=0.01)
    #r = numpy.sqrt(stars.x**2 + stars.y**2 + stars.z**2)
    r = stars.r
    nc = (r<=rcore).sum()
    return  rcore , 3.*nc/4./numpy.pi/rcore**3
