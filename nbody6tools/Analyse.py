from nbody6tools.Reader import read_snapshot,get_number_of_snapshots
from nbody6tools.ext import qparameter
from nbody6tools import Utilities
from nbody6tools.Datamodel import get_binaries_from_files
import sys

import numpy

local_variables = locals()

def compute(folders,function,args=None,output=None,overwrite=False,
            fmt ="%10.7g",doc=False,stdout=None,**kw):
    """ 
    Evaluate function at every snapshot and saves result
    in a file
    folder                 : Simulation folder
    funcion(snapshot,**kw) : the function to evaluate, should accept
                        a snapshot object as input and return
                        a tuple of the resulting values (if more than one)
    args                   : extra arguments (list of strings like: ['arg1=value']
                             'value' can only floats. comma sepparated arguments will be used as lists
                             e.g. 'arg1=1.0,2,3.0,4' is interpreted as 'arg1 = [1.0,2.0,3.0,4.0]'
    output                 : output file with results. Default 'function_name'.dat is given.
    overwrite              : if False raise Error if output file exists. Otherwise overwrite. Default: True
    fmt                    : format for numbers in columns. Default: '%10.7g' 
    """

    function = local_variables[function]

    if doc:
        print(function.__doc__)
        return

    scriptmode = False
    if stdout is not None:
        sys.stdout = open(stdout,"a")
        sys.stderr = open(stdout+"err","a")
        scriptmode=True

    kwargs = dict()
    if args is not None:
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

    for folder in folders:
        #parse the extra arguments if any
        if len(folders) > 1 and not scriptmode:
            print("Working on folder %s"%folder)

        ns = get_number_of_snapshots(folder)
        if output is None:
            output =  function.__name__+".dat"
        outputfile = "%s/%s"%(folder,output)

        mode = "w" if overwrite else "x"
        try:
            resultfile =  open(outputfile,mode,1)   #fail if already exists for safety
        except FileExistsError:
            print("File %s exists, skipping. Try --overwrite option"%outputfile)
            sys.stdout.flush()
            sys.stderr.flush()
            continue

        resultfile.write("# function %s with extra arguments: %s \n"%(
            function.__name__,str(args) ) )
        resultfile.write("# time   time[Myr]  results  \n")
        for i in range(ns):
            if not scriptmode:
                print("Snapshot: %i/%i    "%(i,ns),end="\r")
            try:
                sn = read_snapshot(folder,i,inputfilename="input")
            except:
                break
            t = sn.parameters["time"]
            tscale = sn.parameters["tscale"]
            resultfile.write( (fmt+"  "+fmt)%( t,t*tscale ) )
            fout = function(sn,**kwargs)
            #print(str(fout)+"\r")
            nout = len(fout) if hasattr(fout,"__len__") else 1
            if nout > 1:
                for j in range(nout):
                    resultfile.write( " "+fmt%(fout[j] ) )
            else:
                resultfile.write( " "+fmt%(fout) )
            resultfile.write("\n")
        if not scriptmode:
            print("\n")
        else:
            print("%s Done"%outputfile)
            sys.stdout.flush()
            sys.stderr.flush()
        resultfile.close()

def Qpar(snapshot,average=1,zeroaxis=1,rmax=0.9,**args):
    """ 
    Calculate the Q parameter (##add ref).
    extra arguments:
    average : if 1 average the three projections. Default: 1
    zeroaxis: if average=0, will set the axis which position will be taken as zero. keys are: 1-> x ; 2-> y, 3-> z
    rmax    : Mass fraction radius to use as maximum. Default : 0.9
    
    returns "mean Qparameter","standard error of mean", current rmax
    """
    #snapshot.to_physical() #unnesessary
    rmax = Utilities.get_mass_radius(snapshot.stars,rmax)
    x = snapshot.stars["x"]
    y = snapshot.stars["y"]
    z = snapshot.stars["z"]
    r = numpy.sqrt(x**2+y**2+z**2)
    mask =  r < rmax
    result = qparameter(x[mask],y[mask],z[mask],int(average),int(zeroaxis),rmax)
    return result[0],result[1],rmax

def lrad(snapshot,mfrac = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1.0],**args) :
    """
    Custom Lagrangian radii calculation.
    extra arguments:
    mfrac : list of lagrangian radii to calculate
            default: 0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1.0

    """
    result = []
    for fr in mfrac:
        result.append(Utilities.get_mass_radius(snapshot.stars,fraction=fr) )
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
    snapshot.to_physical()
    snapshot.to_center()
    #if snapshot.parameters["kz"][8]>0 :
    stars = snapshot.unresolve_all()
    #else:
    #    stars = snapshot.stars
    
    bound_set = stars.bound_subset()
    sigma = bound_set.velocity_dispersion()
    rh = bound_set.half_mass_radius()
    return bound_set.mass.sum() / stars.mass.sum(), sigma, rh

def virial_ratio(snapshot,bound = True):
    """ Calculate the virial ratio of snapshot:

    bound : if True only consider bound stars
    smoothing_length : Softhening length for avoid potential get too high.
    """
    #snapshot.to_physical()
    snapshot.to_center()
    #stars = snapshot.stars
    stars = snapshot.unresolve_all()

    cmv = stars.center_of_mass_velocity()
    stars.to_center_velocity(cmv)

    #if bound:
    #    mask = stars.pot + stars.mass*numpy.sqrt(stars.vx**2 + stars.vy**2 +
                                                 #stars.vz**2)*0.5 < 0
        #stars = stars.bound_subset()
    #    stars = stars[mask]
    pot= (stars.pot*stars.mass).sum()*0.5 + (stars.epot*stars.mass).sum()
    return -stars.kinetic_energy()/pot
    #return stars.virial_ratio()

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

def velocity_dispersion(snapshot) :
    snapshot.to_physical()
    stars = snapshot.unresolve_all()
    return stars.velocity_dispersion()

def bound_statistics(snapshot):
    snapshot.to_physical()
    stars = snapshot.unresolve_all()
    boundset = stars.bound_subset()

    fbound = boundset.mass.sum()/stars.mass.sum()
    rh = boundset.half_mass_radius()
    vdisp  =  boundset.velocity_dispersion()
    return fbound,rh,vdisp

def number_density(snapshot,mass_fraction_radii = 0.5) :
    """
    Compute the number density on a radius given by mass_fraction_radius.
    Default : 0.5, i.e. Half mass radius.
    returns : number density at given radii in stars / pc**3
    """
    snapshot.to_physical()
    stars = snapshot.unresolve_all()
    r = stars.mass_radius(fraction=mass_fraction_radii) 
    return len(stars) * 3./4./numpy.pi/r**3 
