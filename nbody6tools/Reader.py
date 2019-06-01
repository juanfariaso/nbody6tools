import subprocess
import numpy
from scipy.io import FortranFile
import os
import glob

class Snapshot(object):
    """
    An object containing the information of a snapshot of Nbody6.
    Contains :
      self.stars           : Containing the stars data. Type = dict.
      self.parameters      : Conatining several from AS variable. Type=dict.
      self.inputfile       : Inputfile of the simulation. Type = dict
      self.nsingles        : Number of singles
      self.npairs          : Number of regularized binaries
      self.physical        : Units of the star set True : Msun,parsec,km/s ; False nbody-units: 

      Methods: 
      self.to_physcal()     : Transform stars to physical units (Msun,pc,kms)
      self.to_nbody()       : Transform stars to nbody units

      self.to_center(center): Move stars to the specified center. If not specified, use the center of density calculated by Nbody6
      self.reorder(isort)   : Sort stars according to the 'isort' list (e.g. result of numpy.argsort(). If nothing specified, sort by name.
    """

    def __init__(self,snapshotfile,inputfile) :
        self._inputfile = parse_inputfile(inputfile)

        record = self._read_snapshot(snapshotfile)
        self._stars = dict()
        self._parameters = dict()
        self._structure(record)

    @property
    def parameters(self):
        """ Dict. containing parameters read from the AS variable.
        """
        return self._parameters

    @property
    def stars(self):
        """ Dict. containing stellar coordinates:
        mass, x,y,z, vx,vy,vz
        """
        return self._stars

    @property
    def inputfile(self):
        """ Dict. containing input parameters from the inputfile
        """
        return self._inputfile

    @property
    def nsingles(self):
        """ Number of single particles
        """
        return self._nsingles

    @property
    def npairs(self):
        """ Number of regularized binaries
        """
        return self._npairs

    @property
    def ntot(self):
        """  Total number of particles including center of mass particles
        """
        return self._ntot

    @property
    def n(self):
        """ Number of particles in the simulation
        """
        return self._n

    @property
    def time(self):
        """ Time of snapshot
        """
        return self._time

    @property
    def physical(self):
        """ bool, True : stars are in phycial units: Msun, parsec, km/s
              False: stars are in code units
        """
        return self._physical


    def _read_snapshot(self,name):
      kz = self.inputfile["KZ"]
      f = FortranFile(name,"r")
      datafile = open(name,"r")
      NTOT,MODEL,NDRUN,NK = tuple(f.read_ints(dtype=numpy.int32))
      self._ntot = NTOT
      if not kz[50] == 1 : 
         record = f.read_record(
                 numpy.dtype( (numpy.float32,NK)   ) #AS
                ,numpy.dtype( (numpy.float32,NTOT) ) #bodys
                ,numpy.dtype( (numpy.float32,NTOT) ) #rhos
                ,numpy.dtype( (numpy.float32,NTOT) ) #xns
                ,numpy.dtype( (numpy.float32,(NTOT*3) ) )#x
                ,numpy.dtype( (numpy.float32,(NTOT*3) ) )#v
                ,numpy.dtype( (numpy.float32,NTOT  ) )#phi
                ,numpy.dtype( (numpy.int32  ,NTOT) ) #name
              )
      if kz[50] == 1 : 
         record = f.read_record(
                 numpy.dtype( (numpy.float32,NK)   ) #AS
                ,numpy.dtype( (numpy.float32,NTOT) ) #bodys
                ,numpy.dtype( (numpy.float32,NTOT) ) #rhos
                ,numpy.dtype( (numpy.float32,NTOT) ) #xns
                ,numpy.dtype( (numpy.float32,(NTOT*3) ) )#x
                ,numpy.dtype( (numpy.float32,(NTOT*3) ) )#v
                ,numpy.dtype( (numpy.float32,NTOT  ) )#phi
                ,numpy.dtype( (numpy.int32  ,NTOT) ) #name
                ,numpy.dtype( (numpy.int32  ,NTOT) ) #name
              )
      return record

    def _structure(self,record):
        self._parameters["time"] = record[0][0]
        self._parameters["npairs"] = int(record[0][1])
        self._parameters["rbar"] = record[0][2]
        self._parameters["zmbar"] = record[0][3]
        self._parameters["rtide"] = record[0][4]
        self._parameters["tidal4"] = record[0][5]
        self._parameters["rdens"] = numpy.array(record[0][6:9])
        self._parameters["ntcr"] = record[0][9] # number of elapsed initial crossing times 
        self._parameters["tscale"] = record[0][10]
        self._parameters["vstar"] = record[0][11]
        self._parameters["rc"] = record[0][12]
        self._parameters["nc"] = record[0][13]
        self._parameters["vc"] = record[0][14]
        self._parameters["rhom"] = record[0][15]
        self._parameters["cmax"] = record[0][16]
        self._parameters["rscale"] = record[0][17] #half mass radius
        self._parameters["rsmin"] = record[0][18]
        self._parameters["dmin1"] = record[0][19]

        # what happen if there are multiples? is it the same?
        self._npairs = self._parameters["npairs"]
        self._nsingles = self.ntot  - 3*self._npairs
        n = 2*self._npairs + self._nsingles
        self._n = n

        names = record[7]
        #select only stars, not centers of mass
        mask = names <= n

        self._stars["mass"] = record[1][mask]
        X = numpy.reshape(record[4],(3,self._ntot),order="F")
        self._stars["x"] = X[0,:][mask]
        self._stars["y"] = X[1,:][mask]
        self._stars["z"] = X[2,:][mask]
        X = numpy.reshape(record[5],(3,self._ntot),order="F")
        self._stars["vx"] = X[0,:][mask]
        self._stars["vy"] = X[1,:][mask]
        self._stars["vz"] = X[2,:][mask]
        self._stars["name"] = record[7][mask]

        self._time = self._parameters["time"]
        self._physical = False

    def to_physical(self):
        "Make sure stars are in physical units. Transform if not."
        if not self._physical :
            self._stars["mass"] *= self.parameters["zmbar"]*self.n
            self._stars["x"]    *= self.parameters["rbar"]
            self._stars["y"]    *= self.parameters["rbar"]
            self._stars["z"]    *= self.parameters["rbar"]
            self._stars["vx"]   *= self.parameters["vstar"]
            self._stars["vy"]   *= self.parameters["vstar"]
            self._stars["vz"]   *= self.parameters["vstar"]
            self._time *= self.parameters["tscale"]
            self._physical = True

    def to_nbody(self):
        "Make sure stars are in physical units. Transform if not."
        if self._physical :
            self._stars["mass"] /= self.parameters["zmbar"]*self.n
            self._stars["x"]    /= self.parameters["rscale"]
            self._stars["y"]    /= self.parameters["rscale"]
            self._stars["z"]    /= self.parameters["rscale"]
            self._stars["vx"]   /= self.parameters["vstar"]
            self._stars["vy"]   /= self.parameters["vstar"]
            self._stars["vz"]   /= self.parameters["vstar"]
            self._time /= self.parameters["tscale"]
            self._physical = False

    def to_center(self,center=None):
        """ Move stars to the specified center. If not specified, use the
            center of density calculated by Nbody6
        """
        if center is None:
            center = self.parameters["rdens"]
        self.stars["x"] -= center[0]
        self.stars["y"] -= center[1]
        self.stars["z"] -= center[2]
        self.parameters["rdens"][0] -= center[0]
        self.parameters["rdens"][1] -= center[1]
        self.parameters["rdens"][2] -= center[2]

    def reorder(self,isort=None):
        """ Sort stars according to the isort list (e.g. result of numpy.argsort(). If nothing specified, sort by name. """
        if isort is None:
           isort = numpy.argsort(self.stars["name"])

        for key in self.stars.keys():
           self.stars[key] = self.stars[key][isort]

def parse_inputfile(inputfilename):
    def parseline(line,resdict,names,dtypes):
        values = line.split()
        for name,val,dtype in zip(names,values,dtypes):
            resdict[name] = dtype(val)
    def parsekzline(line,resdict):
        values = line.split()
        for val in values :
            resdict["KZ"].append(int(val))

    result = dict(KZ=[ None ]) # Initial none is added so e.g. KZ[1] is actually KZ(1) and not KZ(2) as it would be in the python syntax
    inputfile = open(inputfilename,"r")

    parseline(inputfile.readline(),result,
        ["KSTART","TCOMP","TCRITp","isernb","iserreg","iserks"],[int,float,float,int,int,int])
    parseline(inputfile.readline(),result,
        ["N","NFIX","NCRIT","NRAND","NNBOPT","NRUN"],[int]*6)
    parseline(inputfile.readline(),result,
        ["ETAI","ETAR","RS0","DTADJ","DELTAT","TCRIT","QE","RBAR","ZMBAR"],[float]*9)

    #parse KZ options
    for _i in range(5): 
        parsekzline(inputfile.readline(),result)
    KZ = result["KZ"] 

    parseline(inputfile.readline(),result,
        ["DTMIN","RMIN","ETAU","ECLOSE","GMIN","GMAX","SMAX"],[float]*7)
    parseline(inputfile.readline(),result,
        ["ALPHA","BODY1","BODYN","NBIN0","NHI0","ZMET","EPOCH0","DTPLOT"],[float]*3+[int]*2+[float]*3)

    if KZ[5]==2 :
        parseline(inputfile.readline(),result,
            ["AP0","ECC","N2","SCALE"],[float,float,int,float])
    if KZ[5]==3 :
        parseline(inputfile.readline(),result,
            ["AP0","ECC","SCALE"],[float,float,float])
    if KZ[5]==4 :
        parseline(inputfile.readline(),result,
            ["SEMI","ECC","M1","M2"],[float]*4)
    if (KZ[5]==6 and KZ[24]<0) :
        parseline(inputfile.readline(),result,
            ["ZMH","RCUT"],[float]*2)

    parseline(inputfile.readline(),result,
            ["Q","VXROT","VZROT","RTIDE"],[float]*4)

    #Gradual formation KZ50 > 1
    if KZ[50] > 0 :
        parseline(inputfile.readline(),result,
                ["SFR", "N00", "DTGF"],[float]*3)
    if KZ[50] == 2 :
        parseline(inputfile.readline(),result,
                ["NT", "ITSTART", "GRIDSIZE", "DTT","THRESHOLD","TFOLDER"],[int,int,float,float,float,str])
        tfold = result["TFOLDER"] 
        result["TFOLDER"] = tfold.strip("\"").strip("\'")

    if KZ[14]==2 :
        parseline(inputfile.readline(),result,
                ["GMG","RG0"],[float]*2)
    if KZ[14]==3 :
        parseline(inputfile.readline(),result,
                ["GMG","DISK","A","B","VCIRC","RCIRC","GMB","AR","GAM","_R1","_R2","_R3","_V1","_V2","_V3"],[float]*15)
        result["RG"] = [result["_R1"],result["_R2"],result["_R3"]]
        result["VG"] = [result["_V1"],result["_V2"],result["_V3"]]
        toremove = ["_R1","_R2","_R3","_V1","_V2","_V3"]
        for k in toremove:
            result.pop(k,None)
    if (KZ[14] == 3 or KZ[14]==4):
        parseline(inputfile.readline(),result,
            ["MP","AP","MPDOT","TDELAY"],[float]*4)
    if KZ[14] == 5 :
        parseline(inputfile.readline(),result,
            ["KRHO", "MP", "AP", "MPDOT", "TDELAY", "BGSCALE"],[float]*5+[int])

    return result

def get_number_of_snapshots(folder):
    l=[x.replace("%sconf.3_"%folder,"") for x in glob.glob(folder+"conf.3*") ]
    l.sort(key=float)
    if len(l) == 0:
        raise ValueError("No snapshots in this folder.")
    return len(l)

def read_snapshot(folder,snapshot=0,inputfilename="input"):
    if folder[-1] != "/" : folder +="/"
    opt = parse_inputfile(folder+"/"+inputfilename)
    kz = opt["KZ"]
   
    l=[x.replace("%sconf.3_"%folder,"") for x in glob.glob(folder+"conf.3*") ]
    l.sort(key=float)
    if len(l) == 0:
        raise ValueError("No snapshots in this folder.")

    snapshotfile ="%s/conf.3_%s"%(folder,l[snapshot])
    inputfile ="%s/%s"%(folder,inputfilename)
    return Snapshot(snapshotfile,inputfile)
