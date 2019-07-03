import subprocess
import numpy
from scipy.io import FortranFile
import os
import glob

from ._ParticleMethods  import Methods

def parse_inputfile(inputfilename,**kw):
    """
    parse Nbody6 inputfile 'inputfilename' into a dictionary.
    Include options in the customized Nbody6GF version.
    """
    # TODO : this function should be a well defined object able to write and
    # modify an inputfile consistenly, e.g. read an inputfile, modify one
    # parameter (e.g. KZ option) and write a new inputfile taking care and
    # raising exception if modifications are inconsistent.
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
                ["NT", "ITSTART", "GRIDSIZE", "DTT","THRESHOLD","BGVEL","TFOLDER"],[int,int,float,float,float,float,str])
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

def get_binaries_from_files(hardfile,widefile):
    bwdat = numpy.loadtxt(widefile,skiprows=2).T
    bdat = numpy.loadtxt(hardfile,skiprows=4).T
    widebin,hardbin = dict(),dict()
    for d,dat in zip( [widebin,hardbin],[bwdat,bdat] ):
        d["primary"  ] = numpy.array(dat[0,:],dtype=int)
        d["secondary"] = numpy.array(dat[1,:],dtype=int)
        d["ebin"]      = dat[4,:] #nbody units
        d["ecc"]       = dat[5,:] 
        d["period"]    = dat[6,:] #days
        d["semi"]      = dat[7,:] #AU
    return hardbin,widebin

class Snapshot(object):
    """
    An object containing the information of a snapshot of Nbody6.
    Contains :
      self.stars           : Datamodel.Particles object conatining the stellar data.
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
        self._snapshotfile = snapshotfile
        self._inputfile = parse_inputfile(inputfile)

        record = self.__read_snapshot(snapshotfile)
        self.__parameters = dict()
        self.__structure(record)
        self.__record = record

    @property
    def parameters(self):
        """ Dict. containing parameters read from the AS variable.
        """
        return self.__parameters

    @property
    def stars(self):
        """ Dict. containing stellar coordinates:
        mass, x,y,z, vx,vy,vz
        """
        return self.__allstars[self.__allstars.name <= self.n]

    @property
    def unresolved_stars(self):
        """
        Returns Particles object with regularized binaries as single stars
        """
        return self.__allstars[self.__allstars.name > 2*self.npairs]

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


    def __read_snapshot(self,name):
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

    def __structure(self,record):
        self.__parameters["time"] = record[0][0]
        self.__parameters["npairs"] = int(record[0][1])
        self.__parameters["rbar"] = record[0][2]
        self.__parameters["zmbar"] = record[0][3]
        self.__parameters["rtide"] = record[0][4]
        self.__parameters["tidal4"] = record[0][5]
        self.__parameters["rdens"] = numpy.array(record[0][6:9])
        self.__parameters["ntcr"] = record[0][9] # number of elapsed initial crossing times 
        self.__parameters["tscale"] = record[0][10]
        self.__parameters["vstar"] = record[0][11]
        self.__parameters["rc"] = record[0][12]
        self.__parameters["nc"] = record[0][13]
        self.__parameters["vc"] = record[0][14]
        self.__parameters["rhom"] = record[0][15]
        self.__parameters["cmax"] = record[0][16]
        self.__parameters["rscale"] = record[0][17] #half mass radius
        self.__parameters["rsmin"] = record[0][18]
        self.__parameters["dmin1"] = record[0][19]

        # what happen if there are multiples? is it the same?
        self._npairs = self.__parameters["npairs"]
        self._nsingles = self.ntot  - 3*self._npairs
        n = 2*self._npairs + self._nsingles
        self._n = n

        names = record[7]
        #select only stars, not centers of mass
        mask_stars = names <= n
        #mask_centers_of_mass = names > n

        stars_dict = dict()
        stars_dict["name"] = record[7]
        stars_dict["mass"] = record[1]
        X = numpy.reshape(record[4],(3,self._ntot),order="F")
        stars_dict["x"] = X[0,:]
        stars_dict["y"] = X[1,:]
        stars_dict["z"] = X[2,:]
        X = numpy.reshape(record[5],(3,self._ntot),order="F")
        stars_dict["vx"] = X[0,:]
        stars_dict["vy"] = X[1,:]
        stars_dict["vz"] = X[2,:]

        allparticles =  Particles(stars_dict,center=self.parameters["rdens"])

        self.__stars = allparticles[mask_stars]
        self.__allstars = allparticles

        self._time = self.__parameters["time"]
        self._physical = False

    def to_physical(self):
        "Make sure stars are in physical units. Transform if not."
        if not self._physical :
            self.__allstars.mass *= self.parameters["zmbar"]
            self.__allstars.x    *= self.parameters["rbar"]
            self.__allstars.y    *= self.parameters["rbar"]
            self.__allstars.z    *= self.parameters["rbar"]
            self.__allstars.vx   *= self.parameters["vstar"]
            self.__allstars.vy   *= self.parameters["vstar"]
            self.__allstars.vz   *= self.parameters["vstar"]
            self.__allstars.center *= self.parameters["rbar"]
            self._time *= self.parameters["tscale"]
            self.parameters["rdens"]*= self.parameters["rbar"]
            self._physical = True

    def to_nbody(self):
        "Make sure stars are in physical units. Transform if not."
        if self._physical :
            self.__allstars.mass /= self.parameters["zmbar"]
            self.__allstars.x    /= self.parameters["rscale"]
            self.__allstars.y    /= self.parameters["rscale"]
            self.__allstars.z    /= self.parameters["rscale"]
            self.__allstars.vx   /= self.parameters["vstar"]
            self.__allstars.vy   /= self.parameters["vstar"]
            self.__allstars.vz   /= self.parameters["vstar"]
            self._time /= self.parameters["tscale"]
            self.parameters["rdens"] /= self.parameters["rbar"]
            self._physical = False

    def to_center(self,center=None):
        """ Move stars to the specified center. If not specified, use the
            center of density calculated by Nbody6
        """
        if center is None:
            center = self.parameters["rdens"]
        #self.__allstars.x -= center[0]
        #self.__allstars.y -= center[1]
        #self.__allstars.z -= center[2]
        self.__allstars.to_center(center)
        self.parameters["rdens"][0] -= center[0]
        self.parameters["rdens"][1] -= center[1]
        self.parameters["rdens"][2] -= center[2]
        #self.__stars.__center = center
        #self.__allstars.center = self.parameters["rdens"]

    def reorder(self,isort=None):
        """ Sort stars according to the isort list (e.g. result of numpy.argsort(). If nothing specified, sort by name. """
        if isort is None:
           isort = numpy.argsort(self.stars["name"])

        for key in self.stars.keys():
           self.stars[key] = self.stars[key][isort]

    def unresolve_all(self):
        """ 
        Unresolve all binaries using data from bwdat.19_* and bdat.9_*. 
        Only possible if KZ(18) >= 2. 
        returns a particle set with binaries as center of mass.
        If only binary parameters are needed, use the wrapper  Reader.read_binaries.
        """
        if self.inputfile["KZ"][8] < 2:
            raise ValueError("Need to run simulation with KZ[8] >= 2 to use this function")
        widefile = self._snapshotfile.replace("conf.3_","bwdat.19_")
        hardfile = self._snapshotfile.replace("conf.3_","bdat.9_")

        hard,wide = get_binaries_from_files(hardfile,widefile)
        ### construct the particle set
        #all names that belong to a binary
        allmembers = numpy.concatenate((hard["primary"],hard["secondary"],
                                        wide["primary"],wide["secondary"]) )
        primaries = numpy.concatenate((hard["primary"],wide["primary"] ))
        secondaries = numpy.concatenate((hard["secondary"],wide["secondary"] ))
        pairs = numpy.array(list( zip(primaries,secondaries)) ) 
        higher_orders = numpy.where( pairs > self.n  )
        #print("ntot",self.n)
        #print("higher_orders",higher_orders)#sorted(higher_orders[0],reverse=True))
        #print("higher",sorted( set(higher_orders[0]),reverse=True) )
        #print("hpairs",numpy.array(pairs) [higher_orders[0]] )
        #for i in sorted(set(higher_orders[0]),reverse=True)  :
        #    pairs.pop(i)
        pairs = numpy.delete(pairs,higher_orders[0] , axis=0 )
        primaries,secondaries = numpy.array(pairs).T

        #extract binaries in multiples and add to singles for now #TODO handle multiples into one single particle

        #position of singles, primaries and secondaries
        #TODO: handle multiple systems
        i_singles = numpy.invert(numpy.isin(self.stars.name,allmembers))
        singles = self.stars[i_singles]

        i_prim = numpy.isin(self.stars.name,primaries)
        prim_stars = self.stars[i_prim]

        i_sec = numpy.isin(self.stars.name,secondaries)
        sec_stars = self.stars[i_sec]

        #print(i_prim.sum())
        #print(i_sec.sum())



        #Make sure they are in the correct order
        #isort = numpy.argsort(primaries)
        #primaries = primaries[isort]
        #secondaries = secondaries[isort]

        #order = numpy.argsort(prim_stars)
        #TODO: add check that pairs are correct
        bdict = dict()
        bdict["name"] = prim_stars.name + self.n #TODO: check this is the standard
        bdict["mass"] = prim_stars.mass + sec_stars.mass
        bdict["x"]    = (prim_stars.mass*prim_stars.x + sec_stars.mass*sec_stars.x) /bdict["mass"]
        bdict["y"]    = (prim_stars.mass*prim_stars.y + sec_stars.mass*sec_stars.y) /bdict["mass"]
        bdict["z"]    = (prim_stars.mass*prim_stars.z + sec_stars.mass*sec_stars.z) /bdict["mass"]
        bdict["vx"]   = (prim_stars.mass*prim_stars.vx + sec_stars.mass*sec_stars.vx) /bdict["mass"]
        bdict["vy"]   = (prim_stars.mass*prim_stars.vy + sec_stars.mass*sec_stars.vy) /bdict["mass"]
        bdict["vz"]   = (prim_stars.mass*prim_stars.vz + sec_stars.mass*sec_stars.vz) /bdict["mass"]

        result_dict = dict()
        for key in singles.available_attributes():
            sval = singles[key]
            bval = bdict[key]
            result_dict[key] = numpy.concatenate( [sval,bval] )

        return Particles(result_dict)

class Particles(Methods,object):
    """
    Object to store and manipulate set of particles.
    This object is not supposed to be used by itself, but
    wrapped by the Snapshot object.
    """
    def __init__(self,stars_dict,center=[0,0,0]):
        self.__n = len(stars_dict["name"]) #must be first parameter to be setted
        self.__data = stars_dict
        l=[]
        for key in stars_dict:
            if not hasattr(stars_dict[key],"__len__"):
                stars_dict[key] = [stars_dict[key]] 
            setattr(self,key,numpy.array( stars_dict[key] ))

        self.__sanity_check()
        self.__index = 0
        self.__center = numpy.array(center)

    @property
    def center(self):
        return self.__center

    def to_center(self,center=None):
        if center is None:
            center = self.center
        else:
            center = numpy.array(center)
        self.x -= center[0]
        self.y -= center[1]
        self.z -= center[2]
        self.__center -= center


    def set_center(self,center):
        print("setting center",center)
        if len(center) == 3:
            self.__center = numpy.array(center,dtype=numpy.float32)
        else:
            raise ValueError("Particles center must have length 3")

    def __sanity_check(self):
        l=[]
        for key in self.__data:
            l.append(len(self.__data[key]) ) 
        if len(set(l)) != 1 :
            raise ValueError("Not all parameters have the same length")

    def __iter__(self):
        return self

    def __next__(self):
        i = self.__index
        self.__index += 1
        if i > self.__n-1:
            self.__index = 0
            raise StopIteration
        d = dict()
        for k in self.__data:
            d[k] = self.__data[k][i]
        return Particle(d)

    def __getitem__(self,index) :
        if type(index) == str :
            return self.__data[index]
        else:
            d = dict()
            for k in self.__data:
                d[k] = self.__data[k][index]
            return Particles(d,center=self.center)

    def __setitem__(self,index,value):
        vlen = 1 if not hasattr(value,"__len__") else len(value)
        if index in self.__data.keys():
            if vlen == len(self):
                self.__data[index] = value
                self.__dict__[index] = value
            else:
                raise ValueError(" length of value must be the same as the number of particles" )
        else:
            raise KeyError("%s not in storage"%index)

    def __setattr__(self,key,value):
        particle_attribute = (not "_" in key) and ( not "center"  in key)
        if particle_attribute and len(self.__dict__ ) != 0  :
            vlen = 1 if not hasattr(value,"__len__") else len(value)
            if vlen != len(self) :
                raise ValueError("length of value  must be the same as the number of particles" )
        self.__dict__[key] = value
        if particle_attribute :
            self.__data[key] = value

    def __len__(self):
        return self.__n

    def __str__(self):
        #TODO Format this properly
        return str(self.__data)

    def available_attributes(self):
        return list(self.__data.keys() )

class Particle(object):
    def __init__(self,star_dict):
        for key in star_dict:
            setattr(self,key,star_dict[key])
        self.__data = star_dict
    def __str__(self):
        return "Particle Object: %s "% str(self.__data)


        
