# -*- coding: utf-8 -*-
import subprocess
import warnings
import os
import glob

from scipy.io import FortranFile
import numpy

from ._ParticleMethods import Methods

def parse_inputfile(inputfilename):
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

    result = dict(KZ=[ None ]) 
    # Initial none is added so e.g. KZ[1] is actually KZ(1) and not KZ(2) as it
    # would be in the python syntax
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
                ["NT", "ITSTART", "GRIDSIZE","DTT","THRESHOLD","MBG","BGVEL","TFOLDER"],[int,int,float,float,float,float,float,str])
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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bwdat = numpy.loadtxt(widefile,skiprows=2).T
        bdat = numpy.loadtxt(hardfile,skiprows=4).T
    widebin,hardbin = dict(),dict()
    for d,dat in zip( [widebin,hardbin],[bwdat,bdat] ):
        if len(dat) > 0 :
            if dat.ndim == 1:
                dat = numpy.array([ numpy.array(dat,dtype=float) ,]  ).T
            d["primary"  ] = numpy.array(dat[0,:],dtype=int)
            d["secondary"] = numpy.array(dat[1,:],dtype=int)
            d["m1"]      = dat[2,:] #Msun units
            d["m2"]      = dat[3,:] #MSun units
            d["ebin"]    = dat[4,:] #nbody units
            d["ecc"]     = dat[5,:] 
            d["period"]  = dat[6,:] #days
            d["semi"]    = dat[7,:] #AU
            d["kstar1"]  = dat[10,:]
            d["kstar2"]  = dat[11,:]
        else:
            d["primary"  ] = numpy.array([])
            d["secondary"] = numpy.array([])
            d["m1"]      = numpy.array([])
            d["m2"]      = numpy.array([])
            d["ebin"]    = numpy.array([])
            d["ecc"]     = numpy.array([])
            d["period"]  = numpy.array([])
            d["semi"]    = numpy.array([])
            d["kstar1"]  = numpy.array([])
            d["kstar2"]  = numpy.array([])

    return hardbin,widebin

class Snapshot(object):
    """
    WARNING: This description is outdated
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
    # TODO: Update Snapshot class documentation

    def __init__(self,snapshotfile,inputfile,singlefile=False,snapshot=1) :
        self._snapshotfile = snapshotfile
        self._inputfile = parse_inputfile(inputfile)
        self._singlefile = singlefile
        self.snapshot = snapshot
        self.__recordnumber = 0 if singlefile else snapshot
        self.__recordfile = None
        self.__parameters = dict()
        self._physical = False
        self.__read_snapshot()


    @property
    def parameters(self):
        """ Dict. containing parameters read from the AS variable.
        """
        return self.__parameters

    @property
    def stars(self):
        """ Dict. containing stellar coordinates:
        mass, x,y,z, vx,vy,vz,epot
        """
        return self.__allstars[ (self.__allstars.name <= self.n)
                               *(self.__allstars.name > 0)  ]
    @property
    def allparticles(self):
        """ Dict. containing stellar coordinates:
        mass, x,y,z, vx,vy,vz,epot
        """
        return self.__allstars[ (self.__allstars.name > 0)  ]

    @property
    def unresolved_stars(self):
        """
        Returns Particles object with regularized binaries as center of mass
        particles.
        """
        #return self.__allstars[self.__allstars.name > 2*self.npairs]
        return self.__allstars[(2*self.npairs+1):]

    @property
    def bound_stars(self):
        """
        Bound set of particles
        """
        return self.__bound_stars(True)

    @property
    def bound_stars_unresolved(self):
        """
        Bound set of stars, with unresolved hard binaries
        """
        return self.__bound_stars(False)

    @property
    def unbound_stars(self):
        """
        Bound set of particles
        """
        return self.__unbound_stars(True)

    @property
    def unbound_stars_unresolved(self):
        """
        Bound set of stars, with unresolved hard binaries
        """
        return self.__unbound_stars(False)

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
        """ bool, True : stars are in physical units: Msun, parsec, km/s
              False: stars are in code units
        """
        return self._physical

    @property
    def singlefile(self):
        """ bool, True : All snapshot are stored in one file (use step to advance)
                  False: Snapshot in different files. Provide a file to advance.
        """
        return self._singlefile

    def __bound_stars(self,resolved=True):
        """Return bound subset of stars.
        arguments:
            resolved : if true, all binaries are resolved. if False, center of
            mass particles are returned. Default: True
        """
        stars = self.unresolved_stars
        bmask = abs(stars.epot + stars.pot)*stars.mass \
                >0.5*stars.mass*numpy.array(stars.vx**2+stars.vy**2+stars.vz**2) 
        bound_unresolved = stars[bmask]
        if not resolved:
            return bound_unresolved
        else:
            return self.resolve_set(bound_unresolved)

    def __unbound_stars(self,resolved=True):
        """Return unbound subset of stars.
        arguments:
            resolved : if true, all binaries are resolved. if False, center of
            mass particles are returned. Default: True
        """
        stars = self.unresolved_stars
        bmask = abs(stars.epot + stars.pot)*stars.mass \
                <=0.5*stars.mass*stars.v**2 #numpy.array(stars.vx**2+stars.vy**2+stars.vz**2) 
        unbound_unresolved = stars[bmask]
        if not resolved:
            return unbound_unresolved
        else:
            return self.resolve_set(unbound_unresolved)

    def virial_energy_of_set(self,particles):
        """
        Compute the virial energy VIR = SUM F·r. 
        Not all options are implemented here. Please update to your needs.
        If no background potential is used, VIR := POT for an Nbody system,
        however this is not true for other options.

        It relies on the potential calculated by nbody6 and the external
        potential from the function Snapshot.external_potential_at_point() that
        should be well defined for each case.

        Currently only power law potential is implemented (from my customized
        version of Nbody6). 
        TODO: Update and add other backgrounds
        """
        result = -(particles.pot*particles.mass).sum()*0.5
        if self.inputfile["KZ"][14] == 5 :
            krho = self.inputfile["KRHO"] 
            Rcl = self.inputfile["AP"]
            if particles.physical:
                Rcl*=self.parameters["rbar"]
            GMR = self.external_potential_at_point([0],[0],[Rcl],
                                                   physical=particles.physical
                                                   )[0]
            r = numpy.sqrt(particles.x**2 + particles.y**2 + particles.z**2)

            mask = r<=Rcl 
            # K = (3.0-krho) / (2.-krho)
            # result+=numpy.nansum( (2.0-krho) / (5.-2.*krho) *(
                    # (particles.epot+K*GMR)*particles.mass
                    # )[mask]
                    # )
            #result += numpy.nansum( (GMR*(krho-3)*particles.mass[mask] - 
            #         (2.0 - krho)*particles.epot[mask]*particles.mass[mask]
            #         ))  
            result += numpy.nansum(
                  #-GMR*(r[mask]/Rcl)**(2.0-krho)*particles.mass[mask]
                  -GMR*(r[mask]/Rcl)**(2.0-krho*0.5)*particles.mass[mask]
                  )

            mask = r>=Rcl 
            #result += numpy.nansum( -(particles.epot*particles.mass)[mask] )
            result += numpy.nansum( -(particles.epot*particles.mass)[mask] )

        return result

    def step(self,n=1):
        """ 
        Advance to next snapshot. 
        Usefull when there is a single snapshotfile.
        """
        if not self._singlefile:
            raise Exception("Snapshot.step function only works with singlefile=True")
        self.snapshot += n
        self.__read_snapshot()
        pass

    def __read_record(self):
      inttype = numpy.int32
      floattype = numpy.float32
      kz = self.inputfile["KZ"]
      NTOT,MODEL,NDRUN,NK = tuple(self.__recordfile.read_ints(dtype=inttype))
      self._ntot = NTOT
      if not kz[50] == 1 : 
         record = self.__recordfile.read_record(
                 numpy.dtype( (floattype,NK)   ) #AS
                ,numpy.dtype( (floattype,NTOT) ) #bodys
                ,numpy.dtype( (floattype,NTOT) ) #rhos
                ,numpy.dtype( (floattype,NTOT) ) #xns
                ,numpy.dtype( (floattype,(NTOT*3) ) )#x
                ,numpy.dtype( (floattype,(NTOT*3) ) )#v
                ,numpy.dtype( (floattype,NTOT  ) )#phi
                ,numpy.dtype( (inttype  ,NTOT) ) #name
              )
      if kz[50] == 1 : 
         record = self.__recordfile.read_record(
                 numpy.dtype( (floattype,NK)   ) #AS
                ,numpy.dtype( (floattype,NTOT) ) #bodys
                ,numpy.dtype( (floattype,NTOT) ) #rhos
                ,numpy.dtype( (floattype,NTOT) ) #xns
                ,numpy.dtype( (floattype,(NTOT*3) ) )#x
                ,numpy.dtype( (floattype,(NTOT*3) ) )#v
                ,numpy.dtype( (floattype,NTOT  ) )#phi
                ,numpy.dtype( (inttype,NTOT) ) #name
                ,numpy.dtype( (inttype,NTOT) ) #kstar
              )
      return record

    def __read_snapshot(self):
      if  self.__recordfile is None:
          self.__recordfile = FortranFile(self._snapshotfile,"r") 
      for i in range(self.__recordnumber,self.snapshot+1):
          record = self.__read_record()
          self.__recordnumber += 1
      physical = self.physical
      self.__structure(record)
      #self.__record = record
      if physical:
          self.to_physical()

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
        self.__parameters["rc"] = float(record[0][12])
        self.__parameters["nc"] = int(record[0][13])
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
        self._time = self.__parameters["time"]
        self._physical = False

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
        stars_dict["pot"] = record[6]
        if self.inputfile["KZ"][50] == 1:
            stars0 = (stars_dict["name"] <= 150) & (stars_dict["name"] > 0)
            self.__mgas0 = self.inputfile["MP"]-stars_dict["mass"][stars0].sum()
        stars_dict["epot"] = self.external_potential_at_point(X[0,:],X[1,:],X[2,:])
        if self.inputfile["KZ"][50] ==1 :
            stars_dict["kstar"] = record[8]
        else:
            stars_dict["kstar"] = stars_dict["name"] * 0 - 999 #not known

        allparticles =  Particles(stars_dict,center=self.parameters["rdens"])
        allparticles.I = numpy.arange(1,self.ntot+1,1)

        self.__stars = allparticles[mask_stars]
        self.__allstars = allparticles

    def external_potential_at_point(self,x,y,z,physical=None):
        """ Return background potential field at given point.
        physical : if True x,y,z are on pc. And result in km**2 * s**-1
                   if False in nbody units
                   if None : taken from the Parent snapshot class
                   default: None
        
        returns:
          background potential at x,y,z
          if physical is True : result in km**2 * s**-2
          else: result in Nbody units
        """
        if physical is None:
            physical = self.physical
        x = numpy.array(x)
        y = numpy.array(y)
        z = numpy.array(z)

        G = 1.0
        if "MP" in self.inputfile.keys():
            mp = self.inputfile["MP"]
            if self.inputfile["KZ"][50]==1 :
                mp = self.__mgas0 
            Mgas = mp - self.inputfile["MPDOT"]*(self.parameters["time"] - self.inputfile["TDELAY"] )
            Rcore = self.inputfile["AP"]
            if physical : 
                Mgas *= self.parameters["zmbar"]
                Rcore *= self.parameters["rbar"]
                G = 4.3020077853E-3
            if self.inputfile["KZ"][14] == 5 and Mgas > 0 : #TODO implement other external potentials
                #truncated power law potential
                #print(self.parameters["time"],Mgas,end="\n")
                r = numpy.sqrt(x**2 + y**2 + z**2 ) 
                krho = self.inputfile["KRHO"]
                #result = - G*Mgas*(r/Rcore)**(3.0 - self.inputfile["KRHO"] )/r
                result = G*Mgas*((r/Rcore)**(2-krho) - 3.0 + krho)/Rcore/(2-krho)
                result[ r >= Rcore ] = - G * Mgas / r[r >= Rcore]
            else:
                result = x*0.0
        else:
            result = x*0.0
        return result

    def resolve_set_old(self, unresolved_set ):
        """ 
        Old version: Worked and tested but slow
        Resolve binaries from unresolved_set, based on names.
        unresolved_set must be a set of stars taken from the same snapshot
        object containing this function.
        """
        def where(a):
            return numpy.where(a==self.stars.name)[0][0]
        where = numpy.vectorize(where)
        cm_name = (unresolved_set.name - self.n)[unresolved_set.name > self.n]
        iprim = where(cm_name)
        ising = where(unresolved_set[ (unresolved_set.name <= self.n) &
                                       (unresolved_set.name > 0) ].name)
        iall = numpy.concatenate([numpy.dstack( (iprim,iprim+1) ).flatten(),
                                   ising])
        return self.stars[iall]

    def resolve_set(self, unresolved_set, split_set=False ):
        """ 
        New version: under testing. Should be faster.
        Resolve binaries from unresolved_set, based on names.
        unresolved_set must be a set of stars taken from the same snapshot
        object containing this function.

        Input:
        -----
        unresolved_set : Particles object
                         Must be a Particle object contained in this Snapshot.

        split_set: bool (default:False)
                   if True returns tuple of Particles with single, primary and 
                   secondary set of stars

        Output:
        ------
        resolved_particles  : if split_set == False
        single,primary,secondary  : Tuple of Particles set if split_set==True
                                    primary and secondary order is equivalent.
        
        Example:
        -------
        >> sn = Snapshot(folder,input)
        >> single,primary,secondary = sn.resolve_set(sn.stars,split_set=True)
        >> bin_sep = numpy.sqrt((primary.x - secondary.x)**2 
                             +  (primary.y - secondary.y)**2 
                             +  (primary.z - secondary.z)**2
                               )
        Note: if resolved_particles is already resolved, primary and secondary
              sets will result on empty Particles objects.
        """
        primary_names = (unresolved_set.name - self.n)[unresolved_set.name > self.n]
        singles_names = unresolved_set.name[unresolved_set.name <= self.n] 

        main_indexes = numpy.arange(len(self.stars))
        iprim = main_indexes[numpy.isin(self.stars.name,primary_names)]
        ising =  main_indexes[numpy.isin(self.stars.name,singles_names) ]

        if not split_set:
            iall = numpy.concatenate([numpy.dstack( (iprim,iprim+1) ).flatten(),
                                       ising])
            return self.stars[iall]
        else:
            #print(primary_names)
            #print(main_indexes)
            print(singles_names)
            print(ising)
            print(iprim)
            return self.stars[ising], self.stars[iprim], self.stars[iprim+1]

    def to_physical(self):
        "Make sure stars are in physical units. Transform if not."
        if not self._physical :
            self.__allstars.to_physical(self.parameters)
            self._time *= self.parameters["tscale"]
            self.parameters["rdens"]*= self.parameters["rbar"]
            self._physical = True

    def to_nbody(self):
        "Make sure stars are in physical units. Transform if not."
        if self._physical :
            self.__allstars.to_nbody(self.parameters)
            self._time /= self.parameters["tscale"]
            self.parameters["rdens"] /= self.parameters["rbar"]
            self._physical = False

    def to_center(self,center=None):
        """ Move stars to the specified center. If not specified, use the
            center of density calculated by Nbody6
        """
        if center is None:
            center = self.parameters["rdens"]
        self.__allstars.to_center(center)
        self.parameters["rdens"][0] -= center[0]
        self.parameters["rdens"][1] -= center[1]
        self.parameters["rdens"][2] -= center[2]

    def reorder(self,isort=None):
#TODO (check): Snaphsot.reorder should not be needed anymore 
        """ Sort stars according to the isort list (e.g. result of numpy.argsort(). If nothing specified, sort by name. """
        if isort is None:
           isort = numpy.argsort(self.stars["name"])

        for key in self.stars.keys():
           self.stars[key] = self.stars[key][isort]

    def unresolve_all(self):
#TODO : Check Snapshot.unrsolve_all is still necessary
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
        higher_orders = numpy.where( (pairs > self.n) | ( pairs <= 0) )  
        if len(pairs) == 0:
            return self.stars
        pairs = numpy.delete(pairs,higher_orders[0] , axis=0 )
        primaries,secondaries = numpy.array(pairs).T

        #extract binaries in multiples and add to singles for now 

        #position of singles, primaries and secondaries
#TODO: handle multiple systems
        i_singles = numpy.invert(numpy.isin(self.stars.name,allmembers))
        singles = self.stars[i_singles]

        i_prim = numpy.isin(self.stars.name,primaries)
        prim_stars = self.stars[i_prim]

        #i_sec = numpy.isin(self.stars.name,secondaries)
        #sec_stars = self.stars[i_sec]

        #TODO:  improve this, must be a better way
        isec = []
        for prim,sec  in  pairs:
            #print(numpy.where(self.stars.name == sec)[0][0] )
            ii = numpy.where(self.stars.name == sec)[0][0]
            #print("[",prim,sec,"]", self.stars.name[ii],  ii)
            isec.append( ii )
        sec_stars = self.stars[ numpy.array(isec) ] 

        #Make sure they are in the correct order
        #isort = numpy.argsort(primaries)
        #primaries = primaries[isort]
        #secondaries = secondaries[isort]

        #order = numpy.argsort(prim_stars)
        #TODO: add check that pairs are correct
        bdict = dict()
        bdict["name"] = prim_stars.name + self.n #TODO: check this is the standard
        bdict["mass"] = prim_stars.mass + sec_stars.mass
        bdict["x"] = (prim_stars.mass*prim_stars.x 
                      + sec_stars.mass*sec_stars.x) / bdict["mass"]
        bdict["y"] = (prim_stars.mass*prim_stars.y 
                      + sec_stars.mass*sec_stars.y) / bdict["mass"]
        bdict["z"] = (prim_stars.mass*prim_stars.z 
                      + sec_stars.mass*sec_stars.z) / bdict["mass"]
        bdict["vx"] = (prim_stars.mass*prim_stars.vx 
                       + sec_stars.mass*sec_stars.vx) / bdict["mass"]
        bdict["vy"] = (prim_stars.mass*prim_stars.vy
                       + sec_stars.mass*sec_stars.vy) / bdict["mass"]
        bdict["vz"] = (prim_stars.mass*prim_stars.vz 
                       + sec_stars.mass*sec_stars.vz) / bdict["mass"]
        bdict["pot"] = (prim_stars.mass*prim_stars.pot
                        + sec_stars.mass*sec_stars.pot) / bdict["mass"]
        bdict["epot"] = self.external_potential_at_point(bdict["x"],bdict["y"],
                                                         bdict["z"])

        result_dict = dict()
        for key in singles.available_attributes():
            sval = singles[key]
            bval = bdict[key]
            result_dict[key] = numpy.concatenate( [sval,bval] )

        return Particles(result_dict,
                         physical=self.physical)

class Particles(Methods):
    """
    Object to store and manipulate set of particles.
    This object is not supposed to be used by itself, but
    wrapped by the Snapshot object.
    """
    def __init__(self,stars_dict,center=[0.,0.,0.],
                 physical = False):
        self.__n = len(stars_dict["name"])  if hasattr(stars_dict["name"],"__len__" ) else 1  #must be first parameter to be setted
        self.__data = stars_dict
        self.__physical = physical
        l=[]
        for key in stars_dict:
            if not hasattr(stars_dict[key],"__len__"):
                stars_dict[key] = [stars_dict[key]] 
            setattr(self,key,numpy.array( stars_dict[key] ))

        self.__sanity_check()
        self.__index = 0
        self.__center = numpy.array(center)
        self.__center_velocity = numpy.array([0.,0.,0.])
        self.__r = None
        self.__v = None

    @property
    def center(self):
        return self.__center
    @property
    def center_velocity(self):
        return self.__center_velocity

    @property
    def physical(self):
        return self.__physical

    @property
    def r(self):
        if self.__r is None:
            self.__r = numpy.sqrt( (self.x - self.center[0])**2 
                                  +(self.y - self.center[1])**2  
                                  +(self.z - self.center[2])**2)
        return self.__r

    @property
    def v(self):
        if self.__v is None:
            self.__v = numpy.sqrt( (self.vx - self.center_velocity[0])**2 
                                  +(self.vy - self.center_velocity[1])**2  
                                  +(self.vz - self.center_velocity[2])**2)
        return self.__v 


    def to_physical(self,parameters):
        """
        Make sure stars are in physical units. Transform if not.
        input:
            parameters : dictionary with scaling parameters, 
            zmbar for mass, rbar for length and vstar for velocities.
            Normally taken from Snapshot object as:
                sn = Reader.read_snapshot("./",0)
                parameters = sn.parameters
        """
        if not self.physical:
            self.mass *= parameters["zmbar"]
            self.x    *= parameters["rbar"]
            self.y    *= parameters["rbar"]
            self.z    *= parameters["rbar"]
            self.vx   *= parameters["vstar"]
            self.vy   *= parameters["vstar"]
            self.vz   *= parameters["vstar"]
            G = 4.3020077853E-3
            escale = G*parameters["zmbar"]/parameters["rbar"]
            self.pot*=escale
            self.epot*=escale
            self.center *= parameters["rbar"]
            self.__physical = True

    def to_nbody(self,parameters):
        """
        Make sure stars are in physical units. Transform if not.
        input:
            parameters : dictionary with scaling parameters, 
            zmbar for mass, rbar for length and vstar for velocities.
            Normally taken from Snapshot object as:
                sn = Reader.read_snapshot("./",0)
                parameters = sn.parameters
        """
        if self.physical:
            self.mass /= parameters["zmbar"]
            self.x    /= parameters["rscale"]
            self.y    /= parameters["rscale"]
            self.z    /= parameters["rscale"]
            self.vx   /= parameters["vstar"]
            self.vy   /= parameters["vstar"]
            self.vz   /= parameters["vstar"]
            G = 4.3020077853E-3
            escale = G*parameters["zmbar"]/parameters["rbar"]
            self.epot /= escale
            self.pot /=escale
            self.__physical = False

    def to_center(self,center=None):
        if center is None:
            center = self.center
        else:
            center = numpy.array(center)
        self.x -= center[0]
        self.y -= center[1]
        self.z -= center[2]
        self.__center -= center

    def to_center_velocity(self,center_velocity=None):
        if center_velocity is None:
            center_velocity = self.center_velocity
        else:
            center_velocity = numpy.array(center_velocity)
        self.vx -= center_velocity[0]
        self.vy -= center_velocity[1]
        self.vz -= center_velocity[2]
        self.__center_velocity -= center_velocity

    def set_center(self,center):
        print("setting center",center)
        if len(center) == 3:
            self.__center = numpy.array(center,dtype=numpy.float32)
            self.__r = None #need to be recalulated
        else:
            raise ValueError("Particles center must have length 3")

    def copy(self):
        outputdict = self.__data.copy()
        # for key  in self.__data.keys() :
            # outputdict[key] = self.__data[key].copy()
        return Particles(outputdict,center=self.center,physical=self.physical)

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
        return Particle(d,physical=self.physical)

    def __getitem__(self,index) :
        if type(index) == str :
            return self.__data[index]
        else:
            d = dict()
            for k in self.__data:
                d[k] = self.__data[k][index]
            return Particles(
                d,center=self.center,
                physical=self.physical
                )

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
        #TODO Format string representation of Particle Class
        return str(self.__data)


    def available_attributes(self):
        return list(self.__data.keys() )

class Particle(object):
    def __init__(self,star_dict,physical=1):
        for key in star_dict:
            setattr(self,key,star_dict[key])
        self.__data = star_dict
        self.__physical
    @property
    def physical(self):
        return self.__physical
    def __str__(self):
        return "Particle Object: %s "% str(self.__data)
    def __repr__(self):
        return "Particle Object: %s "% str(self.__data)
