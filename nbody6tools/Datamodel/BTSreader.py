import h5py
import numpy 
import os
from nbody6tools import Reader

#def get_scale_factors(inputpar):
class Converter(object):
    def __init__(self,inputfile):
        self.inputpar = Reader.parse_inputfile(inputfile)
        self.initialize()

    def initialize(self):
        #from units.f  of the Nbody6 code

        #Define GM & PC in cgs units and #AU in pc (2008 IAU values).
        GM = 6.6743e-08*1.9884e+33
        PC = 3.0856776e+18
        self.RBAR = self.inputpar["RBAR"]
        self.ZMBAR = self.inputpar["ZMBAR"]
        N = self.inputpar["N"]
        self.ZMASS = self.ZMBAR * N

        #Form time scale in seconds and velocity scale in km/sec.
        self.TSTAR = numpy.sqrt(PC/GM)*PC

        #Convert time scale from units of seconds to million years.
        self.TSTAR = self.TSTAR/(3.15e+07*1.0e+06)
        if (self.inputpar["KZ"][22] != 10) :
             self.VSTAR = 1.0e-05*numpy.sqrt(GM/PC)

        #Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
             if (self.ZMBAR <= 0.0) :
                 raise NotImplemented("ZMASS == 0 not implemented")
                 #ZMBAR = N/ZMASS
             if (self.RBAR <= 0.0):
                 raise NotImplemented("RBAR == 0 not implemented")
                 #RBAR = 1.0
    
        #    Scale to working units of RBAR in pc & ZMBAR in solar masses.
             self.TSTAR = self.TSTAR*numpy.sqrt(self.RBAR**3/(self.ZMASS*self.ZMBAR))
             self.VSTAR = self.VSTAR*numpy.sqrt(self.ZMASS*self.ZMBAR/self.RBAR)
        else :
        #USE scale factor from original input data instead
             self.TSTAR = self.TSTAR*numpy.sqrt(self.RBAR**3/self.ZMBAR) 
        #END IF

class H5nb6xxSnapshot(object):
    """ 
    Implement a reader for the Block Timestep Storage introduced by Cai et al.
    (2015) http://dx.doi.org/10.1088/0067-0049/219/2/31 
    Most of the code below is adapted from a script provided by Maxwell Cai.

    All units within this script are in Nbody units. User must convert to
    physical units using an appropriate converter.

    Parameters
    ----------
    inputfile : string
                Path to file with input parameters.
    snapshotfiles : list
                  List of filenames with the snapshots. This due that is not
                  guaranteed that all stars are stored in the same file or
                  adjacent files.  So, to be sure, the full list should be
                  provided.

    GOALS:
        Main goal:
            1) retrieve a Particles object with positions at an arbitrary time
            2) give an option to obtain resolved or unresolved binaries
            3) retrieve binary properties
            4) Integrate with snapshot object in the main package.
        Secondary goal:
            2) be able to advance step by step retrieving only current active
            particles.

        Note: The current version of Nbody6 do not output the total number of
        particles. I will make this reader to collect stars without knowing the
        total number. But identify them and collect their information in an ID array.

    """
    def __init__(self,snapshotfiles,inputfile,time=0):
        self.inputpar = Reader.parse_inputfile(inputfile)
        if self.inputpar["KZ"][46] != 1:
            raise NotImplemented("Hdf5 only implemented for KZ[46]==1")

        #self.set_h5_filename(snapshotfiles[0])
        #self.set_next_h5_filename(next_snapshotfile)
        self.step_dt = 2**(-self.inputpar["KZ"][47]) 
        self.requested_time = time
        self.snapshotfiles = snapshotfiles
        self.current_time = 0
        self.snapshot_id = 0 
        self.step_id = 0
        self.h5part_file = None

        self.dataset_list = ['KW', 'L', 'M', 'MC', 'NAM', 'POT', 'RC', 'RS', 'TE', 'V1', 'V2', 'V3', 'X1', 'X2', 'X3']
        self.data = None
        self.data_next = None
        self.data_interp = None

        self.initialize_code()

    def set_h5_filename(self, h5fn,snap_pos):
        if os.path.isfile(h5fn):
            self.h5part_file = h5py.File(h5fn, 'r')
        else:
            raise IOError('File %s not exist!' % h5fn)

    def get_step_dt(self):
        if self.step_dt == 0:
            raise ValueError("H5nb6xx input have step_dt == 0, check KZ(47)")
        return self.step_dt

    def get_number_of_steps(self):
        if self.h5part_file is not None:
            return len(self.h5part_file)
        else:
            return 0

    def get_adjacent_steps(self, time):
        return (int(numpy.floor(time / self.step_dt)), 
                int(numpy.ceil(time / self.step_dt)) )

    def initialize_snapshot_list(self):
        tstart,tend = numpy.array([]),numpy.array([])
        for snap in self.snapshotfiles: 
            h5part_file = h5py.File(snap,"r")
            step_l = len(h5part_file) - 1 
            tstart = numpy.append( tstart, 
                                   h5part_file["Step#%d" % 0].attrs["Time"] )
            tend = numpy.append( tend, 
                                 h5part_file["Step#%d" % step_l ].attrs["Time"])
            h5part_file.close()
        iord = numpy.argsort(tstart)
        self.tstart = tstart[iord]
        self.tend = tend[iord]
        self.snapshotfiles = numpy.array(self.snapshotfiles)[iord] 
        self.Nsnap = len(self.snapshotfiles)

    def set_current_snapshot(self,isnap=None) : 
        if isnap is None:
            try:
                isnap = numpy.argwhere( self.tstart <= self.requested_time )[-1][0] 
            except IndexError:
                raise IndexError("Requested time [%f] below first snapshot time [%f] "
                                 "contained in: %s" % (self.requested_time, 
                                                       self.tstart[0], 
                                                       self.snapshotfiles[0]
                                                                                ))
        self.snaphot_id = isnap
        self.close()

        self.h5part_file = h5py.File(self.snapshotfiles[isnap],"r") 

    def initialize_data(self):
        self.__emptyData = dict()

        # initialize empty ID vector and step vector that will update as new
        # stars may appear.
        self.id_vec = numpy.array([])
        #self.dt_vec = numpy.array([])
        self.step_vec = numpy.array([])
        self.step_next_vec = numpy.array([])

        for dset_name in self.dataset_list:
            self.__emptyData[dset_name] = numpy.array([])
        
        self.data = self.__emptyData.copy()
        self.data_next = self.__emptyData.copy()
        self.data_interp = self.__emptyData.copy()
        self.current_step_data = self.__emptyData.copy()
            
    def initialize_code(self):
        self.initialize_snapshot_list()
        tmin = self.tstart[0]
        tmax = self.tend[-1]
        if not tmin <= self.requested_time <= tmax :
            raise ValueError("Requested time [% f] outside snapshot boundaries"
                             " [% f, % f]"% (self.requested_time, tmin, tmax) )

        self.initialize_data()
        self.load_prev_step_data() #this advance the code until requested time
        self.load_next_step_data()

    def load_current_step_data(self):
        self.current_step_data,self.current_time = self.get_data_by_step_id(self.step_id)

    def get_data_by_step_id(self, stepid):
        data_dict = self.__emptyData.copy() 
        Step = self.h5part_file["Step#%d" % stepid ]
        for data_key in self.dataset_list :
            if data_key in Step :
                StepData = Step[data_key][()]
                data_dict[data_key] = StepData
        if "Binaries" in Step:
            Bstep = Step["Binaries"]
            bprim = dict()
            bsec = dict()
            Mcm = Bstep["M1"][()] + Bstep["M2"][()]
            for k in [1,2,3]:
                bprim["X%i" % k] = (Bstep["XC%i" % k][()] +
                                    Bstep["M2"][()]/Mcm*Bstep["XR%i" % k][()] )
                bsec["X%i"%k] = bprim["X%i"%k][()] + Bstep["XR%i"%k][()]
                bprim["V%i"%k] = (Bstep["VC%i"%k][()] 
                        + Bstep["M2"][()]/Mcm*Bstep["VR%i"%k][()] )
                bsec["V%i"%k] = bprim["V%i"%k][()] + Bstep["VR%i"%k][()]

            for k in ["KW","L","M","MC","NAM","RC","RS","TE" ]:
                bprim[k] = Bstep["%s1"%k ][()]
                bsec[k] = Bstep["%s2"%k ][()]

            #TODO: This must be fixed. Calculate the right potential for each member
            bprim["POT"] = Bstep["POT"] 
            bsec["POT"] = Bstep["POT"] 
            self.append_data(data_dict,bprim)
            self.append_data(data_dict,bsec)
        return data_dict,Step.attrs["Time"]

    def load_prev_step_data(self):
        """
        Load data scanning snapshots from t=t0 to t=requested_time.
        """
        snaps = range(self.Nsnap)
        goal_time = self.requested_time
        done = False
        for isnap in snaps :
            self.set_current_snapshot(isnap)
            filename = self.h5part_file.filename
            Nsteps = len(self.h5part_file) 
            step_ids = range(Nsteps)
            for step_id in step_ids:
                id_vec = self.data["NAM"]
                data_new,tstep = self.get_data_by_step_id(step_id)
                self.current_time = tstep
                print("Scanning file: %s - Time: %.3f          "%(
                    filename,tstep   ),end="\r" )
                ids_in_step = data_new["NAM"] 
                if len(id_vec) == 0:
                    self.data = data_new.copy()
                    self.step_vec = tstep*numpy.ones_like(ids_in_step,dtype=numpy.float32)
                    if tstep + self.step_dt >= goal_time :
                        done = True
                        break
                    continue

                #use with id_vec or data[key]
                ids_not_stored = numpy.isin(id_vec,ids_in_step,invert=True)
                not_stored_order = numpy.argsort(id_vec[ids_not_stored] )
                # use with ids_in_step or data_new[key]
                ids_to_update =  numpy.isin(ids_in_step,id_vec)
                to_update_order = numpy.argsort(ids_in_step[ids_to_update] )

                ids_to_add = numpy.invert(numpy.isin(ids_in_step,id_vec) )
                #to_add_order = numpy.argsort(ids_in_step[ids_to_add] )

                n1 = numpy.sum(ids_not_stored)
                n2 = numpy.sum(ids_to_update)
                
                assert  n1+n2 == len(id_vec)
                data = self.__emptyData.copy()
                #not stored first:
                self.append_data(data,self.data,ids_not_stored,not_stored_order)
                #stored, need update
                self.append_data(data,data_new,ids_to_update,to_update_order)
                #new stars
                self.append_data(data,data_new,ids_to_add)
                self.data = data

                new_step_vec = numpy.zeros_like(self.data["NAM"],dtype=numpy.float32)
                new_step_vec[0:n1]= self.step_vec[ids_not_stored][not_stored_order]
                new_step_vec[n1:] = tstep
                self.step_vec = new_step_vec

                if tstep + self.step_dt >= goal_time :
                    done = True
                    break
            if done :
                break

        id_vec = self.data["NAM"]
        order = numpy.argsort(id_vec)
        #self.data = self.__emptyData.copy()
        for key in self.data : 
            self.data[key] = self.data[key][order]
        self.step_vec = self.step_vec[order]
        self.id_vec = self.data["NAM"]

    def load_next_step_data(self):
        id_vec = self.id_vec
        snapid0 = self.snapshot_id
        stepid0 = self.step_id +1
        found_ids = numpy.array([])
        found_tstep = numpy.array([])
        found_data = self.__emptyData.copy() 

        done=False
        for snapid in range(snapid0,self.Nsnap):
            if snapid > snapid0 :
                h5file = h5py.File(self.snapshotfiles[snapid],"r")
                stepid0 = 0
            else:
                h5file = self.h5part_file

            nsteps = len(h5file) 
            for istep in range(stepid0,nsteps):
                idata,tstep = self.get_data_by_step_id(istep)
                #Step =  h5file["Step#%d"%istep]
                #tstep = Step.attrs["Time"]
                sids = idata["NAM"]
                print("Scanning file: %s - Time: %.3f          "%(
                    h5file.filename,tstep   ),end="" )
                already_found = numpy.isin(sids,found_ids) 
                unknown = numpy.isin(sids,id_vec,invert=True)
                # just found = NOT already_found AND NOT unknown
                just_found  =  numpy.invert(already_found) * numpy.invert(unknown)
                for group in [ just_found, unknown ]: 
                    order = numpy.argsort(sids[group])
                    self.append_data(found_data,idata,group,order)
                    found_tstep = numpy.concatenate([
                            found_tstep,
                            numpy.ones_like(sids[group])*tstep
                            ])
                
                found_ids = numpy.concatenate([found_ids,sids[just_found]
                    ] ).astype(int)
                print(" Found: %i/%i"%(
                    numpy.sum(numpy.isin(id_vec,found_ids,invert=False)),
                    len(id_vec)
                    ))
                if numpy.sum(numpy.isin(id_vec,found_ids,invert=True)) == 0:
                    done = True
                    break
            if done:
                if snapid != self.snapshot_id:
                    h5file.close()
                break
            #if snapid != self.snaphot_id:
            #    h5file.close()

        #finally make sure the firts N common ids are in the same order
        # id_vec is assumed to be already sorted 
        # TODO: Remove following assert once code is tested
        # Assert that the id_vec is sorted
        assert numpy.array_equal(numpy.argsort(id_vec),numpy.arange(len(id_vec)))
        found = numpy.isin(found_data["NAM"],id_vec) 
        new = numpy.invert(found)
        ifound = numpy.argsort(found_data["NAM"][found] ) 
        inew = numpy.argsort(found_data["NAM"][new] ) 
        
        data = self.__emptyData.copy()
        self.append_data( data,found_data,found,ifound )
        self.append_data( data,found_data,new,inew )
        # for key in found_data : 
            # found_data[key] = numpy.concatenate([
                # found_data[key][found][ifound],
                # found_data[key][new][inew]
            # ])
        found_tstep = numpy.concatenate([found_tstep[found][ifound],
                                                found_tstep[new][inew]
                                                ])
        self.data_next = found_data
        self.step_next_vec = found_tstep
        self.id_vec_next = found_data["NAM"]
        return 0

    def append_data(self,data,data_new,mask=None,order=None):
        mask = range(len(data_new["NAM"])) if mask is None else mask
        order = range(len(data_new["NAM"])) if order is None else order
        for key in self.dataset_list : 
            key_dtype = data_new[key].dtype
            if key in data_new: 
                data[key] = numpy.concatenate([
                                data[key],
                                data_new[key][mask][order]   
                                ]).astype(key_dtype)
        return

    def evolve_model(self,time):
        self.step_vec,self.data = self.scan_data(time,
                                  self.step_vec, self.data,1)
            
    def evolve_step(self):
        self.step_id += 1 
        if self.current_time + self.step_dt > self.tend[self.snapshot_id] :
            self.snapshot_id += 1 
            self.set_current_snapshot(self.snapshot_id)
            self.step_id = 0
        self.load_current_step_data()

    def close(self):
        if self.h5part_file is None:
            return
        else :
            self.h5part_file.close()
