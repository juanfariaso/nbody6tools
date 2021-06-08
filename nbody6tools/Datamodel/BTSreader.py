import h5py
import numpy 
import os
from nbody6tools import Reader

#def get_scale_factors(inputpar):
class Converter(object):
    def __init__(self,inputfile):
        self.input = Reader.parse_inputfile(inputfile)
        self.initialize()

    def initialize(self):
        #from units.f  of the Nbody6 code

        #Define GM & PC in cgs units and #AU in pc (2008 IAU values).
        GM = 6.6743e-08*1.9884e+33
        PC = 3.0856776e+18
        self.RBAR = self.inputpar["RBAR"]
        self.ZMBAR = self.inputpar["ZMBAR"]
        N = inputpar["N"]

        #Form time scale in seconds and velocity scale in km/sec.
        self.TSTAR = numpy.sqrt(PC/GM)*PC

        #Convert time scale from units of seconds to million years.
        self.TSTAR = self.TSTAR/(3.15e+07*1.0e+06)
        if (inputpar["KZ"][22] != 10) :
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
        self.snapshot_id=0
        self.step_id = 0

        self.dataset_list = ['KW', 'L', 'M', 'MC', 'NAM', 'POT', 'RC', 'RS', 'TE', 'V1', 'V2', 'V3', 'X1', 'X2', 'X3']
        self.data = None
        self.data_next = None
        self.data_interp = None

        self.initialize_code()


    def set_h5_filename(self, h5fn,snap_pos):
        if os.path.isfile(h5fn):
            self.h5part_file[snap_pos] = h5py.File(h5fn, 'r')
        else:
            raise IOError('File %s not exist!' % h5fn)

    def get_step_dt(self):
        if self.step_dt == 0:
            raise ValueError("H5nb6xx input have step_dt == 0, check KZ(47)")
        return self.step_dt

    def get_number_of_steps(self,snap_pos=0):
        "sid : -1 prev, 0 current , 1 next snapshotfile"
        if self.h5part_file[snap_pos] is not None:
            return len(self.h5part_file[ snap_pos ])
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

        self.h5part_file[0] = h5py.File(self.snapshotfiles[isnap],"r") 
        if isnap > 0 :
            self.h5part_file[-1] = h5py.File(self.snapshotfiles[isnap-1],"r")
        else :
            self.h5part_file[-1] = None

        if isnap < (self.Nsnap - 1) :
            self.h5part_file[1] = h5py.File(self.snapshotfiles[isnap+1],"r")
        else:
            self.h5part_file[1] = None

    def initialize_data(self):
        self.data = dict()
        self.data_next = dict()
        self.data_interp = dict()

        # initialize empty ID vector and step vector that will update as new
        # stars may appear.
        self.id_vec = numpy.array([])
        #self.dt_vec = numpy.array([])
        self.step_vec = numpy.array([])
        self.step_next_vec = numpy.array([])

        for dset_name in self.dataset_list:
            self.data[dset_name] = numpy.array([])
            self.data_next[dset_name] = numpy.array([])

    def initialize_code(self):
        self.initialize_snapshot_list()
        #Nsteps = self.get_total_number_of_steps()
        tmin = self.tstart[0]
        tmax = self.tend[-1]
        if not tmin <= self.requested_time <= tmax :
            raise ValueError("Requested time [% f] outside snapshot boundaries"
                             " [% f, % f]"% (self.requested_time, tmin, tmax) )

        #self.h5part_file_prev,self.h5part_file_current,self.h5part_file_next = None,None,None
        # self.h5part_files = { -1 : self.h5part_file_prev, 
                              # 0 : self.h5part_file_current,
                              # 1 : self.h5part_file_next}
        self.h5part_file = { -1 : None, 
                               0 : None,
                               1 : None}
        self.initialize_data()
        self.set_current_snapshot(self.snapshot_id)
        self.load_current_step_data()
        #self.initialize_adjacent_data()

    def load_current_step_data(self):
        self.data,self.current_time = self.get_data_by_step_id(self.step_id)
        

    def get_data_by_step_id(self, stepid, snap_pos=0):
        data_dict = { n : numpy.array([]) for n in self.dataset_list } 
        Step = self.h5part_file[snap_pos]["Step#%d" % stepid ]
        for data_key in self.dataset_list :
            if data_key in Step :
                StepData = Step[data_key][()]
                data_dict[data_key] = StepData

        return data_dict,Step.attrs["Time"]

    #def scan_data(self,time,step_start,step_end,step_vec,data):
    def scan_data_old(self,time,step_vec,data,direction=1):
        #direction = numpy.sign(step_end - step_start)
        #print("scanning",time,step_start,step_end,step_vec)
        if direction > 0:
            snaps = range(self.Nsnap)
        else:
            snaps = range(self.Nsnap-1,0,-1)
        done = False
        for isnap in snaps :
            self.set_current_snapshot(isnap)
            filename = self.h5part_file[0].filename
            Nsteps = len(self.h5part_file[0]) 
            if direction > 0:
                step_ids = range(Nsteps)
            else:
                step_ids = range(Nsteps-1,0,-1)
            for step_id in step_ids:
                id_vec = data["NAM"]
                data_new,tstep = self.get_data_by_step_id(step_id)
                print("Scanning file: %s - Time: %.3f          "%(
                    filename,tstep   ),end="\r" )
                ids_in_step = data_new["NAM"] 
                if len(id_vec) == 0:
                    data = data_new.copy()
                    step_vec = tstep*numpy.ones_like(ids_in_step,dtype=numpy.float32)
                    continue
                #second try

                #use with id_vec or data[key]
                ids_not_stored = numpy.invert(numpy.isin(id_vec,ids_in_step))
                not_stored_order = numpy.argsort(id_vec[ids_not_stored] )
                # use with ids_in_step or data_new[key]
                ids_to_update =  numpy.isin(ids_in_step,id_vec)
                to_update_order = numpy.argsort(ids_in_step[ids_to_update] )

                ids_to_add = numpy.invert(numpy.isin(ids_in_step,id_vec) )
                to_add_order = numpy.argsort(ids_in_step[ids_to_add] )

                n1 = ids_not_stored.sum()
                n2 = (ids_to_update).sum()
                assert  n1+n2 == len(id_vec)
                for key in self.dataset_list :
                    #not stored first:
                    data[key][0:n1] = data[key][ids_not_stored][not_stored_order] 
                    #stored, need update
                    data[key][n1:n1+n2] = data_new[key][ids_to_update][to_update_order]
                    #new stars
                    data[key] = numpy.concatenate((data[key],
                                                   data_new[key][ids_to_add] )
                                                  )
                new_step_vec = numpy.zeros_like(data["NAM"],dtype=numpy.float32)
                new_step_vec[0:n1]= step_vec[ids_not_stored][not_stored_order]
                new_step_vec[n1:] = tstep
                step_vec = new_step_vec

                if direction > 0 :
                    if tstep + self.step_dt >= time :
                        done = True
                        break
                if direction < 0 :
                    if tstep - self.step_dt < time :
                        done = True
                        break
            if done :
                break

        #print(step_vec) #debug
        return step_vec, data

    def scan_data(self,time,step_vec,data,mode = 1,ids_to_find=None):
        """
        Scan data from t=t0 to t=time.
        Mode: 
            0 : set data (previous step) from t=t0 until t=time
            1 : set data_next (next step) from t = time until all NAM in 
                in data is found
        """
        #direction = numpy.sign(step_end - step_start)
        #print("scanning",time,step_start,step_end,step_vec)

        #TODO maybe remove mode and only work with id_to_find
        if ids_to_find is not None:
            mode = 1
        else:
            mode = 0

        if  mode == 0:
            snaps = range(self.Nsnap)
        else:
            first_snap =  numpy.where(self.tstart >= time)[0][0]
            snaps = range(first_snap,self.Nsnap)
        done = False
        for isnap in snaps :
            self.set_current_snapshot(isnap)
            filename = self.h5part_file[0].filename
            Nsteps = len(self.h5part_file[0]) 
            if mode == 0:
                step_ids = range(Nsteps)
            else:
                t0 = self.h5part_file[0]["Step#0"].attrs["Time"]
                if time > t0 :
                    ifirst = int(numpy.ceil(( time -  t0)/self.step_dt )) 
                else:
                    ifirst = 1
                step_ids = range(ifirst,Nsteps)
                print( self.h5part_file[0]["Step#%d"%ifirst].attrs["Time"],time)
                assert self.h5part_file[0]["Step#%d"%ifirst].attrs["Time"]>time 

            for step_id in step_ids:
                id_vec = data["NAM"]
                if mode == 1:
                    id_vec0 = id_vec.copy()

                data_new,tstep = self.get_data_by_step_id(step_id)
                print("Scanning file: %s - Time: %.3f          "%(
                    filename,tstep   ),end="\r" )
                ids_in_step = data_new["NAM"] 
                if len(id_vec) == 0:
                    data = data_new.copy()
                    step_vec = tstep*numpy.ones_like(ids_in_step,dtype=numpy.float32)
                    continue

                #use with id_vec or data[key]
                ids_not_stored = numpy.invert(numpy.isin(id_vec,ids_in_step))
                not_stored_order = numpy.argsort(id_vec[ids_not_stored] )
                # use with ids_in_step or data_new[key]
                ids_to_update =  numpy.isin(ids_in_step,id_vec)
                to_update_order = numpy.argsort(ids_in_step[ids_to_update] )

                ids_to_add = numpy.invert(numpy.isin(ids_in_step,id_vec) )
                to_add_order = numpy.argsort(ids_in_step[ids_to_add] )

                #to use with id_vec and data[key]
                #II have to substract ids_to_add
                ids_already_found = numpy.isin(id_vec,ids_in_step[ids_to_update])
                ids_already_found = numpy.isin( )
                found_order = numpy.argsort(id_vec[ids_already_found] )

                if mode == 1  :
                    print
                    print(id_vec[ids_already_found],sum(ids_already_found))
                    print(ids_in_step[ids_to_update],sum(ids_to_update))
                    print
                n1 = ids_not_stored.sum()
                n2 = (ids_to_update).sum()
                assert  n1+n2 == len(id_vec)
                for key in self.dataset_list :
                    #not stored first:
                    data[key][0:n1] = data[key][ids_not_stored][not_stored_order] 
                    #stored, need update
                    if mode == 0 :
                        data[key][n1:n1+n2] = data_new[key][ids_to_update][to_update_order]
                    else:
                        #do not update if already found in mode=1
                        data[key][n1:n1+n2] = data[key][ids_already_found][found_order]
                    #new stars
                    data[key] = numpy.concatenate((data[key],
                                                   data_new[key][ids_to_add] )
                                                  )
                new_step_vec = numpy.zeros_like(data["NAM"],dtype=numpy.float32)
                new_step_vec[0:n1]= step_vec[ids_not_stored][not_stored_order]
                new_step_vec[n1:] = tstep
                step_vec = new_step_vec

                if mode == 0 :
                    if tstep + self.step_dt >= time :
                        done = True
                        break
                else:
                    missing = numpy.invert(numpy.isin(ids_to_find,data["NAM"])).sum()
                    print
                    print(missing)
                    if (  missing == 0  ) :
                        done = True
                        break
            if done :
                break

        #print(step_vec) #debug
        return step_vec, data

    def initialize_adjacent_data(self):
        #check that time is between current snapshot limits
        #isnap = self.snap_id
        #if not self.tstart[isnap] < self.requested_time < self.tend[isnap]:
        #    raise ValueError("Requested time out of current snapshot bound. "
        #                     "This should not happen here.")

        #Nsteps = self.get_number_of_steps(0)
        #assert Nsteps > 0
        ## Scanning from beggining to requested_time
        ## This should collect all stars, and set step_prev
        #self.step_vec,self.data = self.scan_data(self.requested_time,0,Nsteps,
        #                          self.step_vec,self.data)
        self.step_vec,self.data = self.scan_data(self.requested_time,
                                  self.step_vec, self.data,1)
        
        # Now scan from end to time, to catch any new star, and set step_next
        # self.step_next_vec,self.data_next = self.scan_data(self.requested_time,Nsteps-1,0,
                                  # self.step_next_vec,self.data_next)

        # self.step_next_vec,self.data_next = self.scan_data(self.requested_time,
                                  # self.step_next_vec, self.data_next,-1)

        self.step_next_vec,self.data_next = self.scan_data(self.requested_time,
                                  self.step_next_vec, self.data_next,
                                  ids_to_find=self.data["NAM"])
        #Sort data_next to be in the same order than data, new in data_next at
        # the end 
        common, i, j = numpy.intersect1d(self.data_next["NAM"], self.data["NAM"],
                                         return_indices=True)
        not_in_data = numpy.where(
            numpy.isin(self.data_next["NAM"], self.data["NAM"]) == False
        )[0]
        for key in self.dataset_list :
            self.data[key] = self.data[key][j]
            data_next = self.data_next[key][i] 
            self.data_next[key] = numpy.concatenate((data_next,
                                                     self.data_next[key][not_in_data]
                                                     ))
        self.id_vec = self.data["NAM"]
        self.step_vec = self.step_vec[j] 
        step_next_vec = self.step_next_vec[i]
        self.step_next_vec = numpy.concatenate((step_next_vec,
                                               self.step_next_vec[not_in_data]))
        self.dt = self.step_next_vec[0:len(i)] - self.step_vec


    def evolve_model(self,time):
        self.step_vec,self.data = self.scan_data(time,
                                  self.step_vec, self.data,1)
            
    def evolve_step(self):
        self.step_id += 1 
        if self.current_time + self.step_dt > self.tend[self.snapshot_id] :
            self.snapshot_id +=1 
            self.set_current_snapshot(self.snapshot_id)
            self.step_id = 0
        self.load_current_step_data()

    # def load_data_by_step_id(self, stepid, data_dict, snap_pos=0):
        # """may not be used"""
        # #data_dict = { n : numpy.array([]) for n in self.dataset_list } 
        # Step = self.h5part_file[snap_pos]["Step#%d" % stepid ]
        # h5_id_vec = Step["NAM"] 
        # new_ids_mask = numpy.invert(numpy.isin(h5_id_vec,self.id_vec))
        # stored_ids = h5_id_vec[numpy.isin(h5_id_vec,self.id_vec)]
        # stored_ind = Utilities.index_in_array(h5_id_vec,self.id_vec)

        # for data_key in self.dataset_list :
            # if data_key in Step :
                # StepData = Step[data_key][()]
                # #update existing data first
                # data_dict[data_key][stored_ind] = StepData[stored_ind]
                # #extend for new particles
                # vec = data_dict[data_key] 
                # data_dict[data_key] = numpy.concatenate(
                    # [ vec, StepData[new_ids_mask] ]
                # )

        # return Step["NAM"],Step.attrs["Time"]

    def close(self):
        for i in [-1,0,1]:
            if self.h5part_file[i] is not None :
                self.h5part_file[i].close()

    # def get_adjacent_files(snap,time): 
    # def get_snapName(self,time,dt):
        # tnext = self.time + dt
        # return str(round(tnext,4))
