import os
import atexit
import h5py
import numpy 
from nbody6tools import Reader
import time #debug
from multiprocessing import Process,Queue,Manager
import traceback

class Data(object):
    def __init__(self):
        self.dataset_list = ['KW', 'L', 'M', 'MC', 'NAM', 'POT', 'RC', 'RS',
            'TE', 'V1', 'V2', 'V3', 'X1', 'X2', 'X3','Time']
        self.data = dict()
        for dset_name in self.dataset_list:
            self.data[dset_name] = numpy.array([])
        self.times = numpy.array([])

    def load_from_step(self,Step):
        if "N_SINGLE" in Step.attrs :
            data_dict = dict()
            nsingle = Step.attrs["N_SINGLE"]
            for data_key in self.dataset_list :
                if data_key in Step :
                    StepData = Step[data_key][()]
                    data_dict[data_key] = StepData
            data_dict["Time"] = numpy.full(nsingle, Step.attrs["Time"] )
            self.data = data_dict
        if "Binaries" in Step:
            AUtoPC = 4.8481368111333442e-06 
            Bstep = Step["Binaries"]
            nbin = Bstep.attrs["N_BINARY"]

            bprim = dict()
            bsec = dict()
            Mcm = Bstep["M1"][()] + Bstep["M2"][()]
            for k in [1,2,3]:
                bprim["X%i" % k] = (
                    Bstep["XC%i" % k][()] +
                    Bstep["M2"][()]/Mcm*Bstep["XR%i" % k][()]*AUtoPC
                    )
                bsec["X%i"%k] = bprim["X%i"%k][()] - Bstep["XR%i"%k][()]*AUtoPC
                bprim["V%i"%k] = (
                    Bstep["VC%i"%k][()] 
                    + Bstep["M2"][()]/Mcm*Bstep["VR%i"%k][()]
                    )
                bsec["V%i"%k] = bprim["V%i"%k][()] - Bstep["VR%i"%k][()]

            for k in ["KW","L","M","MC","NAM","RC","RS","TE" ]:
                bprim[k] = Bstep["%s1"%k ][()]
                bsec[k] = Bstep["%s2"%k ][()]

            #TODO: This must be fixed. Calculate the right potential for each member
            bprim["POT"] = Bstep["POT"] 
            bsec["POT"] = Bstep["POT"] 
            
            bprim["Time"] = numpy.full(nbin, Step.attrs["Time"] )
            bsec["Time"] = numpy.full(nbin, Step.attrs["Time"] )
            if ( len(bprim["NAM"]) == len(set(bprim["NAM"]))  or 
                len(bsec["NAM"]) == len(set(bsec["NAM"])) ) :
                self.append(bprim)
                self.append(bsec)
            else:
                print("WARNING: skipping binares on step %f due to duplicated binary"  
                  " member."%Step.attrs["Time"])

        if "Mergers" in Step: 
    # TODO: Implement mergers in BTSreader
            raise NotImplemented("Megers not implemented yet")

        self.sort()
        #self.check_duplicated(data_dict["NAM"],"combined")
        #self.data = data_dict

    def append(self,data_new,mask=None,order=None):
        mask = numpy.ones(len(data_new["NAM"]),dtype=bool) if mask is None else mask
        order = range(mask.sum()) if order is None else order

        if numpy.invert(mask).all() : #if all are False
           return
        if mask.sum() != len(order):
            raise IndexError("mask must contain same number of elements than order")

        for key in self.dataset_list : 
            key_dtype = data_new[key].dtype
            #print(key,data_new[key][:5])
            if key in data_new: 
                self.data[key] = numpy.concatenate([
                                self.data[key],
                                data_new[key][mask][order]   
                                ]).astype(key_dtype)
        return

    def update(self,data_new):
        data0 = self.data.copy()
        id_vec = data0["NAM"]
        new_ids = data_new["NAM"] 
        if len(id_vec) == 0:
            self.append(data_new)
            return
        #use with id_vec or data[key]
        ids_not_stored = numpy.isin(id_vec,new_ids,invert=True)
        not_stored_order = numpy.argsort(id_vec[ids_not_stored] )
        # use with ids_in_step or data_new[key]
        ids_to_update =  numpy.isin(new_ids,id_vec)
        to_update_order = numpy.argsort(new_ids[ids_to_update] )

        ids_to_add = numpy.invert(numpy.isin(new_ids,id_vec) )
        #to_add_order = numpy.argsort(ids_in_step[ids_to_add] )

        n1 = numpy.sum(ids_not_stored)
        n2 = numpy.sum(ids_to_update)
        #n3 = numpy.sum(ids_to_add)
        
        assert  n1+n2 == len(id_vec)
        self.clear()
        #not updated first:
        self.append(data0,ids_not_stored,not_stored_order)
        #to update
        self.append(data_new,ids_to_update,to_update_order)
        #new stars
        self.append(data_new,ids_to_add)
        #keep common data sorted
        order = numpy.concatenate([
            numpy.argsort( self.data["NAM"][:n1+n2] ),
            n1+n2 + numpy.argsort( self.data["NAM"][n1+n2:] ) 
            ])
        self.sort(order)

    def __len__(self):
        return len( self.data["NAM"]  )

    def sort(self,order=None):
        """ sort the data. If no ordered indexes privided will sort by NAM"""
        if order is None:
            order = numpy.argsort(self.data["NAM"])
        for key in self.dataset_list: 
            if key in self.data: 
                self.data[key] = self.data[key][order]
        return

    def clear(self):
        for key in self.dataset_list : 
            key_dtype = self.data[key].dtype
            if key in self.data: 
                self.data[key] = numpy.array([]).astype(key_dtype)

    def copy(self):
        data = Data()
        data_copy = self.data.copy()
        data.append(data_copy)
        return data

    def __getitem__(self,key) :
        return self.data[key]

    def __setitem__(self,key,value) :
        assert len(value) == self.__len__()
        self.data[key] = value

class H5nb6xxSnapshot(object):
    """ 
    Implement a reader for the Block Timestep Storage introduced by Cai et al.
    (2015) http://dx.doi.org/10.1088/0067-0049/219/2/31 

    Only KZ(46) = 1 with KZ(19)!=0 and KZ(12) != -1 is currently implemented.
    This mean that output is in physical units, but time is in Nbody.


    Parameters
    ----------
    inputfile : string
                Path to file with input parameters.
    snapshotfiles : list
                  List of filenames with the snapshots. This due that is not
                  guaranteed that all stars are stored in the same file or
                  adjacent files.  So, to be sure, the full list should be
                  provided.
    time : float
           Requested time.  Default = 0 

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
    def __init__(self,snapshotfiles,inputfile,time=0,Njobs=1,buffer_length=1):
        self.inputpar = Reader.parse_inputfile(inputfile)
        if self.inputpar["KZ"][46] != 1:
            raise NotImplemented("Hdf5 only implemented for KZ[46]==1")

        self.step_dt = 2**(-self.inputpar["KZ"][47]) 
        self.requested_time = time
        self.snapshotfiles = snapshotfiles
        self.current_time = 0
        self.snapshot_id = 0 
        self.step_id = 0
        self.h5part_file = None
        self.skipped_steps = 0 #debug

        self.dataset_list = ['KW', 'L', 'M', 'MC', 'NAM', 'POT', 'RC', 'RS', 'TE', 'V1', 'V2', 'V3', 'X1', 'X2', 'X3']
        self.data = None
        self.data_next = None
        self.data_interp = None
        self.tappend = 0 #debug
        self.Njobs = Njobs
        self.Nbf = buffer_length
        self.closed=False

        #atexit.register(self.close)
        self.initialize_code()

    # def get_step_dt(self):
        # if self.step_dt == 0:
            # raise ValueError("H5nb6xx input have step_dt == 0, check KZ(47)")
        # return self.step_dt

    # def get_number_of_steps(self):
        # if self.h5part_file is not None:
            # return len(self.h5part_file)
        # else:
            # return 0

    def __initialize_snapshot_list(self):
        tstart,tend,Nsteps = numpy.array([]),numpy.array([]), numpy.array([],dtype=int)
        for snap in self.snapshotfiles: 
            h5part_file = h5py.File(snap,"r")
            step_l = len(h5part_file) - 1 
            tstart = numpy.append( tstart, 
                                   h5part_file["Step#%d" % 0].attrs["Time"] )
            tend = numpy.append( tend, 
                                 h5part_file["Step#%d" % step_l ].attrs["Time"])
            Nsteps = numpy.append( Nsteps, 
                                  len(h5part_file))

            h5part_file.close()
        iord = numpy.argsort(tstart)
        self.tstart = tstart[iord]
        self.tend = tend[iord]
        self.Nsteps = Nsteps[iord]
        self.snapshotfiles = numpy.array(self.snapshotfiles)[iord] 
        self.Nsnap = len(self.snapshotfiles)

    def _set_current_snapshot(self,isnap=None) : 
        if isnap is None:
            try:
                isnap = numpy.argwhere( self.tstart <= self.requested_time )[-1][0] 
            except IndexError:
                raise IndexError("Requested time [%f] below first snapshot time [%f] "
                                 "contained in: %s" % (self.requested_time, 
                                                       self.tstart[0], 
                                                       self.snapshotfiles[0]
                                                                                ))
        self.snapshot_id = isnap
        #self.close()
        self.h5part_file = self.open_snapshotfile(isnap)

    def open_snapshotfile(self,snapid):
        if self.snapshotfiles[snapid] in self.__OpenFiles:
            return self.__OpenFiles[self.snapshotfiles[snapid]]

        if os.path.isfile(self.snapshotfiles[snapid]):
            h5file = h5py.File(open(self.snapshotfiles[snapid],"rb"),"r")
            self.__OpenFiles[h5file] = h5file
            return h5file
        else:
            raise IOError('File %s do not exist!' % h5fn)

    def __initialize_data(self):
        self.__OpenFiles = dict()

        self.data = Data()#self.__emptyData.copy() 
        #next step data: FULL list including new particles (appended after known in self.data ) 
        self.data_next = Data() #self.__emptyData.copy() 
        # syncronized data: FULL list (only stars in both data and data_next)
        self.data_interp = Data() #self.__emptyData.copy() 
        #active partilces. Particles stored in the current step
        self.current_step_data = Data()#self.__emptyData.copy() 
            
    def initialize_code(self):
        self.__initialize_snapshot_list()
        tmin = self.tstart[0]
        tmax = self.tend[-1]
        if not tmin <= self.requested_time <= tmax :
            raise ValueError("Requested time [% f] outside snapshot boundaries"
                             " [% f, % f]"% (self.requested_time, tmin, tmax) )

        self.__initialize_data()
        self._BufferDaemon = BufferDaemon(self.snapshotfiles,
                                           self.Nsteps,self.Nbf)
        self._BufferDaemon.start()
        self.load_prev_step_data() #this advance the code until requested time
        self.load_current_step_data() #load current step data and get data_next
        self.interpolate() #

    def load_current_step_data(self):
        #t0 = time.time() #debug
        tstep = self.h5part_file["Step#%d"%self.step_id].attrs["Time"]
        key = "Snapshot%d:Step%d"%(self.snapshot_id,self.step_id)
        tstep,self.current_step_data, data_next = self._BufferDaemon.get_data(key)
        self.current_time = tstep
        #t1 = time.time() #debug
        self.data.update( self.current_step_data.data)
        self.data_next.update(data_next.data)
        #t2 = time.time() #debug
        assert ( (self.data_next["NAM"][:len(self.data)] -
                 self.data["NAM"]).sum() == 0 )
        self.check_duplicated(self.data["NAM"])
        self.check_duplicated(self.data_next["NAM"])
        #t3 = time.time() #debug

        #tt= time.time() - t0
        #print("t3,t2,t1")
        #print((t3-t2)/tt,(t2-t1)/tt,(t1-t0)/tt,tt)

        #self.delta_t  = self.step_next_vec[0:len(self.id_vec)] - self.step_vec
        self.delta_t = (  self.data_next["Time"][:len(self.data)]
                        - self.data["Time"] )

    def check_duplicated(self,ids,label=None):
        label = "" if label is None else label
        if len(ids) > len(set(ids)):
            raise AssertionError("Duplicated in ids: %s"%label)

    def load_prev_step_data(self):
        """
        Load data scanning snapshots from t=t0 to t=requested_time.
        """
        snaps = range(self.Nsnap)
        goal_time = self.requested_time
        done = False
        for isnap in snaps :
            self._set_current_snapshot(isnap)
            filename = self.h5part_file.filename.split("name=")[-1]
            Nsteps = self.Nsteps[isnap]
            step_ids = range(Nsteps)
            for step_id in step_ids:
                Step = self.h5part_file["Step#%d"%step_id]
                tstep = Step.attrs["Time"]
                data_new = Data()
                data_new.load_from_step(Step)
                #step_vec_new = numpy.ones(len(data_new["NAM"]) )*tstep
                #self.step_vec = self.update_data(self.data,self.step_vec,data_new,step_vec_new)
                self.data.update(data_new.data)
                print("Scanning file: %s - Time: %.3f          "%(
                    filename, tstep  ),end="\r" )

                self.current_time = tstep
                self.step_id = step_id
                self.snapshot_id = isnap

                if tstep + self.step_dt >= goal_time :
                    done = True
                    break
            if done :
                break

    def interpolate(self, to_time=None):
# TODO: Add warning and option in BTS interpolation for when accel and derivatives are found
        to_time = self.current_time if to_time is None else to_time
        dt = to_time - self.data["Time"]  #self.step_vec
        tau = dt / self.delta_t
        dt[dt<0] = 0

        dataset_interp = ['X1', 'X2', 'X3']
        vel = ['V1', 'V2', 'V3']
        self.data_interp = self.data.copy()
        self.data_interp["Time"] = numpy.full(len(self.data_interp),to_time)
        #acc = ['AX', 'AY', 'AZ']
        #adot = ['JX', 'JY', 'JZ']
        for dset_id, dset_name in enumerate(dataset_interp):
            X0 = self.data[dset_name][:len(self.data)]
            V0 = self.data[vel[dset_id]][:len(self.data)]
            #A0 = self.data[acc[dset_id]]
            #J0 = self.data[adot[dset_id]]
            X1 = self.data_next[dset_name][:len(self.data)]
            #V1 = self.data_next[vel[dset_id]]
            #A1 = self.data_next[acc[dset_id]]
            #J1 = self.data_next[adot[dset_id]]

            p0 = X0
            p1 = V0 * dt
            #jpf: I dont have A0 .. 
            p2 = 0
            p3 = 0
            #p2 = 0.5 * A0 * dt * dt
            #p3 = 1. / 6 * A0 * dt * dt * dt
            if X0.shape[0] == X1.shape[0]:
                #p4 = -1.0 / 6 * (4 * J0 + J1) * dt * dt * dt - 2.5 * (2 * A0 - A1) * dt * dt - 5 * (4 * V0 + 3 * V1) * dt - 35 * (X0 - X1)
                #p5 = 0.5 * (2 * J0 + J1) * dt * dt * dt + (10 * A0 - 7 * A1) * dt * dt + 3 * (15 * V0 + 13 * V1) * dt + 84 * (X0 - X1)
                #p6 = -1.0 / 6 * (4 * J0 + 3 * J1) * dt * dt * dt - 0.5 * (15 * A0 - 13 * A1) * dt * dt - 2 * (18 * V0 + 17 * V1) * dt - 70 * (X0 - X1)
                #p7 = 1.0 / 6 * (J0 + J1) * dt * dt * dt + 2 * (A0 - A1) * dt * dt + 10 * (V0 + V1) * dt + 20 * (X0 - X1)
                p4,p5,p6,p7 = 0, 0, 0, 0
                pred = p0 + p1 * tau + p2 * pow(tau, 2.0) + p3 * pow(tau, 3.0) + p4 * pow(tau, 4.0) + p5 * pow(tau, 5.0) + p6 * pow(tau, 6.0) + p7 * pow(tau, 7.0)
            else:
                pred = p0 + p1 * tau + p2 * pow(tau, 2.0) + p3 * pow(tau, 3.0)
            self.data_interp[dset_name] = pred

        # for key in self.dataset_list : 
            # if key in self.data and key not in dataset_interp:
                # self.data_interp[key] = self.data[key]
        

    #def evolve_model(self,time):
    #    self.step_vec,self.data = self.scan_data(time,
    #                              self.step_vec, self.data,1)
            
    def evolve_step(self):
        self.step_id += 1 
        if self.current_time + self.step_dt > self.tend[self.snapshot_id] :
            self.snapshot_id += 1 
            self._set_current_snapshot(self.snapshot_id)
            self.step_id = 0
        self.load_current_step_data()
        #print("evolving",len(self.current_step_data["NAM"]))
    
    def synchronize_particles(self):
        self.interpolate()

    def close(self):
        if not self.closed:
            self._BufferDaemon.close()
            if not self.h5part_file is None :
                self.h5part_file.close()
            self.closed = True

#TODO: BufferDaemon need to be closed Explicitly!. Should be automatic
class BufferDaemon(object):
    def __init__(self,snapshotfiles,Nsteps,Njobs=1, Nbf=1):
        #Daemon must be self contained. Should send out just the data buffered
        #Will assume we look forward from whatever snapshotfiles are given
        #snapshot files is assumed to be sorted
        self.sid = 0
        self.snapid = 0
        self.snapshotfiles = snapshotfiles
        self.manager = Manager()
        self.databuffer = self.manager.dict()
        self.JobsStatus = self.manager.dict()
        self.runningJobs = self.manager.Value("i",0)
        self.Master = Process( target = self.main )
        if type(Nsteps) is int:
            self.Nsteps = numpy.full_like(self.snapshotfiles,Nsteps,dtype="int")
        else:
            self.Nsteps = Nsteps
        self.Nbf = Nbf
        self.Njobs = Njobs
        self.Queue = Queue()
        self.JobExitQueue = Queue()
        self.closed = self.manager.Value("i",0)
        atexit.register(self.close)
        self.waitingtime = 0

    @property
    def RunningJobs(self):
        return self.runningJobs.get()

    def start(self):
        self.Master.start()

    def stop(self):
        self.Master.terminate()

    def worker(self,stepid,snapid):
        try : 
            key = "Snapshot%d:Step%d"%(snapid,stepid)
            tstep,data,data_next = next_step_finder(stepid,self.snapshotfiles[snapid:] )
            self.Queue.put( (key,tstep,data,data_next)  )
            self.JobExitQueue.put( (tstep,0,"") )
        except:
            pos = "Snapshot%d:Step%d"%(snapid,stepid)
            err = traceback.format_exc()
            self.JobExitQueue.put((pos,-1,err))
        return

    def collecter(self,Cdict,Queue):
        # get() block the code here until there is something to get
        # it expect a 3-length tuple wher first item is used as key
        while True:
            out = Queue.get() 
            Cdict[out[0]] = out[1:]

    def main(self):
        sid = self.sid
        snapid = self.snapid
        #This main process can not modify class attributes except if they are 
        # created with self.manager
        DataCollecter = Process(target = self.collecter, 
                            args = (self.databuffer,self.Queue),
                            daemon=True)
        StatusCollecter = Process(target = self.collecter, 
                            args = (self.JobsStatus,self.JobExitQueue),
                            daemon=True)
        DataCollecter.start()
        StatusCollecter.start()

        Jobs = []
        while True:
            if self.closed.get() == 1:
                DataCollecter.terminate()
                StatusCollecter.terminate()
                break
            if len(Jobs) < self.Njobs and len(self.databuffer) < self.Nbf:
                p =  Process(target = self.worker,
                             args = (sid,snapid),
                             daemon=True)
                p.start()
                Jobs.append(p) #do not work. It dont update 

                sid += 1
                if sid >= self.Nsteps[snapid]:
                    snapid +=1
                    sid = 0
            self.runningJobs.set(len(Jobs))
            for job in Jobs.copy():
                if not job.is_alive() :
                    #Do i need to do something else to free resources?
                    Jobs.remove(job)

            for key in self.JobsStatus.copy():
                status,err = self.JobsStatus.pop(key)
                if status < 0:
                    raise Exception("\n"+err)


    def get_data(self,key):
        while True:
            t0 = time.time()
            if (not key in self.databuffer) and self.RunningJobs == 0 and len(self.databuffer) >= self.Nbf :
                self.stop()
                raise KeyError("%s not in list"%key)
            if key in self.databuffer:
                result = self.databuffer.pop(key)
                return result
            self.waitingtime += time.time() - t0
    def close(self):
        self.closed.set(1) 
        time.sleep(1)
        self.stop()
        self.Queue.close()
        self.JobExitQueue.close()

def next_step_finder(stepid,snapshotfiles,verbose=False):
    """
    Helper for H5nb6xxSnapshot object.
    Read the current stored data in stepid and their next occurrence,
    collecting any new stars finded after.
    This function is intended to run isolated and being used with multiple 
    processes.

    Input
    -----
    snapshotfiles : sorted list of hdf5 files 
    stepid : stepid of first item of snapshot files.
    
    Output
    -----
    data0 : data at stepid on first snapshot
    data  : data containing next occurence of particles stored in data.
            data.data[ len(data0)::] contains new particles found on the way
    """
    
    found_ids = numpy.array([])
    found_data = Data() 
    data = Data()
    snapid0 = 0

    h5file = h5py.File(open(snapshotfiles[snapid0],"rb"),"r")
    data.load_from_step(h5file["Step#%d"%stepid ]) 
    id_vec = data.data["NAM"]
    tstep = h5file["Step#%d"%stepid ].attrs["Time"]

    stepid0 = stepid + 1
    done=False
    for snapid in range(snapid0,len(snapshotfiles)) :
        if h5file is not None:
            h5file.close()
        #print("opening",snapshotfiles[snapid])
        h5file = h5py.File(open(snapshotfiles[snapid],"rb"),"r")
        nsteps = len(h5file)

        for istep in range(stepid0,nsteps):
            iStep = h5file["Step#%d" % istep ]
            newHere = numpy.isin(get_stored_names(iStep),found_ids,invert=True)
            if newHere.sum() == 0:
                continue
            idata = Data()
            idata.load_from_step(iStep)
            sids = idata.data["NAM"]
            just_found = numpy.isin(sids,found_ids,invert=True)

            for group in [ just_found]: 
                order = numpy.argsort(sids[group])
                found_data.append(idata.data,group,order)
                # found_tstep = numpy.concatenate([
                        # found_tstep,
                        # numpy.full_like(sids[group],tstep)
                        # ])
            found_ids = found_data.data["NAM"]

            if verbose:
                print(" Found: %i/%i        "%(
                        numpy.sum(numpy.isin(id_vec,found_ids,invert=False)),
                        len(id_vec)
                    ),
                    end = "\r")
            if numpy.sum(numpy.isin(id_vec,found_ids,invert=True)) == 0:
                done = True
                break
        if done:
            #print("closing file",h5file)
            h5file.close()
            break

    #finally make sure the firts N common ids are in the same order
    # id_vec is assumed to be already sorted 
# TODO: Add debug mode to avoid unnessesary sorting
    # Assert that the id_vec is sorted
    assert numpy.array_equal(numpy.argsort(id_vec),numpy.arange(len(id_vec)))
    found = numpy.isin(found_data.data["NAM"],id_vec) 
    #self.check_duplicated(found_data["NAM"],"(B) duplicated in found data")
    new = numpy.invert(found)
    ifound = numpy.argsort(found_data.data["NAM"][found] ) 
    inew = numpy.argsort(found_data.data["NAM"][new] ) 

    dataout = Data()
    dataout.append( found_data.data,found,ifound )
    dataout.append( found_data.data,new,inew )
    # found_tstep = numpy.concatenate([found_tstep[found][ifound],
                                            # found_tstep[new][inew]
    return tstep,data,dataout

def get_stored_names(Step):
    names = numpy.array([])
    if "N_SINGLE" in Step.attrs :
        names = numpy.concatenate([names, Step["NAM"] ])
    if "Binaries" in Step :
        names = numpy.concatenate([names,Step["Binaries"]["NAM1"],Step["Binaries"]["NAM2"] ])
    return names
