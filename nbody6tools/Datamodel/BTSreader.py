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
        
        self.data = self.__emptyData.copy() # previous step data (FULL)
        self.data_next = self.__emptyData.copy() #next step data (FULL)
        self.data_interp = self.__emptyData.copy() # syncronized data (FULL)
        self.current_step_data = self.__emptyData.copy() #active partilces only
            
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
        self.load_current_step_data() #active particles
        self.interpolate()

    def load_current_step_data(self):
        self.current_step_data,self.current_time = self.get_data_by_step_id(self.step_id)
        if (self.step_vec != self.current_time).sum() == 0 :
            self.current_step_data = self.data
            return
        nfound = len(self.current_step_data["NAM"])
        if nfound == 0 :
            return
        
        current_step_vec = numpy.ones(nfound)*self.current_time
        self.step_vec = self.update_data(self.data,self.step_vec,self.current_step_data,
                         current_step_vec)
        self.id_vec = self.data["NAM"]

        data_next_update,step_vec_next_update = self.find_next_step_data(self.current_step_data)

        self.step_next_vec = self.update_data(self.data_next,self.step_next_vec,
                             data_next_update,step_vec_next_update)

        self.delta_t  = self.step_next_vec[0:len(self.id_vec)] - self.step_vec

    def get_data_by_step_id(self, stepid,h5file=None):

        data_dict = self.__emptyData.copy() 
        h5part_file = self.h5part_file if h5file is None else h5file

        Step = h5part_file["Step#%d" % stepid ]
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
            try:
                self.check_duplicated(bprim["NAM"],"bprim")
                self.check_duplicated(bsec["NAM"],"bsec")
                self.check_duplicated(numpy.concatenate((bsec["NAM"],bprim["NAM"])),"bcom")
                self.append_data(data_dict,bprim)
                self.append_data(data_dict,bsec)
            except AssertionError:
                #TODO: BTS: handle duplicated binary, ignore only the one causing it.
                print("WARNING: skipping binares on step %f due to duplicated binary"  
                      " member."%Step.attrs["Time"])

        if "Mergers" in Step: 
            raise NotImplemented("Megers not implemented yet")

        
        self.sort_data(data_dict)
        self.check_duplicated(data_dict["NAM"],"combined")

        return data_dict,Step.attrs["Time"]

    def check_duplicated(self,ids,label=None):
        label = "" if label is None else label
        if len(ids) > len(set(ids)):
            raise AssertionError("Duplicated in ids: %s"%label)

    def update_data(self,data,step,data_new,step_new):
        data0 = data.copy()
        id_vec = data0["NAM"]
        new_ids = data_new["NAM"] 
        self.check_duplicated(new_ids)
        if len(id_vec) == 0:
            self.append_data(data,data_new)
            return step_new

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
        n3 = numpy.sum(ids_to_add)
        
        try:
            assert  n1+n2 == len(id_vec)
        except AssertionError:
            print(numpy.isin(id_vec,new_ids).sum(),len(id_vec),len(new_ids))
            print( numpy.isin(new_ids,id_vec ).sum() )
            print(len(set(new_ids)),"set")
            print(len(id_vec))
            print(len(new_ids))
            print(n1,n2,n3)
            print(self.data[ids_not_stored])
            print(self.data[ids_to_update])
        self.clear_data(data)
        #not stored first:
        self.append_data(data,data0,ids_not_stored,not_stored_order)
        #stored, need update
        self.append_data(data,data_new,ids_to_update,to_update_order)
        #new stars
        self.append_data(data,data_new,ids_to_add)

        step = numpy.concatenate([ step[ids_not_stored][not_stored_order],
                                   step_new[ids_to_update][to_update_order],
                                   step_new[ids_to_add]] ) 
        #keep common data sorted
        order = numpy.concatenate([
            numpy.argsort( data["NAM"][:n1+n2] ),
            n1+n2 + numpy.argsort( data["NAM"][n1+n2:] ) 
            ])

        self.sort_data(data,order = order  ) 
        step = step[order]

        return step


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
                data_new,tstep = self.get_data_by_step_id(step_id)
                step_vec_new = numpy.ones(len(data_new["NAM"]) )*tstep
                # print(self.data["NAM"])
                self.step_vec = self.update_data(self.data,self.step_vec,data_new,step_vec_new)
                # print(self.data["NAM"])
                # print
                print("Scanning file: %s - Time: %.3f          "%(
                    filename,tstep   ),end="\r" )
                # print("advanced")
                # id_vec = self.data["NAM"]
                # self.current_time = tstep
                # print("Scanning file: %s - Time: %.3f          "%(
                    # filename,tstep   ),end="\r" )
                # ids_in_step = data_new["NAM"] 
                # if len(id_vec) == 0:
                    # self.data = data_new.copy()
                    # self.step_vec = tstep*numpy.ones_like(ids_in_step,dtype=numpy.float32)
                    # if tstep + self.step_dt >= goal_time :
                        # done = True
                        # break
                    # continue

                # #use with id_vec or data[key]
                # ids_not_stored = numpy.isin(id_vec,ids_in_step,invert=True)
                # not_stored_order = numpy.argsort(id_vec[ids_not_stored] )
                # # use with ids_in_step or data_new[key]
                # ids_to_update =  numpy.isin(ids_in_step,id_vec)
                # to_update_order = numpy.argsort(ids_in_step[ids_to_update] )

                # ids_to_add = numpy.invert(numpy.isin(ids_in_step,id_vec) )
                # #to_add_order = numpy.argsort(ids_in_step[ids_to_add] )

                # n1 = numpy.sum(ids_not_stored)
                # n2 = numpy.sum(ids_to_update)
                
                # assert  n1+n2 == len(id_vec)
                # data = self.__emptyData.copy()
                # #not stored first:
                # self.append_data(data,self.data,ids_not_stored,not_stored_order)
                # #stored, need update
                # self.append_data(data,data_new,ids_to_update,to_update_order)
                # #new stars
                # self.append_data(data,data_new,ids_to_add)
                # self.data = data

                # new_step_vec = numpy.zeros_like(self.data["NAM"],dtype=numpy.float32)
                # new_step_vec[0:n1]= self.step_vec[ids_not_stored][not_stored_order]
                # new_step_vec[n1:] = tstep
                # self.step_vec = new_step_vec

                if tstep + self.step_dt >= goal_time :
                    done = True
                    break
            if done :
                break

        # id_vec = self.data["NAM"]
        # order = numpy.argsort(id_vec)
        # #self.data = self.__emptyData.copy()
        # for key in self.data : 
            # self.data[key] = self.data[key][order]
        # self.step_vec = self.step_vec[order]
        # self.id_vec = self.data["NAM"]

    def load_next_step_data(self):
        found_data,found_tstep = self.find_next_step_data(self.data)
        self.data_next = found_data
        self.step_next_vec = found_tstep
        self.id_vec_next = found_data["NAM"]
        self.id_vec = self.data["NAM"]
        self.delta_t  = self.step_next_vec[0:len(self.id_vec)] - self.step_vec

    def find_next_step_data(self,data):
        id_vec = data["NAM"]
        snapid0 = self.snapshot_id
        stepid0 = self.step_id + 1
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
                idata,tstep = self.get_data_by_step_id(istep,h5file)
                self.check_duplicated(idata["NAM"],"idata")
                #Step =  h5file["Step#%d"%istep]
                #tstep = Step.attrs["Time"]
                sids = idata["NAM"]
                print("Scanning file: %s - Time: %.3f          "%(
                    h5file.filename,tstep   ),end="" )
                already_found = numpy.isin(sids,found_ids) 
                unknown = numpy.isin(sids,id_vec,invert=True)
                # just found = NOT already_found AND NOT unknown
                #just_found  =  numpy.invert(already_found) * numpy.invert(unknown)
                just_found = numpy.isin(sids,found_ids,invert=True)

                #for group in [ just_found, unknown ]: 
                for group in [ just_found]: 
                    order = numpy.argsort(sids[group])
                    self.append_data(found_data,idata,group,order)
                    found_tstep = numpy.concatenate([
                            found_tstep,
                            numpy.ones_like(sids[group])*tstep
                            ])
                found_ids = found_data["NAM"]
                try:
                    self.check_duplicated(found_data["NAM"],"found_data")
                    self.check_duplicated(found_ids,"found ids")
                except AssertionError :
                    print(idata["NAM"][just_found]  )
                    print(idata["NAM"][unknown]  )
                    print(numpy.isin(idata["NAM"][just_found],
                          idata["NAM"][unknown]  ).sum() )
                    print(numpy.isin(idata["NAM"][unknown],
                          idata["NAM"][just_found]  ).sum() )
                    raise


                
                # found_ids = numpy.concatenate([found_ids,sids[just_found]
                    # ] ).astype(int)
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
        #print(id_vec)
        # Assert that the id_vec is sorted
        assert numpy.array_equal(numpy.argsort(id_vec),numpy.arange(len(id_vec)))
        found = numpy.isin(found_data["NAM"],id_vec) 
        self.check_duplicated(found_data["NAM"],"found data after loop")
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
        return data,found_tstep
        self.data_next = found_data
        self.step_next_vec = found_tstep
        self.id_vec_next = found_data["NAM"]
        self.delta_t  = self.step_next_vec[0:len(self.id_vec)] - self.step_vec
        return 0

    def append_data(self,data,data_new,mask=None,order=None):
        mask = numpy.ones(len(data_new["NAM"]),dtype=bool) if mask is None else mask
        order = range(mask.sum()) if order is None else order

        if numpy.invert(mask).all() : #if all are False
           return
        if mask.sum() != len(order):
            raise IndexError("mask must contain same number of elements than order")

        for key in self.dataset_list : 
            key_dtype = data_new[key].dtype
            if key in data_new: 
                data[key] = numpy.concatenate([
                                data[key],
                                data_new[key][mask][order]   
                                ]).astype(key_dtype)
        return

    def clear_data(self,data):
        for key in self.dataset_list : 
            key_dtype = data[key].dtype
            if key in data: 
                data[key] = numpy.array([]).astype(key_dtype)

    def sort_data(self,data,order=None):
        """ sort the data. If no ordered indexes privided will sort by NAM"""
        #order = range(len(data_new["NAM"])) if order is None else order
        if order is None:
            order = numpy.argsort(data["NAM"])
        for key in self.dataset_list : 
            if key in data: 
                data[key] = data[key][order]
        return



    def interpolate(self, to_time=None):
        to_time = self.current_time if to_time is None else to_time
        #if self.flag_interp is True:
        if True:
            # dt = to_time - self.dt_vec
            dt = to_time - self.step_vec
            tau = dt / self.delta_t
            #if dt < 0:
            #    dt = 0
            dt[dt<0] = 0

            dataset_interp = ['X1', 'X2', 'X3']
            vel = ['V1', 'V2', 'V3']
            #acc = ['AX', 'AY', 'AZ']
            #adot = ['JX', 'JY', 'JZ']
            for dset_id, dset_name in enumerate(dataset_interp):
                X0 = self.data[dset_name]
                V0 = self.data[vel[dset_id]]
                #A0 = self.data[acc[dset_id]]
                #J0 = self.data[adot[dset_id]]
                X1 = self.data_next[dset_name]
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

            for key in self.dataset_list : 
                if key in self.data and key not in dataset_interp:
                    self.data_interp[key] = self.data[key]
            # self.data_interp['Mass'] = self.data['Mass']
            # self.data_interp['V1'] = self.data['V1']
            # self.data_interp['V2'] = self.data['V2']
            # self.data_interp['V3'] = self.data['V3']
            # self.data_interp['A1'] = self.data['A1']
            # self.data_interp['A2'] = self.data['A2']
            # self.data_interp['A3'] = self.data['A3']
            # self.data_interp['J1'] = self.data['J1']
            # self.data_interp['J2'] = self.data['J2']
            # self.data_interp['J3'] = self.data['J3']
            # self.data_interp['NAM'] = self.data['NAM']
        else:
            # No interpolation, just return the data as loaded from the HDF5
            for dset_name in self.dataset_list:
                self.data_interp[dset_name] = self.data[dset_name]

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
        #self.interpolate()
        #X = self.data_interp["X1"].mean()
        #print("Xmean",X)
        #For performance, only interpoate when nessesary
    
    def synchronize_particles(self):
        self.interpolate()

    def close(self):
        if self.h5part_file is None:
            return
        else :
            self.h5part_file.close()