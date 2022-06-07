
import os
import numpy
from scipy.interpolate import interp1d


class ClusterOrbitInterpolator():
    def __init__(self,folder,outputfile = None):
        self.folder = folder
        self.outputfile =  self.__setup_output_file(outputfile)
        self.__interpolators = dict()
        self.__setup_data()
        self.__output_times = None


    def __setup_output_file(self,outputfile = None):
        if outputfile is None:
            possible_output = ['output','output.txt','out.txt','out']
            for out in possible_output:
                if os.path.isfile("/".join([self.folder,out]) ) :
                    outputfile = "/".join([self.folder,out] )
                    break
                else:
                    continue
            if outputfile is None:
                raise Exception("output file not found, tryed: %s"%possible_output )
        return outputfile

    def __setup_data(self):
        lines = open(self.outputfile,'r').readlines()
        result = []
        for l in lines:
            if "CLUSTER ORBIT" in l:
                result.append( numpy.array(l.split()[7:],dtype=numpy.float32 ) )
        self.data  = numpy.array(result)
        self.columns= {  'RGX' : 1, 'RGY' : 2, 'RGZ' : 3,
                'VGX' : 4, 'VGY' : 5, 'VGZ' : 6,
                'JZ' : 7 , 'ET' : 8
                }
        self.__columns = self.columns
        for k in self.keys():
            self.__interpolators[k] = None

    def keys(self):
        return self.columns.keys()

    def __setup_interpolator(self,key):
        self.__interpolators[key] = interp1d(self.data[:,0],
                self.data[:,self.columns[key]] )


    @property
    def times(self):
        if self.__output_times is None:
            raise Exception('Must setup the output times first')
        else:
            return self.__output_times
        
    @times.setter
    def times(self,times):
        self.__output_times = times

    @property
    def RG(self):
        ''' return the cluster position at times in pc
            times must be contained within simulation time range
            
            > orbit = ClusterOrbitInterpolator(folder)
            > orbit.times = [ 1,2,3,4,5]
            > x,y,z =  orbit.RG


        '''
        return self.RG_at_time(self.__output_times)


    def RG_at_time(self,t):
        result = []
        for rg in ['RGX','RGY','RGZ']:
            if self.__interpolators[rg] is None:
                self.__setup_interpolator(rg)
            result.append( self.__interpolators[rg](t) )

        return numpy.array(result)*1000 #to pc
            
    @property
    def VG(self):
        ''' return the cluster velocity at times in km/s
            times must be contained within simulation time range
            
            > orbit = ClusterOrbitInterpolator(folder)
            > orbit.times = [ 1,2,3,4,5]
            > vx,vy,vz =  orbit.VG
        '''
        self.VG_at_time(self.__output_times)

    def VG_at_time(self,t):
        result = [ ]
        for vg in ['VGX','VGY','VGZ']:
            if self.__interpolators[vg] is None:
                self.__setup_interpolator(vg)
            result.append( self.__interpolators[vg](t) )

        return numpy.array(result) #km/s

    @property
    def JZ(self):
        #TODO: ClusterOrbitInterpolator.JZ fix and check units add documentation 
        if self.__interpolators['JZ'] is None:
           self.__setup_interpolator('JZ')
        result = self.__interpolators['JZ'](self.__output_times) 
        return result

    def JZ_at_time(self,t):
        return self.__interpolators['JZ'](t)

    @property
    def ET(self):
        #TODO: ClusterOrbitInterpolator.ET fix and check units add documentation 
        if self.__interpolators['ET'] is None:
           self.__setup_interpolator('ET')
        result = self.__interpolators['ET'](self.__output_times) 
        return result
