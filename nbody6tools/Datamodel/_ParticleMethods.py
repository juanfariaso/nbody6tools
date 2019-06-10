from nbody6tools.snowballing import snowballing_method
from nbody6tools import Utilities
import numpy


class Methods():
    """
    This class does not work by its own. It defines the methods for the Particle class
    defined in Datamodel.__init__.py
    """

    def half_mass_radius(self):
        return Utilities.get_mass_radius(self,fraction=0.5)

    def bound_indexes(self,verbose = False):
        """ Return bound positions using the snowballing method. Only works if stars are in physical units"""
        epot = self.x*0.0 # background potential #TODO implement
        parts=numpy.array([self.x,self.y,self.z,self.vx,self.vy,self.vz,self.mass,epot],dtype=numpy.float64) #must be float64
        rcores = numpy.array([ 2*self.half_mass_radius() ],dtype=numpy.float64 )
        clumps =  numpy.array(self.center,dtype=numpy.float64)  #should be use center of density. Maybe initialize this value at start
        cl_flags=snowballing_method(parts,clumps,rcores,verbose) 
        return cl_flags

    def bound_subset(self):
        """ Returns a bound subset of particles using the snowballing method"""
        return stars[stars.bound_indexes]



