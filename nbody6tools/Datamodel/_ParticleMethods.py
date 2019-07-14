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
        #epot = self.x*0.0 # background potential #TODO implement
        parts=numpy.array([self.x,self.y,self.z,self.vx,self.vy,self.vz,self.mass,self.epot],dtype=numpy.float64) #must be float64
        rcores = numpy.array([ 2*self.half_mass_radius() ],dtype=numpy.float64 )
        clumps =  numpy.array(self.center,dtype=numpy.float64)  #should be use center of density. Maybe initialize this value at start
        cl_flags=snowballing_method(parts,clumps,rcores,verbose) 
        return numpy.array(cl_flags[0],  dtype=bool )

    def bound_subset(self,verbose = False):
        """ Returns a bound subset of particles using the snowballing method"""
        return self[self.bound_indexes(verbose)]

    def velocity_dispersion(self,direction="average"):
        """ return one dimensional velocity dispersion (average of the three projections)"""
        return Utilities.get_velocity_dispersion(self,direction)


