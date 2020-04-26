import numpy
from nbody6tools import Utilities
from nbody6tools.snowballing import snowballing_method


class Methods():
    """
    This class does not work by its own. It defines the methods for the
    Particle class defined in Datamodel.__init__.py self is the Particle
    Object.
    """

    def half_mass_radius(self,direction="average"):
        return Utilities.get_mass_radius(self,direction=direction,fraction=0.5)

    def mass_radius(self,direction="average",fraction=0.5):
        return Utilities.get_mass_radius(self,fraction=fraction,
                                         direction=direction)

    def bound_indexes(self,verbose = False):
        """ 
        Return bound positions using the snowballing method. Only works if
        stars are in physical units
        """
        #epot = self.x*0.0 # background potential #TODO implement
        parts=numpy.array([self.x,self.y,self.z,self.vx,self.vy,self.vz,
                           self.mass,self.epot],dtype=numpy.float64) 
                           #must be float64
        rcores = numpy.array([ 2*self.half_mass_radius() ],dtype=numpy.float64)
        #should be use center of density. Maybe initialize this value at start
        clumps =  numpy.array(self.center,dtype=numpy.float64)  
        cl_flags=snowballing_method(parts,clumps,rcores,verbose,
                self.gravitational_constant) 
        return numpy.array(cl_flags[0],  dtype=bool )

    def bound_subset(self,verbose = False,gravitational_constant=1):
        """ Returns a bound subset of particles using the snowballing method"""
        return self[self.bound_indexes(verbose)]

    def velocity_dispersion(self,direction="average"):
        """ 
        return one dimensional velocity dispersion (average of the three
        projections)
        """
        return Utilities.get_velocity_dispersion(self,direction)

    def potential_energy(self):
        """
        Returns the total potential energy of the particles.
        Uses internal PHI calculated by Nbody6
        """
        return (self.pot*self.mass).sum()/2.


    def kinetic_energy(self):
        """
        Returns kinetic energy of particles
        """
        ke = 0.5*self.mass*(self.vx**2 + self.vy**2 +self.vz**2 )
        return ke.sum()

    def virial_ratio(self,G = 4.3020077853E-3,smoothing_length=0.01):
        """
        Returns virial ratio. It follows routine scale.F in Nbody6
        TODO: Implement other external potentials, e.g. Plummer
        """
        epot = (self.epot*self.mass).sum()
        vir = self.potential_energy() - 2*epot
        return -self.kinetic_energy()/vir

    def center_of_mass(self):
        mtot = self.mass.sum()
        x = (self.mass*self.x).sum()/mtot
        y = (self.mass*self.y).sum()/mtot
        z = (self.mass*self.z).sum()/mtot
        return x,y,z

    def center_of_mass_velocity(self):
        mtot = self.mass.sum()
        vx = (self.mass*self.vx).sum()/mtot
        vy = (self.mass*self.vy).sum()/mtot
        vz = (self.mass*self.vz).sum()/mtot
        return vx,vy,vz

        
