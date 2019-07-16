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

    def potential_energy(self, G = 4.3020077853E-3,smoothing_length=0.01):
        """
        Returns the total potential energy of the self in the self set.
        """

        if len(self) < 2:
            return 0

        mass = self.mass
        x_vector = self.x
        y_vector = self.y
        z_vector = self.z
        epot = self.epot

        sum_of_energies = 0.0

        for i in range(len(self) - 1):
            x = x_vector[i]
            y = y_vector[i]
            z = z_vector[i]
            dx = x - x_vector[i+1:]
            dy = y - y_vector[i+1:]
            dz = z - z_vector[i+1:]
            dr = numpy.sqrt( (dx * dx) + (dy * dy) + (dz * dz) + smoothing_length*smoothing_length)
            m_m = mass[i] * mass[i+1:]

            energy_of_this_particle = (m_m / dr).sum() + epot[i]*mass[i]
            sum_of_energies -= energy_of_this_particle
        return G * sum_of_energies

    def kinetic_energy(self):
        ke = 0.5*self.mass*(self.vx**2 + self.vy**2 +self.vz**2 )
        return ke.sum()

    def virial_ratio(self,G = 4.3020077853E-3,smoothing_length=0.01):
        return -self.kinetic_energy()/self.potential_energy(G=G,smoothing_length=smoothing_length) 

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

        


