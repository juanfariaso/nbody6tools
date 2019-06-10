import numpy

def get_mass_radius(stars,fraction=0.5,direction="average"):
    """ Calculate mass fraction radius (default: half mass) of stars  
        input : stars    : Dictionary containing mass and x,y,z
                fraction : Mass fraction to calculate. 
    """

    mcut = stars["mass"].sum()*fraction
    if direction == "average" :
        r = numpy.sqrt(stars["x"]**2 + stars["y"]**2 + stars["z"]**2)
    elif direction == "x":
        r = numpy.abs(stars["x"]) * numpy.sqrt(3.0)
    elif direction == "y":
        r = numpy.abs(stars["y"]) * numpy.sqrt(3.0)
    elif direction == "z":
        r = numpy.abs(stars["z"]) * numpy.sqrt(3.0)
    isort  = numpy.argsort(r)

    if mcut < stars["mass"][isort[0]]:
        return r[isort[0]]

    csum = numpy.cumsum(stars["mass"][isort])
    rh = r[isort][numpy.where(csum <= mcut)[0][-1]]
    return rh

def get_velocity_dispersion(stars,direction = "average"):
    """ Calculate velocity dispersion in a given direction 
        input: stars     : Dictionary containing vx,vy,vz
               direction : 'x','y' or 'z'. If not specified calculate all and take average.
    """
    if direction == "average":
        sx = numpy.std(stars["vx"])
        sy = numpy.std(stars["vy"])
        sz = numpy.std(stars["vz"])
        return (sx+sy+sz)/3.0
    elif direction == "x" : 
        return numpy.std(stars["vx"])
    elif direction == "y" : 
        return numpy.std(stars["vy"])
    elif direction == "z" : 
        return numpy.std(stars["vz"])

