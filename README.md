# nbody6tools
A set of tools for Nbody6 simulations.

## Install

After cloning, from inside clonning folder just do:
```
pip install ./
```
This package only works with python3


## Quick Tutorial:

Once installed you only need to do:

`from nbody6tools import Reader`

The Reader module contains functions that work with the output data from Nbody6. The most useful function for now is the `read_snapshot` function.

So if you do:
```
snapshot = Reader.read_snapshot("relative_path_to_simulation_folder/", snapshot=0 )
```
This returns a Snapshot object were you can see the properties in the `Datamodel` module. 
You can also do:

```
snapshot = Reader.read_snapshot("relative_path_to_folder/", time=10 )
```
where time is in Myr.

If you check the `Datamodel.__init__` you can see how snapshot works, but for now the most important thing is that:

```
snapshot.to_physical()
```
Transforms the internal values from Nbody units to astrophysical units: km/s, parsecs, Msun, Msun km²/s² for energy.
You can come back to nbody if you like with the `.to_nbody()` function.

Then the particles are obtained with:
`stars = snapshot.stars`

`stars` is a Particles object defined in `Datamodel`. It have several useful functions, but basically give you access to the stellar quantities like for instance:
```
x = stars.x  #position in x
vx = stars.vx #velocity on x
r = stars.r  # distance to the center
v = stars.v  # velocity module
```
All of these are `numpy.arrays`

Check for other properties in `Datamodel.__init__` it also contains some functions like:
```
rh = stars.half_mass_radius() #
cmx,cmy,cmz = stars.center_of_mass() #
cmvx,cmvy,cmvz = stars.center_of_mass_velocity() #
vd = stars.velocity_dispersion() # 
```
For other functions check the Particles object definition, but for extra available functions check the file `Datamodel/_ParticleMethods` 


## Available set of stars:
Depending if you have binaries you may want them resolved or bound or both. 
Here a list of available subset of particles:

#### Resolved sets
Each star in these sets represents ONE star:

`snapshot.stars`  :  All stars in the set

`snapshot.bound_stars`  : Stars bound 

`snapshot.unbound_stars` : Stars unbound

#### Unresolved sets:
Some stars represent a binary system. These sets are usefull to calculate global quantities:

`snapshot.unresolved_stars` : All stars in the set unresolved

`snapshot.bound_stars_unresolved` : Stars bound unresolved

`snapshot.unbound_stars_unresolved` : Stars unbound unresolved

##### Binary components:
You may want to distinguish binaries and singles from the above unresolved groups.
You can either get the full set together with resolved stars with:

`resolved_set = snapshot.resolve_set(snapshot.unresolved_stars)`

Or obtain three separate set of stars with:

`singles,primaries,secondaries = snapshot.resolve_set(snapshot.unresolved_stars,split_set=True)`

`primaries` and `secondaries` have data in the same order for easy manipulation.


## Available scripts

Outside a python script you also should have access to the standalone programs
`nbtools` and  `nbt-movie`.

`nbtools`  is meant to generate quick graphs to check out the simulation. 
`nbt-movie` is meant to quickly generate an animation

Do `nbtools -h` and `nbt-movie -h` for help.


# Disclaimer
I made this just for things I used for my PhD. So there are many options from Nbody6 that may not be implemented, and also options here that only work for my custom versions of Nbody6.
Please let me know if there is any feature you need and I'll implement it.


