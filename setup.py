#from skbuild import setup
#from setuptools import setup
import setuptools

setuptools.find_packages(),
with open("README.md", "r") as fh:
     long_description = fh.read()

#extensions = [ setuptools.Extension("Qparameter",["Qpar.f90"])]
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

#ext = Extension(name = 'nbody6tools.ext',
#                 extra_compile_args = ['-O3'],
#                 sources = ['nbody6tools/ext/src/Qpar.f90']
#                 )
snwbl = Extension(name = 'nbody6tools.snowballing',
                 sources = ['nbody6tools/ext/src/snowballing.f',
                            'nbody6tools/ext/src/snowballing.pyf'], # you may add several modules files under the same extension
                 extra_compile_args = ['-O3'],
                )

dependencies = ["scipy","h5py","matplotlib","tqdm",'numpy']

setup(
     name='nbody6tools',  
     version='0.1',
     scripts=['nbody6tools/scripts/nbtools','nbody6tools/scripts/nbt-movie'] ,
     author="Juan P. Farias",
     author_email="juanfariaso@gmail.com",
     description="Private tools for Nbody simulation research",
     long_description=long_description,
     long_description_content_type="text/markdown",
     packages=setuptools.find_packages(),
     #ext_modules = [ext,ext2],
     ext_modules = [snwbl],
     install_requires = dependencies,
 )

#setup(ext_modules = [flib] )
