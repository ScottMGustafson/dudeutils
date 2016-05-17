from setuptools import setup, find_packages
from setuptools.extension import Extension
import os


datadir = os.path.join(os.path.realpath(__file__),'data')
datafiles = [(d, [os.path.join(d,f) for f in files])
             for d, folders, files in os.walk(datadir)]


ext = Extension('_spectrum', 
                 sources=['dudeutils/cext/_spectrum.c', 
                          'dudeutils/cext/spectrum.c', 
                          'dudeutils/cext/get_absorption.c', 
                          'dudeutils/cext/get_continuum.c',
                          'dudeutils/cext/voigt.c', 
                          'dudeutils/cext/util.c'],
               )
 
setup(name='dudeutils', 
    version='2.0', 
    description='utilities for dude', 
    packages=find_packages(),
    package_dir={'dudeutils': 'dudeutils', 'dudeutils.tests':'tests'},
    ext_modules=[ext],

    #data_files=datafiles,
    #scripts=['scripts/run_random_sampling']
    #entry_points={
    #    'console_scripts': [
    #        'run_random_sampling=dudeutils.random_sampling:main',
    #    ],
    #},
    )

