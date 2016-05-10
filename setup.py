#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import os


datadir = os.path.join(os.path.realpath(__file__),'Data')
datafiles = [(d, [os.path.join(d,f) for f in files])
    for d, folders, files in os.walk(datadir)]


ext = Extension('_spectrum', sources=['cext/_spectrum.c', 'cext/spectrum.c', 
                                      'cext/get_absorption.c', 'cext/get_continuum.c',
                                      'cext/voigt.c', 'cext/util.c'])
 
setup(name='dudeutils', 
    version='2.0', 
    description='utilities for dude', 
    ext_modules=[ext],
    package_dir = {'','dudeutils'},
    data_files=datafiles,
    #test_suite='nose.collector',
    #tests_require=['nose'],



    )
