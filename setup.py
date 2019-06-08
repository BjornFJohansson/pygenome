#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer

# Read author etc. from __init__.py
for line in open('pygenome/__init__.py'):
    if line.startswith('__') and not line.startswith('__version__'):
        exec(line.strip()) 

from setuptools import setup
    
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(  name='pygenome',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        author          =__author__,
        author_email    =__email__,
        packages=['pygenome', "pygenome.tests"],
        package_data={'pygenome': [ 'Saccharomyces_cerevisiae.zip']},
        url='http://pypi.python.org/pypi/pygenome/',
        license='LICENSE.txt',
        description='''Accessing the Saccharomyces cerevisiae genome from Python''',
        long_description=long_description,
        install_requires =[ "pydna", "appdirs", "tqdm", "requests", "prettytable" ],
        zip_safe = False,
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 3.6',
                       'Programming Language :: Python :: 3.7',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',]
        )