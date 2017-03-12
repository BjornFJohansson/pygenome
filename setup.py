#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer

# Read author etc. from __init__.py
for line in open('pygenome/__init__.py'):
    if line.startswith('__') and not line.startswith('__version__'):
        exec(line.strip()) 

from setuptools import setup

with open("README.md") as f:
    long_description = f.read()
    
try:
    from pypandoc import convert_file
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    with open("README.md", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = convert_file("README.md", 'rst') 

setup(  name='pygenome',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        author          =__author__,
        author_email    =__email__,
        packages=['pygenome'],
        url='http://pypi.python.org/pypi/pygenome/',
        license='LICENSE.txt',
        description='''Accessing the Saccharomyces cerevisiae genome from Python''',
        long_description=long_description,
        install_requires =[ "pydna",        "percache" ],
        setup_requires =  ['pytest-runner', "percache"],
        tests_require  =  ['pytest'         "percache"],
        zip_safe = False,
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 3.5',
                       'Programming Language :: Python :: 3.6',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])

#versioneer.VCS = 'git'
#versioneer.versionfile_source = 'pygenome/_version.py'
#versioneer.versionfile_build = 'pygenome/_version.py'
#versioneer.tag_prefix = '' # tags are like 1.2.0
#versioneer.parentdir_prefix = '' # dirname like 'myproject-1.2.0'
#https://pypi.python.org/pypi?%3Aaction=list_classifiers
