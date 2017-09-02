#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer

# Read author etc. from __init__.py
for line in open('pygenome/__init__.py'):
    if line.startswith('__') and not line.startswith('__version__'):
        exec(line.strip()) 

from setuptools import setup
    
try:
    from pypandoc import convert_file
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    with open("README.md", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "\n"+convert_file("README.md", 'rst')

setup(  name='pygenome',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        author          =__author__,
        author_email    =__email__,
        packages=['pygenome'],
        package_data={'pygenome': [ 'Saccharomyces_cerevisiae.zip']},
        url='http://pypi.python.org/pypi/pygenome/',
        license='LICENSE.txt',
        description='''Accessing the Saccharomyces cerevisiae genome from Python''',
        long_description=long_description,
        #install_requires =[ "pydna",         "appdirs"],
        #setup_requires =  [ "pytest-runner", "appdirs"],
        #tests_require  =  [ "pytest",        "appdirs"],
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