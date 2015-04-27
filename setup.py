#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = 'pygenome/_version.py'
versioneer.versionfile_build = 'pygenome/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = '' # dirname like 'myproject-1.2.0'

# Read version numbers, author etc..
__version__ = "Undefined"
for line in open('pygenome/__init__.py'):
    if line.startswith('__') and not line.startswith('__version__'):
        exec(line.strip())

# Change line ending to windows for all text files
import os
for root, dirs, files in os.walk(os.path.abspath(os.path.dirname(__file__))):
    for name in files:
        if not name.lower().endswith(".txt"):
            continue
        filename = os.path.join(root, name)
        with open(filename, "rb") as f:
            data = f.read()
        temp = data.replace('\r\n', '\n')
        temp = temp.replace('\r', '\n')
        temp = temp.replace('\n', '\r\n')
        if not data == temp:
            with open(filename, "wb") as f:
                f.write(temp)
                print "changed", filename

from setuptools import setup

setup(  name='pygenome',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        author          =__author__,
        author_email    =__email__,
        packages=['pygenome'],
        url='http://pypi.python.org/pypi/pygenome/',
        license='LICENSE.txt',
        description='''Accessing the Saccharomyces cerevisiae genome from Python''',
        long_description=open('README.rst').read(),
        install_requires =[ "biopython", "appdirs", "percache", "pydna"],
        test_suite = 'nose.collector',
        zip_safe = False,
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 2.7',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])
