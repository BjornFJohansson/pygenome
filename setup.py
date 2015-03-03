#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Read version numbers, author etc..
__version__ = "Undefined"
for line in open('pygenome/__init__.py'):
    if line.startswith('__'):
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
        version         =__version__,
        author          =__author__,
        author_email    =__email__,
        packages=['pygenome'],
        url='http://pypi.python.org/pypi/pygenome/',
        license='LICENSE.txt',
        description='''Accessing the Saccharomyces cerevisiae genome from Python''',
        long_description=open('README.rst').read(),
        install_requires =[ "biopython>=1.6.5", "appdirs>=1.4.0", "percache>=0.3.0",],
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
