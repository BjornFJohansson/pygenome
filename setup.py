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
    long_description = "\n"+convert_file("README.md", 'rst')


    
import os
from setuptools.command.install import install as _install

class install(_install):
    def run(self):
        _install.run(self)
        import appdirs 
        data_dir = appdirs.user_data_dir(os.path.join("pygenome", "Saccharomyces_cerevisiae"))
        try:
            os.makedirs( data_dir )
        except OSError:
            if os.path.isdir( data_dir ):
                import glob
                chfiles = glob.glob(os.path.join(data_dir, "chr*.gb"))
                if len(chfiles) == 16:
                    return
            else:
                raise            
        import shutil
        shutil.copy("Saccharomyces_cerevisiae.zip", data_dir)
        import zipfile
        with zipfile.ZipFile(os.path.join(data_dir, "Saccharomyces_cerevisiae.zip"), "r") as z:
            z.extractall( data_dir )

cmdclass={'install': install}
cmdclass.update( versioneer.get_cmdclass() )

setup(  name='pygenome',
        version=versioneer.get_version(),
        cmdclass=cmdclass,
        author          =__author__,
        author_email    =__email__,
        packages=['pygenome'],
        package_data={'zip': [ os.path.join('data','*') ]},
        url='http://pypi.python.org/pypi/pygenome/',
        license='LICENSE.txt',
        description='''Accessing the Saccharomyces cerevisiae genome from Python''',
        long_description=long_description,
        #install_requires =[ "pydna",         "percache", "appdirs"],
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