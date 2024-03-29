#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Read __author__, __email__. from __init__.py
__author__ = "__author__"
__email__ = "__email__"
for line in open("src/pygenome/__init__.py"):
    if line.startswith("__author__") or line.startswith("__email__"):
        exec(line.strip())

from setuptools import setup
from setuptools import Command
from setuptools import find_packages

from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


class PyTest(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys

        errno = subprocess.call([sys.executable, "run_test.py"])
        raise SystemExit(errno)


setup(  name='pygenome',
        author          =__author__,
        author_email    =__email__,
        zip_safe = False,
        cmdclass={"test": PyTest},
        packages=find_packages("src"),
        include_package_data=True,
        package_dir={"": "src"},
        url='https://github.com/BjornFJohansson/pygenome',
        license='LICENSE.txt',
        package_data={'src/pygenome': ['Saccharomyces_cerevisiae.zip',
                                       'pFA6a-kanMX4.gb',
                                       'pFA6a-GFPS65T-kanMX6.gb' ]},


        description='''Accessing the Saccharomyces cerevisiae genome from Python''',
        long_description=long_description,
        long_description_content_type='text/markdown',
        setup_requires=["pytest-runner",
                        "setuptools_scm"],
        install_requires =[ "pydna",
                            "appdirs",
                            "tqdm",
                            "requests",
                            "prettytable" ],
        tests_require=["pytest"],
        use_scm_version={"write_to": "src/pygenome/_version.py"},
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       "Intended Audience :: Developers",
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 3.7',
                       'Programming Language :: Python :: 3.8',
                       'Programming Language :: Python :: 3.8',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',]
        )
