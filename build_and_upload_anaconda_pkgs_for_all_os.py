#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "1.0.3"

print('''

This script automatically builds, converts and uploads conda packages.

It is meant to be placed in the same folder as the meta.yaml
''')
from termcolor import colored
import subprocess
import pathlib

bashCommand = "conda build . --output"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)

pkgpth   = pathlib.Path(process.communicate()[0].strip())
pkgname  = pkgpth.name
builddir = pkgpth.parent.parent

print(bashCommand)
print("Package  : ",colored( pkgpth, "blue"))
print("Build dir: ",colored( builddir, "blue"))

response = input("build (y/n) ? ")

if response.lower().startswith("y"):
    bashCommand = "conda build . --no-anaconda-upload"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)
    output = process.communicate()[0]
    print(colored(output, "blue"))
else:
    print(colored("build skipped.", "red")) 

if "/noarch/" in str(pkgpth): # no convert
    print(colored("no convert, since package is NOARCH", "red"))
    pkgdirs = ["noarch"]
else: # ask to convert
    response = input("convert (y/n) ? ")
    if response.lower().startswith("y"):
        bashCommand = "conda convert {} -p all -o {}".format(str(pkgpth), str(builddir))
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)
        output = process.communicate()[0]
        print(colored(output, "blue"))
        pkgdirs = ["win-64", "osx-64", "win-32", "linux-64", "linux-32"]
    else:
        print(colored("convert skipped.", "red"))
        pkgdirs = ["linux-64"]


pts=[]

print()
print("The following packages have been created:")
print()

for dir_ in pkgdirs:
    upload_pth = builddir.joinpath(dir_, pkgname)
    print(colored(upload_pth, "blue"))
    pts.append(upload_pth)
    
response = input("upload (y/n) ? ")

if response.lower().startswith("y"):
    for p in pts:
        bashCommand = "anaconda upload {}".format(p)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)
        output = process.communicate()[0]
        print(colored(output, "blue"))
else:
    print("upload skipped.")
    
process = subprocess.Popen("conda build purge".split(), stdout=subprocess.PIPE, universal_newlines=True) 
output = process.communicate()[0]
print(output)
print("done!")