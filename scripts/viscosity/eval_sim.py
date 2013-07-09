#! /usr/bin/env python2.3
#
# This file is part of the SYMPLER package.
# https://github.com/kauzlari/sympler
#
# Copyright 2002-2013, 
# David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
# and others authors stated in the AUTHORS file in the top-level 
# source directory.
#
# SYMPLER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SYMPLER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SYMPLER. If not, see <http://www.gnu.org/licenses/>.
#
# Please cite the research papers on SYMPLER in your own publications. 
# Check out the PUBLICATIONS file in the top-level source directory.
#
# You are very welcome to contribute extensions to the code. Please do 
# so by making a pull request on https://github.com/kauzlari/sympler
# 
#

import os
import glob
from scipy import *
import xml.dom.minidom


### General initialization

basedir = os.getcwd()
PROFILEDIR = "profiles"
VELOCITYDIR = "velocities"
#VTK2PROFILE = "/nfs/simulation/simstud/pastewka/src/sympler/src/tools/vtk2profile"
VTK2PROFILE = "/Volumes/Daten/HiWi/dpd/sympler/src/tools/vtk2profile"


### Parsing input file

def filename2glob(fn):
    pos = fn.rfind('.')
    return fn[0:pos] + "_*" + fn[pos:]


class SymplerInputParser:
    def __init__(self, filename):
        self.doc = xml.dom.minidom.parse(filename)
        self.vtks = self.doc.getElementsByTagName("OutputVTK")

        for i in self.vtks:
            self.vtk_filenames = filename2glob(i.getAttribute("fileName"))

        # Extract information
        el = self.doc.getElementsByTagName("Controller")
        self.dt = float(el[0].getAttribute("dt"))
        self.timesteps = float(el[0].getAttribute("timesteps"))

        el = self.doc.getElementsByTagName("ThermostatPeters")
        self.dissipation = float(el[0].getAttribute("dissipation"))

        el = self.doc.getElementsByTagName("ParticleCreatorRandom")
        self.n_particles = float(el[0].getAttribute("nParticles"))
        self.density = float(el[0].getAttribute("density"))

        el = self.doc.getElementsByTagName("OutputVTK")
        self.filenames = filename2glob(el[0].getAttribute("fileName"))

        el = self.doc.getElementsByTagName("GridAveragerStructured")
        self.avg_over = float(el[0].getAttribute("avgOver"))
        
        self.width = pow(self.n_particles/self.density, 1./3)


### Conversion utilities

def createProfiles(curpath, filenames):
    try:
        os.mkdir(curpath + "/" + PROFILEDIR)
    except:
        pass

    try:
        os.mkdir(curpath + "/" + PROFILEDIR + "/" + VELOCITYDIR)
    except:
        pass

    allprofiles = []
    for i in range(len(filenames)):
        curfile = curpath + "/" + PROFILEDIR + "/" + VELOCITYDIR + "_%i.out" % i
        os.system(VTK2PROFILE + " " + filenames[i] + " fluid:velocity:mean > %s" % curfile)
        allprofiles += [curfile]

    return allprofiles
    

### First harmonic calculation

def calculateFirstHarmonics(profilename, width):
    # Use only z-velocity
    profile = io.read_array(profilename)[:,2]
    dx = width/len(profile)
    xvals = arange(dx/2, width, dx)

    # Multiply profile with sin and integrate
    mult = sin(xvals*2*pi/width)*profile
    return integrate.simps(mult, xvals)


### Exponential fitting

def residuals(p, y, x):
    (A, t) = p
#    gplt.plot(x, y, x, A*exp(-t*x))
    return y - A*exp(-t*x)

def fitExponential(filenames, dt, maxt, width):
    first_harmonics = []
    for i in filenames[0:len(filenames)-1]:
        first_harmonics += [calculateFirstHarmonics(i, width)]
    x = arange(dt/2, maxt, dt)
#    gplt.plot(x, first_harmonics)
    (p, cov_x, infodict, ier, mesg) = optimize.leastsq(residuals, (1, 1), args=(first_harmonics, x), full_output=True)
    return (p[1], cov_x[1][1])
    

subdirs = glob.glob(basedir + "/d*")
subdirs.sort()
viscosities = []

for curpath in subdirs:
    parser = SymplerInputParser(curpath + "/sim.in")

    print "Evaluating simulation for dissipation = %f" % parser.dissipation

    # Convert VTKs to space separated files
    print "--- Converting VTKs to space separated files"
    filenames = glob.glob(curpath + "/" + parser.filenames)
    filenames.sort()

    try:
	os.stat(curpath + "/" + PROFILEDIR)
	profiles = []
	for i in range(len(filenames)):
	    profiles += [curpath + "/" + PROFILEDIR + "/" + VELOCITYDIR + "_%i.out" % i]
    except OSError:
	profiles = createProfiles(curpath, filenames)
    
#    profiles = glob.glob(curpath + "/" + PROFILEDIR + "/" + VELOCITYDIR + "_*.out")

    # Calculate the first harmonics and fit them to an exponential
    print "--- Calculating first harmonics and fitting decay to an exponential"
    e = fitExponential(profiles, parser.dt*parser.avg_over, parser.dt*parser.timesteps, parser.width)
    viscosity = parser.density*e[0]*((parser.width/(2*pi))**2)
    error = parser.density*sqrt(e[1])*((parser.width/(2*pi))**2)

    print "=== viscosity = %f, error = %f" % (viscosity, error)

    viscosities += [(parser.dissipation, viscosity, error)]

viscosities.sort()
io.write_array("viscosities2.dat", viscosities)

