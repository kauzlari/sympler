#! /usr/bin/env python
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

def create_input_file_from(template, dict, filename):
    f = open(template)
    lines = f.readlines()
    f.close()
    
    f = open(filename, "w")

    for i in lines:
        f.write(i % dict)

    f.close()


# We want to run simulation with variations of
#    dissipation  (Peters thermostat)
#    kappa        (heat conduction)
# and determine macroscopic parameters
#    viscosity
#    heat conduction


dict = {}
#SYMPLER = "~/sympler-050315"
SYMPLER = "~/sympler-050401"
basepath = os.getcwd()

# paras contains tuples (dissipation, kappa)

paras = []
for dissipation in [0.3, 0.4, 0.5, 1.0, 2.0, 4.0, 8.0, 12, 15, 20, 25, 30, 35, 40, 45, 50, 100, 500, 1000, 5000]:
    paras += [dissipation]

densities = [2, 6, 10, 14]
n_runs = 5

for density in densities:
    try:
        os.mkdir("rho%i" % density)
    except:
        print("--- Path already existed: %s" % path)
    os.chdir("rho%i" % density)
    
    for dissipation in paras:
        #for i in range(n_runs):
            dict["density"] = density
            dict["dissipation"] = dissipation

            path = "d%f" % dissipation

            try:
                os.mkdir(path)
                os.mkdir("%s/results" % path)
                os.mkdir("%s/results/grid" % path)        
            except:
                print("--- Path already existed: %s" % path)
            
            os.chdir(path)

            create_input_file_from("%s/sim.in.template" % basepath, dict, "sim.in")

            print("=== Running sympler for dissipation = %f" % dissipation)
            
            try:
                about = os.stat("OUT")
                print("There is an OUT file in this directory. sympler is not being run.")
            except OSError:
                os.system("%s sim.in > OUT" % SYMPLER)

            os.chdir("..")

    os.chdir("..")
        

