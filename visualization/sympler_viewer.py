#!/usr/bin/env python
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

import sys

if len(sys.argv) != 2:
    print "Syntax: sympler_viewer.py <input-filename>"
    sys.exit(1)

print "Loading os, gtk, gtk.glade, gobject and glob."

import os
import gtk
import gtk.glade
import gobject

import glob


print "Loading vtk. This may take a while."

import vtk
from vtk.util.colors import *
from GtkGLExtVTKRenderWindow import *


print "Loading support modules."

from GladeApp import *
from SymplerInputParser import *



class SymplerRenderer:
    def __init__(self):
        self.files = []
        self.timestep = -1

        self.geometry = vtk.vtkUnstructuredGridReader()
        self.geometry_mapper = vtk.vtkDataSetMapper()
        self.geometry_actor = vtk.vtkActor()

        self.positions = vtk.vtkUnstructuredGridReader()
        self.positions_mapper = vtk.vtkDataSetMapper()
        self.positions_actor = vtk.vtkActor()

        self.geometry_mapper.SetInput(self.geometry.GetOutput())
        self.positions_mapper.SetInput(self.positions.GetOutput())
        self.positions_mapper.ScalarVisibilityOff()

        self.geometry_actor.SetMapper(self.geometry_mapper)
        self.geometry_actor.GetProperty().SetOpacity(0.2)

        self.positions_actor.SetMapper(self.positions_mapper)
        self.positions_actor.GetProperty().SetColor(light_grey)

        self.ren = vtk.vtkRenderer()
        self.ren_win = GtkGLExtVTKRenderWindow()
        self.ren_win.GetRenderWindow().AddRenderer(self.ren)

        self.ren.AddActor(self.geometry_actor)
        self.ren.AddActor(self.positions_actor)
        self.ren.SetBackground(0.1, 0.2, 0.4)
#        self.ren_win.SetSize(700, 700)

    def setGeometry(self, fn):
        self.geometry.SetFileName(fn)


    def setInteractorStyle(self, style):
        self.iren.SetInteractorStyle(style)


    def setTimesteps(self, pos_files):
        self.files = pos_files
        print self.files
        self.selectTimestep(0)


    def size(self):
        return len(self.files)


    def selectTimestep(self, i):
        self.timestep = i
        self.positions.SetFileName(self.files[i])
        self.ren_win.Render()


    def render(self):
        self.ren_win.Render()


    def writeJPEG(self, fn):
        filter = vtk.vtkWindowToImageFilter()
        filter.SetInput(self.ren_win)
        writer = vtk.vtkJPEGWriter()
        writer.SetInput(filter.GetOutput())
        writer.SetFileName(fn)
        writer.Write()


#    def eventLoop(self):
#        self.iren.Initialize()
#        self.ren_win.Render()
#        self.iren.Start()



class SymplerViewer(GladeApp):
    def __init__(self, input_fn):
        GladeApp.__init__(self, "sympler_viewer.glade", "view_window")
        self.input_parser = SymplerInputParser(input_fn)
        
        self.renderer = SymplerRenderer()
        self.renderer.setGeometry("geometry.vtk")
        self.renderer.setTimesteps(self.input_parser.positionFilenames())
        self.renderer.selectTimestep(0)

        self.scale.set_range(0, self.renderer.size()-1)
    
        self.vtk_view = self.renderer.ren_win
        self.vtk_view.show()
        self.vbox_view.pack_start(self.vtk_view)
        self.vbox_view.reorder_child(self.vtk_view, 0)


    def on_scale_value_changed(self, range, *args):
        if round(range.get_value()) != self.renderer.timestep:
            self.renderer.selectTimestep(int(round(range.get_value())))


    def on_play_button_clicked(self, button, *args):
        gobject.timeout_add(1000, self.animation_timer)


    def on_stop_button_clicked(self, button, *args):
        pass


    def animation_timer(self):
        if self.renderer.timestep == self.renderer.size()-1:
            return gtk.FALSE
        else:
            self.scale.set_value(self.renderer.timestep+1)
            return gtk.TRUE
        


if __name__ == "__main__":
    mv = SymplerViewer(sys.argv[1])
    mv.run()
