/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2013, 
 * David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
 * and others authors stated in the AUTHORS file in the top-level 
 * source directory.
 *
 * SYMPLER is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SYMPLER is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SYMPLER.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Please cite the research papers on SYMPLER in your own publications. 
 * Check out the PUBLICATIONS file in the top-level source directory.
 *
 * You are very welcome to contribute extensions to the code. Please do 
 * so by making a pull request on https://github.com/kauzlari/sympler
 * 
 */



#include <math.h>

#include "first_harmonic.h"

#define IDX_TIME 0
#define IDX_INTENSITY 1

/* Register this Postprocessor with the factory. */
const Postprocessor_Register<FirstHarmonic> first_harmonic("FirstHarmonic");



//---- Constructors/Destructor ----

FirstHarmonic::FirstHarmonic(Node *parent, Simulation *simulation)
	: Postprocessor(parent, simulation)
{
    init();
}


FirstHarmonic::~FirstHarmonic()
{
}



//---- Methods ----

void FirstHarmonic::setup()
{
  Postprocessor::setup();
  m_idx_y = m_input_format->attrByName(m_species + ":velocity:mean").index;
}


void FirstHarmonic::describeInput(DataFormat *input_format)
{
    // fixme!!! check format
    m_input_format = input_format;

    m_idx_time = m_input_format->attrByName("time").index;
    m_idx_x = m_input_format->attrByName("cell_positions").index;
    // I try to do next in setup()
//     m_idx_y = m_input_format->attrByName(m_species + ":velocity:mean").index;

/*    m_idx_left = m_input_format->attrByName("left").index;
    m_idx_right = m_input_format->attrByName("right").index;*/
//     m_idx_grid_size = m_input_format->attrByName("grid_size").index; 
//     -> length of positions array now
}


void FirstHarmonic::push(data_sp data)
{
  if(!negativeFound)
    {  
      data_sp new_data = m_output.newData();
      
      double /*left, right,*/ len, step, pos, mul;
      double value = 0;
      vector<point_t> *vp;
      vector<point_t> *cp;
      
      new_data->doubleByIndex(IDX_TIME) = data->doubleByIndex(m_idx_time);
      cp = data->vectorPointByIndex(m_idx_x).value();      
      vp = data->vectorPointByIndex(m_idx_y).value();

      int grid_size = (*cp).size();
                  
//       left = data->doubleByIndex(m_idx_left);
//       right = data->doubleByIndex(m_idx_right);
//       grid_size = data->intByIndex(m_idx_grid_size);
      
      assert(grid_size > 1);
      step = (*cp)[1].x - (*cp)[0].x;
      // are the cells aligned in x-direction?
      if(step > g_geom_eps)
      {
        // make sure there is not more than one cell in y-direction
        if ((*cp)[1].y - (*cp)[0].y > g_geom_eps) 
          throw gError("FirstHarmonic::push", "multiple cells in more than one direction" 
            " found. This is not allowed with FirstHarmonic");
        len = (*cp)[grid_size - 1].x - (*cp)[0].x + step;
        pos = (*cp)[0].x - step/2;
        mul = 2*M_PI/len;
      }
      else // so it's the y-direction
      {
        step = (*cp)[1].y - (*cp)[0].y;
        assert(step > g_geom_eps);
        len = (*cp)[grid_size - 1].y - (*cp)[0].y + step;
        pos = (*cp)[0].y - step/2;
        mul = 2*M_PI/len;
      
      }
        
      for (int i = 0; i < grid_size; i++) {
        /* Currently, we integrate approximating the profile by step functions.
           The profile should be interpolated linearly, quadratically? */
        value += 0.5*(sin((pos + step/2) * mul) /*+ sin(pos * mul)*/) * (*vp)[i].z;
MSG_DEBUG("FirstHarmonic::push", "pos=" << pos << "\nstep=" << step << "\nmul=" << mul << "\nsin((pos+step)*mul) = sin((" << pos << "+" << step << ")*" << mul << ")=" << sin((pos+step)*mul) << "\nsin(pos*mul)=" << sin(pos*mul) << "\nvalue = " << value << "\nfac=" << 0.5*(sin((pos + step) * mul) + sin(pos * mul)) << "\nvel=" << (*vp)[i].z);	
        pos += step;
      }
      // first check is for abort if noise to large
      // second check avoids logarithms of negative numbers      
      if(value < oldValue && value > 0)
    	{
	      new_data->doubleByIndex(IDX_INTENSITY) = value;
        distribute(new_data);
        oldValue = value;
      }
      else negativeFound = true;
    }
}

void FirstHarmonic::init()
{
  negativeFound = false;
  oldValue = HUGE_VAL; 
  
  m_properties.setClassName("FirstHarmonic");

  m_properties.setDescription(
    "Determines the first coefficient of a Fourier series of the velocity profile. "
    "Can only be used with GridAveragerStructured."
  );
  
  m_output.addAttribute("time", DataFormat::DOUBLE, false);
  m_output.addAttribute("intensity", DataFormat::DOUBLE, false);
  
  STRINGPC
    (species, m_species,
     "Defines the species for the case that"
     " the Meter, this Postprocessor belongs to, considers only one type of particles. The"
     " Meter can measure either all or exactly one species.");

  m_species = "UNDEF";

}
