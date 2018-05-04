/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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



#ifndef __METER_LIVE_COLOURED_H
#define __METER_LIVE_COLOURED_H

#ifdef HAVE_SDL

#include <SDL/SDL.h>

#include "meter.h"

using namespace std;



/* template<typename T> */
/* class math_matrix_t { */
/*  protected: */
/*   int m_n, m_m; */
/*   T *m_data; */
	
/*  public: */
/*   math_matrix_t(int n = SPACE_DIMS, int m = SPACE_DIMS) : m_n(n), m_m(m) { */
/*     m_data = new T[n * m]; */
/*   } */
/*     ~math_matrix_t() { */
/*       delete m_data; */
/*     } */

/*     T &operator()(int y, int x) { */
/*       return m_data[y * m_m + x]; */
/*     } */
/*     const T &operator()(int y, int x) const { */
/*       return m_data[y * m_m + x]; */
/*     } */

/*     point_t operator*(const point_t &p) { */
/*       point_t r; */
	
/*       for (int j = 0; j < m_n; j++) { */
/* 	r[j] = 0; */
		
/* 	for (int i = 0; i < m_m; i++) { */
/* 	  r[j] += operator()(j, i)*p[i]; */
/* 	} */
/*       } */
	
/*       return r; */
/*     } */

/*     math_matrix_t<T> operator*(const math_matrix_t<T> &a) { */
/*       math_matrix_t<T> r(m_n, a.m_m); */
	
/*       for (int j = 0; j < m_n; j++) { */
/* 	for (int i = 0; i < a.m_m; i++) { */
/* 	  for (int k = 0; k < m_m; k++) { */
/* 	    r(j, i) += operator()(j, k)*a(k, i); */
/* 	  } */
/* 	} */
/*       } */
	
/*       return r; */
/*     } */
/* }; */


/* struct disp_point_t { */
/*   point_t r; */
/*   double p; */
/* }; */


/* typedef math_matrix_t<double> matrix_t; */


//---- Classes ----

class Simulation;

class MeterLiveColoured : public Meter
{
protected:
  /*!
   * The SDL display
   */
  SDL_Surface *m_display;

  /*!
   * The x-axis on screen is this axis in reality
   */
  size_t m_x_axis;

  /*!
   * The y-axis on screen is this axis in reality
   */
  size_t m_y_axis;

  /*!
   * x-resolution
   */
  int m_x_res;

  /*!
   * y-resolution
   */
  int m_y_res;
  
  /*!
   * Phi Euler angle for rotation of the scene
   */
  double m_phi;

  /*!
   * Theta Euler angle for rotation of the scene
   */
  double m_theta;

  /*!
   * Psi Euler angle for rotation of the scene
   */
  double m_psi;

  /*!
   * Scaling factor for coordinate to screen transformation
   */
  double m_factor;

  /*!
   * The specie to display
   */
  string m_species;

  /*!
   * The particle color to display
   */
  size_t m_colour;

  /*!
   * Color corresponds to this degree of freedom
   */
  string m_colour_is;
  
  /*!
   * Tag offset for the color information
   */
  size_t m_colour_is_offset;

  /*!
   * Smallest displayable value
   */
  double m_colour_range_min;

  /*!
   * Largest displayable value
   */
  double m_colour_range_max;
	
  /*!
   * Rotation matrix
   */
  matrix_t m_rotation;

  /*!
   * Points on screen
   */
  vector<disp_point_t> m_points;
    
  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Draw a pixel on the screen
   * @param screen SDL surface to draw on
   * @param x Pixel x position
   * @param y Pixel y position
   * @param R Red value
   * @param G Green value
   * @param B Blue value
   */
  static void DrawPixel(SDL_Surface *screen, int x, int y, Uint8 R, Uint8 G, Uint8 B);
    
public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   */
  MeterLiveColoured(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~MeterLiveColoured();

  /*!
   * Initialize rotation matrix
   */
  virtual void setup();

  /*!
   * Draw picture on screen
   */
  virtual void measureNow(const double& time);
};

#endif

#endif
