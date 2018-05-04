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



#ifdef HAVE_SDL

#include <string>

#include "meter_live_coloured.h"

#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"

#define XRES 640
#define YRES 480

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterLiveColoured> meter_live_coloured("MeterLiveColoured");


struct Color {
    int R, G, B;
};

// white, red, green
const Color c_color_map[3] = {
    {255, 0, 0}, {255, 255, 255}, {255, 0, 0}
};


Color HSV2RGB(double h, double s, double v)
{
  Color c;
  int hi;
  double f, p, q, t;

  hi = (int) (h*6);

  if (hi >= 6)
    hi = 5;

  f = h*6 - (double) hi;
  p = v*(1-s);
  q = v*(1-f*s);
  t = v*(1-(1-f)*s);

  switch (hi) {
  case 0:
    c.R = (int) (255*v); c.G = (int) (255*t); c.B = (int) (255*p);
    break;
  case 1:
    c.R = (int) (255*q); c.G = (int) (255*v); c.B = (int) (255*p);
    break;
  case 2:
    c.R = (int) (255*p); c.G = (int) (255*v); c.B = (int) (255*t);
    break;
  case 3:
    c.R = (int) (255*p); c.G = (int) (255*q); c.B = (int) (255*v);
    break;
  case 4:
    c.R = (int) (255*t); c.G = (int) (255*p); c.B = (int) (255*v);
    break;
  case 5:
    c.R = (int) (255*v); c.G = (int) (255*p); c.B = (int) (255*q);
    break;
  default:
    abort();
  }

  return c;
}


//---- Constructors/Destructor ----

MeterLiveColoured::MeterLiveColoured(Simulation *simulation)
: Meter(simulation), m_display(NULL), m_x_axis(0), m_y_axis(1), m_x_res(XRES), m_y_res(YRES)
{
    init();
}


MeterLiveColoured::~MeterLiveColoured()
{
}



//---- Methods ----

void MeterLiveColoured::init()
{
  m_properties.setClassName("MeterLiveColoured");

  m_properties.setDescription("Graphical view of the simulation at runtime. The colour of the particles represents a scalar particle property chosen by 'colourIs'.");

  STRINGPC
    (species, m_species,
     "Species to display.");
    
  INTPC
    (xresolution, m_x_res, 0,
     "x-Resolution of the window.");
  INTPC
    (yresolution, m_y_res, 0,
     "y-Resolution of the window.");
  m_properties.addProperty
    ("xaxis", PropertyList::INT, &m_x_axis,
     new PLCIntList3(0, 1, 2),
     "Defines the spatial direction (0, 1 or 2) to be mapped to the x-axis.");
  m_properties.addProperty
    ("yaxis", PropertyList::INT, &m_y_axis,
     new PLCIntList3(0, 1, 2),
     "Defines the spatial direction (0, 1 or 2) to be mapped to the y-axis.");

  // NULL-pointer means, no constraint on the value
  m_properties.addProperty
    ("phi", PropertyList::DOUBLE, &m_phi, NULL,
     "Phi Euler angle for rotation of view.");
  m_properties.addProperty
    ("theta", PropertyList::DOUBLE, &m_theta, NULL,
     "Theta Euler angle for rotation of view.");
  m_properties.addProperty
    ("psi", PropertyList::DOUBLE, &m_psi, NULL,
     "Psi Euler angle for rotation of view.");

  STRINGPC
    (colourIs, m_colour_is,
     "The particles colour displays the scalar property given here.");

  m_properties.addProperty
    ("colourRangeMin", PropertyList::DOUBLE, &m_colour_range_min, NULL,
     "Lower bound of the colour map.");
	
  m_properties.addProperty
    ("colourRangeMax", PropertyList::DOUBLE, &m_colour_range_max, NULL,
     "Upper bound of the colour map.");

  m_species = "fluid";

  m_phi = m_theta = m_psi = 0;

  m_colour_is = "none";
  m_colour_range_min = 0;
  m_colour_range_max = 1;
}


struct z_sorting : public binary_function<double, double, bool> {
  size_t m_z_coord;
  z_sorting(size_t z_coord): m_z_coord(z_coord) { }
  bool operator()(const disp_point_t &a, const disp_point_t &b) {
    return a.r[m_z_coord] < b.r[m_z_coord];
  }
};


void MeterLiveColoured::measureNow(const double &time)
{
  Phase *phase = ((Simulation*) m_parent)->phase();

  if (!m_display) {
    double factor_x, factor_y;
    int x_res, y_res;
    point_t size;

    size = phase->boundary()->boundingBox().size();

    factor_x = (double) m_x_res/size[m_x_axis];
    factor_y = (double) m_y_res/size[m_y_axis];

    m_factor = min(factor_x, factor_y);

    x_res = (int) (m_factor*size[m_x_axis])+1;
    y_res = (int) (m_factor*size[m_y_axis])+1;

    MSG_DEBUG("MeterLiveColoured::measureNow", "m_x_res=" << m_x_res << ", x_res=" << x_res << ", m_y_res=" << m_y_res << ", y_res=" << y_res);
    m_x_res = x_res;
    m_y_res = y_res;
				
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
      throw gError("MeterLiveColoured::measureNow",
		   string("Could not initialize video: ") + SDL_GetError());

    MSG_DEBUG("MeterLiveColoured::measureNow", "Resolution: " << x_res << "x" << y_res);

    m_display = SDL_SetVideoMode(x_res+1, y_res+11, 0, SDL_SWSURFACE);
    if (!m_display)
      throw gError("MeterLiveColoured::measureNow",
		   string("Could not create window: ") + SDL_GetError());        
  }

  SDL_Rect r;
  point_t offset;

  offset = phase->boundary()->boundingBox().corner1;

  r.x = 0;
  r.y = 0;
  r.w = m_x_res + 1;
  r.h = m_y_res + 11;

  SDL_FillRect(m_display, &r, 0/*0xffffff*/);

  for (int x = 0; x < r.w; x++) {
    Color col = HSV2RGB((double) (4.*(r.w-x-1))/(6.*r.w), 1, 1);

    for (int y = m_y_res+2; y < m_y_res+11; y++) {
      DrawPixel
	(m_display,
	 x,
	 y,
	 col.R,
	 col.G,
	 col.B);
    }

    DrawPixel
      (m_display,
       x,
       m_y_res+1,
       255,
       255,
       255);
  }

  m_points.clear();
  FOR_EACH_PARTICLE_C
    (phase, m_colour, 
     disp_point_t p;

     p.r = __iSLFE->r;
     p.p = __iSLFE->tag.doubleByOffset(m_colour_is_offset);

     m_points.push_back(p);
     );

  size_t other = 0;

  while (other == m_x_axis || other == m_y_axis)
    other++;

  sort(m_points.begin(), m_points.end(), z_sorting(other));

  FOR_EACH
    (vector<disp_point_t>,
     m_points,
     point_t &p = __iFE->r;
     double &prop = __iFE->p;

     double h = 4.*(m_colour_range_max - prop)/(6.*(m_colour_range_max-m_colour_range_min));

     if (h < 0)
       h = 0;
     else if (h > 1)
       h = 1;

     Color col = HSV2RGB(h, 1, 1);
				
     int x_coord = (int) (m_factor*(p[m_x_axis]-offset[m_x_axis]));
     int y_coord = (int) (m_factor*(p[m_y_axis]-offset[m_y_axis]));
					
     if (x_coord >= m_x_res ||
	 y_coord >= m_y_res) {
     } else {
       DrawPixel(m_display,
		 x_coord,
		 y_coord,
		 col.R,
		 col.G,
		 col.B);
       DrawPixel(m_display,
		 1+x_coord,
		 y_coord,
		 col.R,
		 col.G,
		 col.B);
       DrawPixel(m_display,
		 x_coord,
		 1+y_coord,
		 col.R,
		 col.G,
		 col.B);
       DrawPixel(m_display,
		 1+x_coord,
		 1+y_coord,
		 col.R,
		 col.G,
		 col.B);
     }
     );


  SDL_Flip(m_display);
}


void MeterLiveColoured::setup()
{
  if(((Simulation*) m_parent)->phase()->nColours() > 3) 
    throw gError("MeterLiveColoured cannot handle more than 3 colours");

  Meter::setup();
	
  m_rotation(0, 0) = cos(m_psi)*cos(m_phi) - cos(m_theta)*sin(m_phi)*sin(m_psi);
  m_rotation(0, 1) = cos(m_psi)*sin(m_phi) + cos(m_theta)*cos(m_phi)*sin(m_psi);
  m_rotation(0, 2) = sin(m_psi)*sin(m_theta);
  m_rotation(1, 0) = - sin(m_psi)*cos(m_phi) - cos(m_theta)*sin(m_phi)*cos(m_psi);
  m_rotation(1, 1) = - sin(m_psi)*sin(m_phi) + cos(m_theta)*cos(m_phi)*cos(m_psi);
  m_rotation(1, 2) = cos(m_psi)*sin(m_theta);
  m_rotation(2, 0) = sin(m_theta)*sin(m_phi);
  m_rotation(2, 1) = - sin(m_theta)*cos(m_phi);
  m_rotation(2, 2) = cos(m_theta);

  m_colour = M_MANAGER->getColour(m_species);

  m_colour_is_offset = Particle::s_tag_format[m_colour].attrByName(m_colour_is).offset;
}

//---- static Functions ----

void MeterLiveColoured::DrawPixel(SDL_Surface *screen, int x, int y, Uint8 R, Uint8 G, Uint8 B)
{
  //  for(size_t i = 0; i < 20000000; ++i);  

  // MSG_DEBUG("MeterLiveColoured::DrawPixel", "begin");

  Uint32 color = SDL_MapRGB(screen->format, R, G, B);

  if (SDL_MUSTLOCK(screen))
    if (SDL_LockSurface(screen) < 0)
      return;

  switch (screen->format->BytesPerPixel) {
  case 1: {
    Uint8 *bufp;

    bufp = (Uint8*) screen->pixels + y*screen->pitch +x;
    *bufp = color;
  }
    break;
            
  case 2: {
    Uint16 *bufp;

    bufp = (Uint16*) screen->pixels + y*screen->pitch/2 +x;
    *bufp = color;
  }
    break;

  case 4: {
    Uint32 *bufp;

    bufp = (Uint32*) screen->pixels + y*screen->pitch/4 +x;
    *bufp = color;
  }
    break;
  }

  if (SDL_MUSTLOCK(screen))
    SDL_UnlockSurface(screen);
				
  // MSG_DEBUG("MeterLiveColoured::DrawPixel", "end");
}

#endif
