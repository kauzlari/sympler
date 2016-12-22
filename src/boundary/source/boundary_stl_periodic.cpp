/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2016, 
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
#include "boundary_stl_periodic.h"

#include "stl_ascii.h"
#include "stl_binary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"


/* Register this Boundary with the factory. */
const Boundary_Register<BoundaryStlPeriodic> boundary_stl_periodic("BoundaryStlPeriodic");

//---- Constructors/Destructor ----

BoundaryStlPeriodic::BoundaryStlPeriodic(Phase *phase): BoundaryArbitrary(phase)
{
  init();
}


BoundaryStlPeriodic::~BoundaryStlPeriodic()
{
}
#define STL_ASCII "ascii"
#define STL_BINARY "binary"

//---- Methods ----

void BoundaryStlPeriodic::init()
{
  MSG_DEBUG("BoundaryStlPeriodic", "Initializing boundary.");  
  m_properties.setClassName("BoundaryStlPeriodic");
  m_properties.setName("BoundaryStlPeriodic");
  m_properties.setDescription("A Boundary with periodic stl geometry. ");
  
  for (int i = 0; i < 3; i++) {
    /* Register all properties */
    m_properties.addProperty
      ("scale" + string(1, 'X'+i), PropertyList::DOUBLE, 
       &m_scale[i],
       new PLCDoubleGreater(0),
       string(1, 'x'+i) + "-scaling factor for the geometry.");
    
    m_properties.addProperty
      ("periodic" + string(1, 'X'+i), PropertyList::BOOLEAN, &m_periodic[i], NULL,
       "Simulation box is periodic in " + string(1, 'x'+i) + "-direction.");
    
    //Defaults
    m_periodic[i] = false;
    m_scale[i] = 1;
  }
  
  STRINGPCINF
    (name, m_filename,
     "Name of the .STL geometry file to be loaded. Geometry must completely lie "
     "in positive XYZ directions. Geometry must have inflow and outflow defined "
     "by extremal points. Inflow and outflow must be at same height and have"
     " same shape.");

  BOOLPC
    (invertNormals, m_invert_normals,
     "Invert the direction of the surface normal vectors.");

  m_properties.addProperty
    ("format", PropertyList::STRING,
     &m_format,
     new PLCStringList2(STL_ASCII, STL_BINARY),
     "Type of the STL: 'ascii' or 'binary'.");
  
  
  
  //Defaults
  
  m_filename = "default.stl";
  m_invert_normals = true;
  m_format = "ascii";  
}

    /*!
     * Read and scale the boundary
    */
void BoundaryStlPeriodic::setup()
{
  
    BoundaryArbitrary::setup();

    m_container = new WallContainer(m_reflectors[0]); 
    
    if (m_format == STL_ASCII) {
    STLAsciiParser loader(m_container);

    loader.setInvertNormals(m_invert_normals);

    m_container = loader.read(m_filename);
    } 
    else {
    STLBinaryParser loader(m_container);

    loader.setInvertNormals(m_invert_normals);

    m_container = loader.read(m_filename);
    }
 

    m_container->setPeriodicityFront(m_periodic);
    m_container->setPeriodicityBack(m_periodic);  
    m_container->stretchBy(m_scale);
    m_container->updateBoundingBox();

  


/*ALGORITHM
 * 1.- Find inflow outflow surfaces. 
 *Loop over all walls and find minimum and maximum corners along the axis defined as periodic.
 *If a wall is on the inlet or outlet all corners will share the same distance in this direction. 
 * for example c1(x1,y1,z1) c2(x1,y2,z2) c3(x1,y3,z3)    : All corners same x1!
 * 2.-
 * If this condition is met the wall will be deleted and with it, its associated reflectors.  
 */
    for (int j=0 ; j<3; j++){
   
        if(m_periodic[j] == true){
      
            double min, max, temp;
            list<Wall*>::iterator d; 
            int countermax, countermin; 
               
            for (list<Wall*>::iterator i = m_container->walls().begin(); i != m_container->walls().end(); i++) { 
                for(int p =0; p<3;p++){  
                    temp= (*i)->returnCorner(p)[j];
                    
                    if(temp > max)
                    //if((temp-max)>=g_geom_eps)
                        max= temp;                        
                    
                    if(temp< min)
                    //if((temp-min) <= g_geom_eps)
                        min= temp;                       
                     
                }   
            }
     
            for (list<Wall*>::iterator i = m_container->walls().begin(); i != m_container->walls().end(); ) { 
                countermin =0; countermax=0; 
                
                for(int p =0; p<3;p++){ 
                    temp= (*i)->returnCorner(p)[j];
                                      
                    if(temp==min) 
                    //if((temp-min)<= g_geom_eps)  
                         ++countermin;                
                
                    if(temp==max)
                    //if((temp-max)<=g_geom_eps)
                        ++countermax;    
                }

                if (countermin == 3 || countermax ==3 ){  
                    d=i; 
                    i++;
                    m_container->deleteWall(d);
              
                }
                else i++;
            }
       
  
        }  
    } 
     MSG_DEBUG("BoundaryStlPeriodic::setup", "Setup completed.");
}





//---- Cell subdivision ----


void BoundaryStlPeriodic::setup(Simulation* sim, ManagerCell *mgr)
{
  BoundaryArbitrary::setup(sim, mgr);
   MSG_DEBUG
    ("BoundaryStlPeriodic::setup",
     "corner1 = " << m_container->boundingBox().corner1 << ", " <<
     "corner2 = " << m_container->boundingBox().corner2);

  MSG_DEBUG("BoundaryStlPeriodic::setup", "Initializing cells.");
  
  region_t *r;
  point_t c1, c2 ; 
  
  c1 = m_container->boundingBox().corner1;
  c2 = m_container->boundingBox().corner2; 
  
  r = mgr->cellSubdivide(sim->maxCutoff, c1, c2, m_periodic, 1, NULL); // 
    
  mgr->cellSubdivisionFinished();
  
  MSG_DEBUG("BoundaryStlPeriodic::setup", "Cell subdivision Finished.");
 
  mgr->assignContainer(m_container);
  
  m_container->toVTK(m_geometry_filename);
}