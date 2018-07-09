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



#include "boundary_stl.h"

#include "pc_inlet.h"
#include "stl_ascii.h"
#include "stl_binary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"

#include "threads.h"

/* Register this Boundary with the factory. */
const Boundary_Register<BoundarySTL> boundary_stl("BoundarySTL");


#define STL_ASCII "ascii"
#define STL_BINARY "binary"


//---- Constructors/Destructor ----

BoundarySTL::BoundarySTL(Phase *phase): BoundaryWithInlet(phase)
{
	init();
}


BoundarySTL::~BoundarySTL()
{
}

//---- Methods ----

void BoundarySTL::init()
{
  m_properties.setClassName("BoundaryWithInlet");
  m_properties.setName("BoundarySTL");

  m_properties.setDescription("Load a .STL stereolithography file.");

  for (int i = 0; i < 3; i++) {
    /* Register all properties */
    m_properties.addProperty
      ("scale" + string(1, 'X'+i), PropertyList::DOUBLE, 
       &m_scale[i],
       new PLCDoubleGreater(0),
       string(1, 'x'+i) + "-scaling factor for the geometry.");
  }


  STRINGPCINF
    (name, m_filename,
     "Name of the .STL geometry file to be loaded.");

  BOOLPC
    (invertNormals, m_invert_normals,
     "Invert the direction of the surface normal vectors.");

  m_properties.addProperty
    ("format", PropertyList::STRING,
     &m_format,
     new PLCStringList2(STL_ASCII, STL_BINARY),
     "Type of the STL: 'ascii' or 'binary'.");

  m_properties.addProperty
    ("inlet", PropertyList::STRING,
     &m_inlet_str,
     new PLCStringList4("none", "x", "y", "z"),
     "Direction of the flow through the inlet. Possible values: 'none' = no inlet, "
     "'x', 'y' and 'z'.");

  m_properties.addProperty
    ("inletPosition", PropertyList::STRING,
     &m_inlet_pos_str,
     new PLCStringList2("top", "bottom"),
     "Location of the inlet. Possible values are 'top' or 'bottom'.");

//   DOUBLEPC
//     (inletLength, m_inlet_length, 0,
//      "Length of the inlet.");

  m_filename = "default.stl";
  m_invert_normals = true;
  m_format = "ascii";

  m_inlet_str = "none";
  m_inlet_pos_str = "bottom";

//   m_inlet_length = 1;

  for (int i = 0; i < 3; i++)
    m_scale[i] = 1;
}


void BoundarySTL::setup()
{
  bool_point_t periodicFront = { { { false, false, false } } };
  bool_point_t periodicBack = { { { false, false, false } } };

  BoundaryWithInlet::setup();

//   if(m_inlet_str != "none")
//   {
//     int dir = m_inlet_str[0] - 'x';
//     if(m_inlet_pos_str == "bottom")
//     {
//       periodicBack[dir] = true;
//     }
//     else
//     {
//       assert(m_inlet_pos_str == "top");
//       periodicFront[dir] = true;
//     }
//   }

  m_container = new WallContainer(m_reflectors[0]);

  if (m_format == STL_ASCII) {
    STLAsciiParser loader(m_container);

    loader.setInvertNormals(m_invert_normals);

    m_container = loader.read(m_filename);
  } else {
    STLBinaryParser loader(m_container);

    loader.setInvertNormals(m_invert_normals);

    m_container = loader.read(m_filename);
  }

  m_container->setPeriodicityFront(periodicFront);
  m_container->setPeriodicityBack(periodicBack);

  m_container->stretchBy(m_scale);

  m_container->updateBoundingBox();

	MSG_DEBUG
    ("BoundarySTL::read",
     "box: " << m_container->boundingBox().corner1 << " - " << m_container->boundingBox().corner2);
}



//---- Cell subdivision ----

#define OUTLET 0
#define NEIGHBOR 1

void BoundarySTL::setup(Simulation* sim, ManagerCell *mgr)
{
    //    BoundaryWithInlet::setup(sim, mgr);
	
    
    /* Cell subdivision */
	
    region_t *inlet, *r;
    bool_point_t periodic = { { { false, false, false } } };

    /* Do we have to construct an inlet??? */
    if (m_inlet_str != "none") {
        int dir;
        double coord;
        ParticleCreatorInlet *pc = NULL;

        MSG_DEBUG("BoundarySTL::setup", "Adding inlet.");

        for(vector<ParticleCreator*>::iterator pcIter = m_pcList.begin();
            pcIter != m_pcList.end(); ++pcIter)
        {
          if(typeid(**pcIter) == typeid(ParticleCreatorInlet))
          {
            if(pc) 
              throw gError("BoundarySTL::setup", "Please define only one ParticleCreatorInlet");
            pc = (ParticleCreatorInlet*) (*pcIter);
          }
        }
        if(!pc) 
          throw gError("BoundarySTL::setup", "Since 'inlet' is not \"none\", you must define a ParticleCreatorInlet.");
        
//         if (typeid(*pcList[0]) != typeid(ParticleCreatorInlet))
//           throw gError
//             ("BoundarySTL::setup"
//              "Please do only use ParticleCreatorInlet with BoundarySTL if you wish to specify an inlet.");
// 
//         pc = (ParticleCreatorInlet*) pcList[0];
// 
// 
//         MSG_DEBUG("BoundarySTL::setup", "Using particle creator: " << pc->className());

        dir = m_inlet_str[0] - 'x';

        m_inlet_normal.assign(0);
        if (m_inlet_pos_str == "bottom") {
            m_inlet_normal[dir] = 1;
            coord = m_container->boundingBox().corner1[dir];
        } else {
            m_inlet_normal[dir] = -1;
            coord = m_container->boundingBox().corner2[dir];
        }

        
        //        MSG_DEBUG("BoundartSTL::setup", "m_inlet_normal = " << m_inlet_normal);

        /* First, delete all walls that might hinder flow through the inlet. */
        list<Wall*> &walls = m_container->walls();
        for (list<Wall*>::iterator i = walls.begin(); i != walls.end(); ) {
            if (((*i)->normal() - m_inlet_normal).abs() < g_geom_eps) {
                /* Fixme!!! Does only work right now. */
                WallTriangle *wt = (WallTriangle*) (*i);
                bool del = true;

                for (int j = 0; j < 3; j++) {
                    //                    MSG_DEBUG("BoundarySTL::setup", "... = " << wt->corner(j)[dir] - coord);
                    if (fabs(wt->corner(j)[dir] - coord) >= g_geom_eps)
                        del = false;
                }

                if (del) {
                    list<Wall*>::iterator d = i;
                    int vs[3];
                    i++;

                    //                    MSG_DEBUG("BoundarySTL::setup", "Deleting wall.");

                    for (int j = 0; j < 3; j++) {
                      vs[j] = m_inlet_base_surface.needVertex(wt->corner(j));
                    }

                    m_container->deleteWall(d);

                    m_inlet_base_surface.addWall
                      (new WallTriangle(&m_inlet_base_surface, NULL, vs[0], vs[1], vs[2]));
                } else
                    i++;
            } else
                i++;
        }

        /* Complicated version, add surfaces */

        MSG_DEBUG("BoundarySTL::setup", "Finding pairs.");
        vector< pair<int, int> > pairs;

        point_t middle_point = { { { 0, 0, 0 } } };
        size_t n_points = 0;
        for (list<Wall*>::iterator i = walls.begin(); i != walls.end(); i++) {
            WallTriangle *w = ((WallTriangle*) (*i));
            vector<int> corners;

            for (int j = 0; j < 3; j++) {
                if (fabs(w->corner(j)[dir] - coord) < g_geom_eps)
                    corners.push_back(w->cornerVertex(j));
            }

            //            MSG_DEBUG("BoundarySTL::setup", "corners.size() = " << corners.size());
            if (corners.size() == 2) {
                pair<int, int> p;

                p.first = corners[0];
                p.second = corners[1];

                middle_point += m_container->vertex(p.first);
                middle_point += m_container->vertex(p.second);

                n_points += 2;

                pairs.push_back(p);
            } else if (corners.size() == 3)
                throw gError("BoundarySTL::setup", "FATAL: corners.size() = 3 -> should not happpen.");
        }

        middle_point /= n_points;

        MSG_DEBUG("BoundarySTL::setup", "pairs.size() = " << pairs.size());
        MSG_DEBUG("BoundarySTL::setup", "Creating walls.");
        vector<int> partners;
        partners.resize(m_container->vertices().size()+pairs.size(), -1);
        for (vector< pair<int, int> >::iterator i = pairs.begin(); i != pairs.end(); i++) {
            int a, b;
            point_t n1, n2, mv;

            a = i->first;
            b = i->second;

            if (partners[a] == -1) {
                point_t p = m_container->vertex(a);
                p += -m_inlet_length*m_inlet_normal;
                partners[a] = m_container->addVertex(p);
            }
            if (partners[b] == -1) {
                point_t p = m_container->vertex(b);
                p += -m_inlet_length*m_inlet_normal;
                partners[b] = m_container->addVertex(p);
            }

            n1 = (m_container->vertex(b)-m_container->vertex(a)).cross
              (m_container->vertex(partners[a])-m_container->vertex(b));
            n2 = (m_container->vertex(partners[a])-m_container->vertex(b)).cross
              (m_container->vertex(partners[b])-m_container->vertex(partners[a]));
            mv = middle_point-m_container->vertex(b);

            /* Does the vector n1 point towards the middle? */
            if (n1 * mv > 0) {
              m_container->addWall
                (new WallTriangle(m_container, m_container->reflector(),
                                  a, b, partners[a]));
            } else {
              m_container->addWall
                (new WallTriangle(m_container, m_container->reflector(),
                                  b, a, partners[a]));
            }

            /* Does the vector n2 point towards the middle? */
            if (n2 * mv > 0) {
              m_container->addWall
                (new WallTriangle(m_container, m_container->reflector(),
                                  b, partners[a], partners[b]));
            } else {
              m_container->addWall
                (new WallTriangle(m_container, m_container->reflector(),
                                  b, partners[b], partners[a]));
            }
        }

        /* Now, readd the closing surface */
        MSG_DEBUG("BoundarySTL::setup", "Readding closing surface.");

        for (list<Wall*>::iterator i = m_inlet_base_surface.walls().begin(); i != m_inlet_base_surface.walls().end(); ++i) {
          WallTriangle *wall = (WallTriangle*) *i;
          int v1, v2, v3;

          v1 = m_container->findVertex(wall->corner(0) - m_inlet_length * m_inlet_normal, g_geom_eps);
          v2 = m_container->findVertex(wall->corner(1) - m_inlet_length * m_inlet_normal, g_geom_eps);
          v3 = m_container->findVertex(wall->corner(2) - m_inlet_length * m_inlet_normal, g_geom_eps);

          m_container->addWall
            (new WallTriangle(m_container, m_container->reflector(),
                              v1, v2, v3));
        }

        /* Easy version */
/*        vector<point_t> &vertices = m_container->vertices();
        for (vector<point_t>::iterator i = vertices.begin(); i != vertices.end(); i++) {
            if (fabs((*i)[dir] - coord) < g_geom_eps) {
                (*i)[dir] -= m_inlet_normal[dir]*m_inlet_length;
            }
            }*/

        m_container->vertexChanged();
        m_container->updateBoundingBox();
    
        m_proposedSize = m_container->boundingBox().corner2;
        
        // next will also call Boundary::setup => this may lead to addition of a frame
        BoundaryWithInlet::setup(sim, mgr);

        point_t inlet_c1, inlet_c2, r_c1, r_c2;

	Phase* phase = (Phase*) m_parent;

        if (m_inlet_pos_str == "bottom") {
            inlet_c1 = m_container->boundingBox().corner1;
            inlet_c2 = m_container->boundingBox().corner2;
            inlet_c2[dir] = inlet_c1[dir] + m_inlet_length;

            r_c1 = m_container->boundingBox().corner1;
//             r_c1[dir] += m_inlet_length;
            r_c2 = m_container->boundingBox().corner2;
            
            // take the frame into account if existent
            if(m_frontFrame[dir]) inlet_c2[dir] += phase->pairCreator()->interactionCutoff();
//             if(m_frontFrame[dir]) inlet_c2[dir] += sim->maxCutoff;

            r_c1[dir] = inlet_c2[dir];
                        
        } else {
            inlet_c1 = m_container->boundingBox().corner1;
            inlet_c2 = m_container->boundingBox().corner2;
            inlet_c1[dir] = inlet_c2[dir] - m_inlet_length;

            r_c1 = m_container->boundingBox().corner1;
            r_c2 = m_container->boundingBox().corner2;
//             r_c2[dir] -= m_inlet_length;
            
            // take the frame into account if existent
            if(m_endFrame[dir]) inlet_c1[dir] -= phase->pairCreator()->interactionCutoff();
//             if(m_endFrame[dir]) inlet_c1[dir] -= sim->maxCutoff;
            
            r_c2[dir] = inlet_c1[dir];
        }

                
        /* *New style* inlets with constant density */
//        periodic[dir] = true;
        periodic[dir] = false;
        inlet = mgr->cellSubdivide(sim->maxCutoff, inlet_c1, inlet_c2, periodic, 1, pc, 
          dir, (int) m_inlet_normal[dir], P_CREATE);

        pc->setRegion(inlet);

//        mgr->setPeriodicity(periodic);

        periodic[dir] = false;
        r = mgr->cellSubdivide(sim->maxCutoff, r_c1, r_c2, periodic, 0, pc, dir, 
          (int) -m_inlet_normal[dir], P_DELETE);
            
        MSG_DEBUG("BoundarySTL::setup", "r->n_cells = " << r->n_cells <<
                  ", inlet->n_cells = " << inlet->n_cells);

//        ManagerCell::connect
//            ((int) m_inlet_normal[dir]*dir-(m_inlet_normal[dir] < 0 ? 1 : 0), inlet, r, periodic, OUTLET);
        ManagerCell::connect
            ((int) m_inlet_normal[dir]*dir-(m_inlet_normal[dir] < 0 ? 1 : 0), inlet, r, periodic, NEIGHBOR);
    } else {
        m_proposedSize = m_container->boundingBox().corner2;
        
        // next could add a frame
        BoundaryWithInlet::setup(sim, mgr);

        r = mgr->cellSubdivide
          (sim->maxCutoff, m_container->boundingBox().corner1, m_container->boundingBox().corner2, periodic, 0);
    }

    mgr->cellSubdivisionFinished();
    
    mgr->assignContainer(m_container);
    
    int counter = 0;

    MSG_DEBUG("BoundarySTL::setup", "Writing geometry information to '" << m_geometry_filename << "'.");
    m_container->toVTK(m_geometry_filename);
}
