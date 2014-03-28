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



#include <fstream>

#include "wall_container.h"
#include "wall_triangle.h"

#include "misc.h"

#define FRAME_OFF 0

//---- Constructors/Destructor ----

WallContainer::WallContainer(Reflector *reflector): m_reflector(reflector)
{
  for (int i = 0; i < SPACE_DIMS; i++)
    m_periodicFront[i] = m_periodicBack[i] = false;
}


WallContainer::~WallContainer()
{
  FOR_EACH
    (list<Wall*>,
     m_walls,
     delete *__iFE;
    );
}



//---- Methods ----

void WallContainer::vertexChanged()
{
	for (list<Wall*>::iterator i = m_walls.begin(); i != m_walls.end(); i++) {
		(*i)->vertexChanged();
	}
}


void WallContainer::stretchBy(const point_t &factor)
{
	MSG_DEBUG("WallContainer::stretchBy", "factor = " << factor);
	MSG_DEBUG("WallContainer::stretchBy", "m_bounding_box.corner1 = " << 
    m_bounding_box.corner1);
	
	for (vector<point_t>::iterator i = m_vertices.begin(); i != m_vertices.end(); i++) {
		for (int j = 0; j < SPACE_DIMS; j++)
			(*i)[j] *= factor[j];
	}

    vertexChanged();
	
	for (int i = 0; i < SPACE_DIMS; i++) {
		m_bounding_box.corner1[i] *= factor[i];
		m_bounding_box.corner2[i] *= factor[i];
	}
}

void WallContainer::moveVertices(const point_t& dist)
{
	MSG_DEBUG("WallContainer::moveVertices", "dist = " << dist);
	
	for (vector<point_t>::iterator i = m_vertices.begin(); i != m_vertices.end(); i++) {
		for (int j = 0; j < SPACE_DIMS; j++)
			(*i)[j] += dist[j];
	}

    vertexChanged();
	
	// for adding of a frame, the following is NOT enough, but additionally
	// two vertices have to be added and updateBoundingBox() has to be called
	for (int i = 0; i < SPACE_DIMS; i++) {
        m_bounding_box.corner1[i] += dist[i];
		m_bounding_box.corner2[i] += dist[i];
	}
}

#define ODD(x)  (x & 1)


bool WallContainer::isInside(const point_t &pos)
{  
  line_t tester;
  int c_up, c_down, direction;

  if (m_walls.empty())
    return true;

  if (!m_bounding_box.isInside(pos))
	return false;
  tester.from = pos;
  tester.to = pos;

  // currently (02/04/05), it does not matter, for which periodicity we check
  bool oddIsIn = true;
	if(m_periodicFront.x || m_periodicBack.x) {
		if(m_periodicFront.y || m_periodicBack.y) {
			if(m_periodicFront.z || m_periodicBack.z) {
        // we are in the case of an interior obstacle combined with PBC in all directions
        oddIsIn = false;
        direction = 0;
      } else direction = 2;
		} else direction = 1;
	} else direction = 0;
			
	
	tester.to[direction] = m_bounding_box.corner2[direction] + 10;

	c_up = intersectionsGeneral(tester);
	
	tester.to[direction] = m_bounding_box.corner1[direction] - 10;
	c_down = intersectionsGeneral(tester);

	// fixme!!! this is a preliminary fix for the case that the ray travels in a plane,
	// which is bound to two other planes so that the ray intersects with these two planes;
	// the fix solves the problem for, e.g.,  BoundaryStepStoch but it could fail for more 
	// complicated geometries, e.g., if the ray suddenly goes through a periodic area that 
	// it did not touch before
	// fixme!!! is the increment of 2*c_wt_dist_eps appropriate?
	// the fix tilts the ray iteratively
/*
	while (ODD(c_up) != ODD(c_down))
	{
 	  MSG_DEBUG("WallContainer::isInside", "in loop: c_wt_dist_eps=" << c_wt_dist_eps << endl << ", pos = "  << pos << ", direction = " << direction << endl << "up=" << c_up << ", down=" << c_down << endl << "tester.to=" << tester.to[direction] << "tester.from=" << tester.from[direction]);
		for(size_t i = 1; i < SPACE_DIMS; ++i) 
			tester.to[(direction+i)%SPACE_DIMS] += 2*c_wt_dist_eps;
		
		tester.to[direction] = m_bounding_box.corner2[direction] + 10;
		c_up = intersectionsGeneral(tester);
		
		tester.to[direction] = m_bounding_box.corner1[direction] - 10;
		c_down = intersectionsGeneral(tester);
		//if(ODD(c_up) == ODD(c_down)) MSG_DEBUG("WallContainer::isInside", "in loop: done");

	}
*/
	
  /* Note: The assertion fails when the point sits on a surface,
     so we're just checking if one of the test got a positive
     result. fixme!!! */
  //    assert(ODD(c_up/*+add*/) == ODD(c_down/*+add*/));
	
	if ((ODD(c_up) && oddIsIn) || !(ODD(c_up) || oddIsIn))
		return true;

  return false;
}

#undef ODD

bool WallContainer::isInside(const cuboid_t &cuboid, const double& range)
{
  //      MSG_DEBUG("WallContainer::isInside(cuboid)", "CALLED");
    /* Quick check. */
 
    // this check shifts the corners a little up and down and tests for both
    // if one is inside, the whole cuboid is treated as inside
    // this assures correct treatement of a cuboid with a corner lying exactly 
    // on a wall  
  point_t a1, b1, c1, d1, e1, f1, g1, h1;
  point_t a2, b2, c2, d2, e2, f2, g2, h2;

    a1 = a2 =  b1 = b2 = c1 = c2 = g1 = g2 = cuboid.corner1;
    d1 = d2 = e1 = e2 = f1 = f2 = h1 = h2 = cuboid.corner2;

//     if(g1.x > 16.1 && g1.x < 17.9 &&
//        g1.y > 1.6 && g1.y < 3.4 &&
//        g1.z > 0.5 && g1.z < 1.625)
//       MSG_DEBUG("WallContainer::isInside(cuboid)", "FOUND corner1: " << g1);

//     if(h1.x > 16.1 && h1.x < 17.9 &&
//        h1.y > 1.6 && h1.y < 3.4 &&
//        h1.z > 0.5 && h1.z < 1.625)
//       MSG_DEBUG("WallContainer::isInside(cuboid)", "FOUND corner2: " << h1);

    double shift = 100*c_wt_dist_eps;

    a1.x = a2.x = cuboid.corner2.x;
    b1.y = b2.y = cuboid.corner2.y;
    c1.z = c2.z = cuboid.corner2.z;

    d1.x = d2.x = cuboid.corner1.x;
    e1.y = e2.y = cuboid.corner1.y;
    f1.z = f2.z = cuboid.corner1.z;
   
    for(size_t i = 0; i < SPACE_DIMS; ++i)
      {
	a1[i] -= shift;
	a2[i] += shift;
	b1[i] -= shift;
	b2[i] += shift;
	c1[i] -= shift;
	c2[i] += shift;
	d1[i] -= shift;
	d2[i] += shift;
	e1[i] -= shift;
	e2[i] += shift;
	f1[i] -= shift;
	f2[i] += shift;
	g1[i] -= shift;
	g2[i] += shift;
	h1[i] -= shift;
	h2[i] += shift;
      }

// if ((isInside(a1) || isInside(a2) ||
//         isInside(b1) || isInside(b2) ||
//         isInside(c1) || isInside(c2) ||
//         isInside(d1) || isInside(d2) ||
//         isInside(e1) || isInside(e2) ||
//         isInside(f1) || isInside(f2) ||
//         isInside(g1) || isInside(g2) ||
//      isInside(h1) || isInside(h2)) && g1.x > 16.1 && g1.x < 17.9 &&
//        g1.y > 1.6 && g1.y < 3.4 &&
//        g1.z > 0.5 && g1.z < 1.625 )
//   MSG_DEBUG("WallContainer::isInside(cuboid)", g1 << " (POINT-CHECK1) true");

// if ((isInside(a1) || isInside(a2) ||
//         isInside(b1) || isInside(b2) ||
//         isInside(c1) || isInside(c2) ||
//         isInside(d1) || isInside(d2) ||
//         isInside(e1) || isInside(e2) ||
//         isInside(f1) || isInside(f2) ||
//         isInside(g1) || isInside(g2) ||
//      isInside(h1) || isInside(h2)) && h1.x > 16.1 && h1.x < 17.9 &&
//        h1.y > 1.6 && h1.y < 3.4 &&
//        h1.z > 0.5 && h1.z < 1.625 )
//   MSG_DEBUG("WallContainer::isInside(cuboid)", h1 << " (POINT-CHECK2) true");


    
    if (isInside(a1) || isInside(a2) ||
        isInside(b1) || isInside(b2) ||
        isInside(c1) || isInside(c2) ||
        isInside(d1) || isInside(d2) ||
        isInside(e1) || isInside(e2) ||
        isInside(f1) || isInside(f2) ||
        isInside(g1) || isInside(g2) ||
        isInside(h1) || isInside(h2) )
      return true;


    a1 = b1 = c1 = g1 = cuboid.corner1;
    d1 = e1 = f1 = h1 = cuboid.corner2;

    a1.x = cuboid.corner2.x;
    b1.y = cuboid.corner2.y;
    c1.z = cuboid.corner2.z;

    d1.x = cuboid.corner1.x;
    e1.y = cuboid.corner1.y;
    f1.z = cuboid.corner1.z;

//     if(g1.x > 4 && g1.x < 5 && g1.y > -1 && g1.y < 1 && g1.z > 0 && g1.z < 1) MSG_DEBUG("WallContainer::isInside(cuboid)", g1 << ", NOW ");

    // this check is always done, but in principle only necessary with 
    //wall particles
    // fixme!!! this check is not 100% safe (it tests only if corners of 
    //cuboid in range)
    if(isInWallRange(range, g1) || isInWallRange(range, h1) || 
       isInWallRange(range, a1) || isInWallRange(range, b1) ||
       isInWallRange(range, c1) || isInWallRange(range, d1) || 
       isInWallRange(range, e1) || isInWallRange(range, f1) )
      {
// 	if(g1.x > 4 && g1.x < 5 && g1.y > -1 && g1.y < 1 && g1.z > 0 && g1.z < 1) MSG_DEBUG("WallContainer::isInside(cuboid)", g1 << ", TRUE ");
	return true;
      }
    /* Long and slow check. */

    //    bool is_inside = false;

    for (list<Wall*>::iterator i = m_walls.begin(); i != m_walls.end(); i++) {
      if ((*i)->intersects(cuboid))
	{
// 	  if(g1.x > 16.1 && g1.x < 17.9 &&
// 	     g1.y > 1.6 && g1.y < 3.4 &&
// 	    g1.z > 0.5 && g1.z < 1.625)
// 	MSG_DEBUG("WallContainer::isInside(cuboid)", g1 << "(INTERSECT1) true ");
	 
// 	  if(h1.x > 16.1 && h1.x < 17.9 &&
// 	     h1.y > 1.6 && h1.y < 3.4 &&
// 	    h1.z > 0.5 && h1.z < 1.625)
// 	MSG_DEBUG("WallContainer::isInside(cuboid)", h1 << "(INTERSECT2) true ");
	 
	  return true;
	}
    }
    
    return false;
}


bool WallContainer::isInWallRange(const double& range, const point_t& point)
{
	for (list<Wall*>::iterator i = m_walls.begin(); i != m_walls.end(); i++)
		if((*i)->isInRange(range, point)) {
			return true;
		}

	return false;
}


inline bool hasCloseP(const list<point_t> &p, const point_t &h)
{
  for (list<point_t>::const_iterator i = p.begin(); i != p.end(); i++) {
    point_t diff = (*i) - h;
    if (diff.absSquare() < g_geom_eps)
      return true;
  }

  return false;
}


int WallContainer::intersectionsGeneral(const line_t &l)
{
  int n = 0;
  list<point_t> pos;
  point_t hit_pos;

  FOR_EACH
    (list<Wall*>,
     m_walls,
     if ((*__iFE)->intersects(l, hit_pos)) {
       /*The next avoids to count the same hit_point twice, which might 
	 happen due to triangle overlap when using small error-epsilons*/
       if (!hasCloseP(pos, hit_pos)) {
         pos.push_back(hit_pos);
         n++;
       }
     }
     );
	
  return n;
}


int WallContainer::intersections(list<Wall*> walls, const line_t &l, int dir)
{
  int n = 0;
  list<point_t> pos;
  point_t hit_pos;

  FOR_EACH
    (list<Wall*>,
     walls,
     if ((*__iFE)->intersectsParallelToDir(l, dir, hit_pos)) {
       if (!hasCloseP(pos, hit_pos)) {
         pos.push_back(hit_pos);
         n++;
       }
     }
     );
	
  return n;
}


void WallContainer::updateBoundingBox()
{  
   	if (m_vertices.empty())
        	throw gError("WallContainer::updateBoundingBox", "Need at least two vertices.");
   	vertexChanged();
	for (int i = 0; i < SPACE_DIMS; i++) {
		m_bounding_box.corner1[i] = HUGE_VAL;
		m_bounding_box.corner2[i] = -HUGE_VAL;
	}

    	for (vector<point_t>::iterator v = m_vertices.begin(); v != m_vertices.end(); v++) 
	{
		for (int i = 0; i < SPACE_DIMS; i++) 
		{
			m_bounding_box.corner1[i] = min(m_bounding_box.corner1[i], (*v)[i]);
			m_bounding_box.corner2[i] = max(m_bounding_box.corner2[i], (*v)[i]);
		}
    	}
	m_cuboidVolume = 
		(m_bounding_box.corner2.x - m_bounding_box.corner1.x) * 
		(m_bounding_box.corner2.y - m_bounding_box.corner1.y) * 
		(m_bounding_box.corner2.z - m_bounding_box.corner1.z);
// 	MSG_DEBUG("WallContainer::updateBoundingBox", "end: corner1=" << m_bounding_box.corner1 << ", corner2=" << m_bounding_box.corner2);

 for (list<Wall*>::iterator i = m_walls.begin(); i != m_walls.end(); i++)
 {
   (*i)->setBoundaryData(m_periodicFront, m_bounding_box.corner2 - m_bounding_box.corner1);
 }

}


/*--- VTK export ---*/

void WallContainer::toVTK(string filename)
{
    ofstream s;

    s.open(filename.c_str());
    toVTK(s);
    s.close();
}


void WallContainer::toVTK(ostream &s)
{
    int n_ints;

    s << "# vtk DataFile Version 2.0" << endl;
    s << "SYMPLER simulation geometry information " << endl;
    s << "ASCII" << endl;
    s << "DATASET UNSTRUCTURED_GRID" << endl;

    s << "POINTS " << m_vertices.size() << " double" << endl;

    for (vector<point_t>::iterator i = m_vertices.begin(); i != m_vertices.end(); i++) {
        s << i->x << " " << i->y << " " << i->z << endl;
    }
    s << endl;

    n_ints = 0;
    FOR_EACH
        (list<Wall*>,
         m_walls,
         n_ints += (*__iFE)->vtkNVertices()+1;
         );

    s << "CELLS " << m_walls.size() << " " << n_ints << endl;
    FOR_EACH
        (list<Wall*>,
         m_walls,
         (*__iFE)->toVTK(s);
         );
    s << endl;

    s << "CELL_TYPES " << m_walls.size() << endl;
    FOR_EACH
        (list<Wall*>,
         m_walls,
         s << (*__iFE)->vtkCellType() << endl;
         );
    s << endl;

    s << "CELL_DATA " << m_walls.size() << endl;
    s << "VECTORS normal double" << endl;

    FOR_EACH
        (list<Wall*>,
         m_walls,
         point_t normal;
         
         normal = ((Wall*) (*__iFE))->normal();

         s << normal.x << " " << normal.y << " " << normal.z << endl;
         );
    s << endl;
}
