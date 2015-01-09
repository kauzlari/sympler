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



#include <sstream>

#include "wall_container.h"
#include "wall_triangle.h"
#include "integrator_position.h"
#include "integrator_velocity_verlet.h"
#include "integrator_velocity_verlet_disp.h"
#include "integrator_static.h"

double c_wt_dist_eps = -1e-5;
double c_wt_time_eps = 0;

//---- Constructors/Destructor ----

WallTriangle::WallTriangle(WallContainer *parent, Reflector *reflector, int c1, int c2, int c3)
    : Wall(parent, reflector)
{
    m_vtk_cell_type = 5;
    m_vtk_n_vertices = 3;

    m_corners_v[0] = c1;
    m_corners_v[1] = c2;
    m_corners_v[2] = c3;
	
    initHelpers();
}

WallTriangle::~WallTriangle()
{
}



//---- Methods ----

void WallTriangle::initHelpers()
{
    for (int i = 0; i < 3; i++) {
        m_corners[i] = ((WallContainer*) m_parent)->vertex(m_corners_v[i]);
    }

    for (int i = 0; i < SPACE_DIMS; i++) {
        m_sides[i] = m_corners[(i+1)%SPACE_DIMS] - m_corners[i];
    }
	
    m_surface_normal = m_sides[0].cross(m_sides[1]);
    m_surface_normal /= m_surface_normal.abs();

    //    MSG_DEBUG("WallTriangle::initHelpers", "m_surface_normal = " << m_surface_normal);

    m_ndotr = m_corners[0]*m_surface_normal;
	
    /* m_side_normals point INWARDS into the triangle. */
    for (int i = 0; i < SPACE_DIMS; i++) {
        m_side_normals[i] = m_surface_normal.cross(m_sides[i]);
	m_side_normals[i] /= m_side_normals[i].abs();
    }
		
/*		point_t test0 = {{{0,0,0}}};
		point_t test1 = {{{1,1,1}}};		
MSG_DEBUG("WallTriangle::initHelpers", "corner1 = " << m_corners[0] << endl 
	<< "corner2 = " << m_corners[1]
	<< endl << "corner3 = " << m_corners[2] << endl << "surface normal = " << m_surface_normal 
	<< endl << "distance to (0, 0, 0) = " << distanceTo(test0) << endl 
	<< "distance to (1, 1, 1) = " << distanceTo(test1) << endl);
if(isInRange(0.5, test0)) cout << "WallTriangle::initHelpers: (0, 0, 0) is in range" << endl;
if(isInRange(0.5, test1)) cout << "WallTriangle::initHelpers: (1, 1, 1) is in range" << endl;	*/

}


cuboid_t WallTriangle::boundingBox()
{
	cuboid_t c;
	
	for (int i = 0; i < SPACE_DIMS; i++) {
		c.corner1[i] = min(m_corners[0][i], min(m_corners[1][i], m_corners[2][i]));
		c.corner2[i] = max(m_corners[0][i], max(m_corners[1][i], m_corners[2][i]));
	}
	
	return c;
}


/*void WallTriangle::stretchBy(const point_t &factor)
{
	for (int i = 0; i < SPACE_DIMS; i++) {
		for (int j = 0; j < 3; j++) {
			m_corners[j][i] *= factor[i];
		}
	}
	
	initHelpers();
}*/


string WallTriangle::toString() const
{
	stringstream s;
	
	s << "WallTriangle" << endl;
    s << "Corners (vertices): " << m_corners_v[0] << ", " << m_corners_v[1] << ", " << m_corners_v[2] << endl;
    s << "Corners (positions): " << m_corners[0] << ", " << m_corners[1] << ", " << m_corners[2] << endl;
	s << "Sides: " << m_sides[0] << ", " << m_sides[1] << ", " << m_sides[2] << endl;
	s << "Side normals: " << m_side_normals[0] << ", " << m_side_normals[1] << ", " << m_side_normals[2] << endl;
    s << "Surface normal: " << m_surface_normal << endl;
	
	return s.str();
}


// FIXME: inline?
bool WallTriangle::intersects(const line_t &l, point_t &hit_pos) const
{
  point_t d = l.to - l.from;
  double np = m_surface_normal*d;
	
  if (!np)
    return false;
	
  double a = (m_ndotr - m_surface_normal*l.from)/np;
	
	// MSG_DEBUG("WallTriangle::intersects", "a = " << a);
	
  if (a < 0 || a > 1)
    return false;
	
  /* Vector to the intersection point relativ to the origin of the plane. */
  hit_pos = l.from + d*a;
	// MSG_DEBUG("WallTriangle::intersects", "hit_pos = " << hit_pos);	
//	cout << s << endl;
	
  return reallyInPlane(hit_pos);
}

// FIXME: inline?
bool WallTriangle::intersects(const point_t& from, const point_t& dir, /*point_t &hit_pos,*/ double& dist) const
{
  double np = m_surface_normal*dir;
	
  if (!np)
    return false;
	
  // absolute non-normalised (compare to intersects() above) normal distance of 'from' to the plane
  // notice that m_surface_normal points INSIDE, but the normal vector we need here points in fact OUTSIDE !
  dist = (-m_surface_normal*from + m_ndotr);
	
//   MSG_DEBUG("WallTriangle::intersects", "dist = " << dist << ", np = " << np << " for: from = " << from << ", m_ndotr = " << m_ndotr << ", m_surface_normal = " << m_surface_normal << ", dir = " << dir);
	
  double a = dist/np;
  
  if (a < 0 || a > 1)
    return false;
	
  // if a is correct, we can now compute the absolute value of dist
  dist = fabs(dist);
  
  point_t hit_pos = from + dir*a;
//    MSG_DEBUG("WallTriangle::intersects", "hit_pos = " << hit_pos << ", m_periodicity = " << m_periodicity);	
 
  // check for periodicity
  
//   const bool_point_t& periodicity = ((WallContainer*) m_parent)->periodicityFront();
//   cuboid_t* tempBox = &(((WallContainer*) m_parent)->boundingBox());
//   boxSize = tempBox->corner1 - tempBox->corner2;

  
  for(size_t i = 0; i < SPACE_DIMS; ++i)
  {
    if(m_periodicity[i])
    {
//        MSG_DEBUG("WallTriangle::intersects", "periodicity-check for dir = " << i << ", hit_pos[i] = " << hit_pos[i] << ", m_boxSize[i] = " << m_boxSize[i] << ", hit_pos[i]-box_size[i] = " << hit_pos[i]-m_boxSize[i]);
       
       // we add the constant since we hope that reallyInPlane(..) still reports an intersection since it also uses the constant in a conservative way, i.e., if the particle is slightly outside, it doesn't matter
       if(hit_pos[i] > m_boxSize[i] - c_wt_dist_eps /*the constant is < 0*/)
      {
        hit_pos[i] -= m_boxSize[i];  
//          MSG_DEBUG("WallTriangle::intersects", "periodicity-check: > TRUE");
      }
      else if(hit_pos[i] < c_wt_dist_eps /*the constant is < 0*/)
      {
        hit_pos[i] += m_boxSize[i];  
//          MSG_DEBUG("WallTriangle::intersects", "periodicity-check: < TRUE");
      }
    }    
  }
//    MSG_DEBUG("WallTriangle::intersects", "hit_pos after periodicity-check = " << hit_pos);	
 
  /* Vector to the intersection point relativ to the origin of the plane. */
//   return reallyInPlane(from + dir*a);
//	cout << s << endl;
	
  return reallyInPlane(hit_pos);
}

// FIXME: inline again?
bool WallTriangle::reallyInPlane(const point_t &s) const {
  /* If the normal product of the side normal with the intersection point
  minus a side-point becomes smaller than zero the point is outside. */
  for (int i = 0; i < 3; i++)
    if (m_side_normals[i] * (s - m_corners[i]) < c_wt_dist_eps)
  {
//      MSG_DEBUG("WallTriangle::reallyInPlane", "NOT IN PLANE: i = " << i << ", m_side_normals[i] * (s - m_corners[i]) = " << m_side_normals[i] * (s - m_corners[i]) << " with c_wt_dist_eps = " << c_wt_dist_eps << "\nside normal = " << m_side_normals[i] << "\ncorner = " << m_corners[i] << "\ns = " << s << "\ns - m_corners[i] = " << s - m_corners[i]);				
    return false;
  }
  // MSG_DEBUG("WallTriangle::reallyInPlane", "it is really in plane");				
  return true;
}


#define TI(a, b)  if (intersects(a, b, dummy)) {/* if(g1.x > 16.1 && g1.x < 17.9 && g1.y > 1.6 && g1.y < 3.4 && g1.z > 0.5 && g1.z < 1.625) MSG_DEBUG("WallTriangle::intersects(cuboid)", "TI: " << g1);*/ return true;}


inline bool __isInside(const cuboid_t &c, const point_t &pos) {
    for (int i = 0; (i < SPACE_DIMS); i++) {
            if ((pos[i] < c.corner1[i]) || (pos[i] > c.corner2[i]))
                return false;
    }
    return true;
}


bool WallTriangle::intersects(const cuboid_t &c) const
{
/* Quick check: does at least one of the triangle corners lie in the cuboid? */
//   point_t g1 = c.corner1;
	if (__isInside(c, m_corners[0]) || __isInside(c, m_corners[1]) || __isInside(c, m_corners[2]))
	  {

// 	    if(g1.x > 16.1 && g1.x < 17.9 &&
// 	     g1.y > 1.6 && g1.y < 3.4 &&
// 	    g1.z > 0.5 && g1.z < 1.625)
// 	  MSG_DEBUG("WallTriangle::intersects(cuboid)", "INSIDE: " << g1);

		return true;
	  }


#if 0


    /* Fixme!!! Quick and dirty and not exact. */
    for (int i = 0; i < SPACE_DIMS; i++) {
        if ((m_corners[0][i] > c.corner2[i]+g_geom_eps &&
             m_corners[1][i] > c.corner2[i]+g_geom_eps &&
             m_corners[2][i] > c.corner2[i]+g_geom_eps) ||
            (m_corners[0][i] < c.corner1[i]-g_geom_eps &&
             m_corners[1][i] < c.corner1[i]-g_geom_eps &&
             m_corners[2][i] < c.corner1[i]-g_geom_eps))
            return false;
    }

// if
//   (g1.x > 16.1 && g1.x < 17.9 &&
//    g1.y > 1.6 && g1.y < 3.4 &&
// 	    g1.z > -1.4 && g1.z < 0.6)
// 	      MSG_DEBUG("WallTriangle::intersects(cuboid)", "END = true" << endl << "corner0 = " << m_corners[0] << endl << "corner1 = " << m_corners[1] << endl << "corner2 = " << m_corners[2]);

//return true;
#endif

//#if 0
    /* Check if one of the edges of the cuboid cuts through
       the triangle. */
	point_t pb, pc, pd, pe, pf, pg, dummy;

    /*   pf         c2
          ----------
         /|        /|              y
        / |       / |
       /  |      /  |              |
    pc ---------- pe|              *- x
      |   /pd    |  |             /
      |  /       |  /pg          z  
      | /        | /           
      |/         |/
       ----------
      c1        pb
     */
	
	pb = pc = pd = pe = pf = pg = c.corner1;
	
	pb.x = c.corner2.x;
	pc.y = c.corner2.y;
	pd.z = c.corner2.z;
	
	pe.x = c.corner2.x;
	pe.y = c.corner2.y;
	pf.y = c.corner2.y;
	pf.z = c.corner2.z;
	pg.x = c.corner2.x;
	pg.z = c.corner2.z;
	
	TI(c.corner1, pb);
	TI(c.corner1, pc);
	TI(c.corner1, pd);
	
	TI(c.corner2, pe);
	TI(c.corner2, pf);
	TI(c.corner2, pg);
	
	TI(pg, pd);
	TI(pg, pb);
	TI(pf, pd);
	TI(pf, pc);
	TI(pe, pb);
	TI(pe, pc);


    /* Check if one of the edges of the triangle cuts
       through the cuboid. */
    if (c.intersects(m_corners[0], m_corners[1]) ||
        c.intersects(m_corners[1], m_corners[2]) ||
        c.intersects(m_corners[2], m_corners[0]))
      {
/*if
  (g1.x > 16.1 && g1.x < 17.9 &&
   g1.y > 1.6 && g1.y < 3.4 &&
	    g1.z > 0.5 && g1.z < 1.625)
  {
	      MSG_DEBUG("WallTriangle::intersects(cuboid)", "END = true" << endl << "corner0 = " << m_corners[0] << endl << "corner1 = " << m_corners[1] << endl << "corner2 = " << m_corners[2]);
	      if(c.intersects(m_corners[0], m_corners[1])) cout << "0_1 true" << endl; 
	      if(c.intersects(m_corners[1], m_corners[2])) cout << "1_2 true" << endl; 
	      if(c.intersects(m_corners[2], m_corners[1])) cout << "2_0 true" << endl; 
}    */
	return true;
      }

	return false;
//#endif
}


void WallTriangle::toVTK(ostream &s)
{
    s << "3 " << m_corners_v[0] << " " << m_corners_v[1] << " " << m_corners_v[2] << endl;
}


// void WallTriangle::solveIntegratorEquation(const Particle *p, const point_t &force, vector<double>* results, IntegratorPosition* integratorP)
// {
//   integratorP->solveHitTimeEquation(this, p, force, results);
// }
// 
// void WallTriangle::integratorHitPos(const Particle *p, const point_t &force, double dt, IntegratorPosition* integratorP, point_t &hit_pos)
// {
//   integratorP->hitPos(this, dt, p, hit_pos, force);
// }

bool WallTriangle::operator==(const Wall& wall) const
{
  bool same = true;
  if(typeid(*this) != typeid(wall))
    return false;
  else
  {
    for(size_t i = 0; i < SPACE_DIMS; ++i)
    {
      if(m_corners_v[i] != ((WallTriangle&) wall).m_corners_v[i])
      {
        same = false;
        i = 2;
      } 
    }
  }
  return same; 
}

