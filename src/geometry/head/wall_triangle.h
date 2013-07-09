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



#ifndef __WALL_TRIANGLE_H
#define __WALL_TRIANGLE_H

#include <gsl/gsl_poly.h>

#include "wall.h"
#include "integrator_position.h"

//#define TRACK_PARTICLE 183

class Cell;
// class IntegratorPosition;

extern double c_wt_dist_eps;
extern double c_wt_time_eps;

/*!
 * A wall piece in the shape of a triangle
 */
class WallTriangle: public Wall
{
 protected:
  /*!
   * Vertices of the corners.
   */
  int m_corners_v[3];

  /*!
   * Position of the corners
   */
  point_t m_corners[3];
	
  /*!
   * Vector connecting the corners.
   */
  point_t m_sides[3];
	
  /*! 
   * Normal vector to the triangle.
   */
  point_t m_surface_normal;

  /*!
   * Normal vectors to the triangles sides lying in the surface of the triangle.
   */
  point_t m_side_normals[3];
	
  /*!
   * Precomputed values for collision tests.
   * r means m_corners[0], "origin" of the triangle, n is m_normal
   */
  double m_ndotr;
	
	
  /*!
   * Checks whether a point in the plane of the triangle is really inside
   * the triangle.
   * @param s Check for this point
   */
  bool reallyInPlane(const point_t &s) const; 
/*  {
    // If the normal product of the side normal with the intersection point
    //   becomes smaller than zero the point is outside. 
    for (int i = 0; i < 3; i++)
      if (m_side_normals[i] * (s - m_corners[i]) < c_wt_dist_eps)
	{
	  //MSG_DEBUG("WallTriangle::reallyInPlane", "m_side_normals[i] * (s - m_corners[i]) = " 
	  //  << m_side_normals[i] * (s - m_corners[i]) << " with c_wt_dist_eps = " << c_wt_dist_eps);				
	  return false;
	}
    // MSG_DEBUG("WallTriangle::reallyInPlane", "it is really in plane");				
    return true;
  }*/
	
 public:
  /*!
   * Constructor
   * @param parent Pointer to the parent container this \a Wall belongs to
   * @param reflector Pointer to the \a Reflector to be used for this \a Wall
   * @param c1 First vertex
   * @param c2 Second vertex
   * @param c3 Third vertex
   */
  WallTriangle(WallContainer *parent, Reflector *reflector, int c1, int c2, int c3);

  /*!
   * Destructor
   */
  virtual ~WallTriangle();
	
  /*!
   * Bounding box for this triangle
   */
  virtual cuboid_t boundingBox();
	
  /*!
   * One of the vertices might have change, recalculate quantities
   * such as the face normal, etc.
   */
  virtual void vertexChanged() {
    initHelpers();
  }
	
  
  /*!
   * Intersection with a line
   * @param l Line to intersect with
   * @param hit_pos Resulting point of intersection
   */
  virtual bool intersects(const line_t &l, point_t &hit_pos) const;
  
  /*!
   * Intersection with a line
   * @param from Starting point of the line to intersect with
   * @param dir direction of the line
   * @param dist resulting perpendicular distance of 'from' to the wall 
   */
  virtual bool intersects(const point_t& from, const point_t& dir, double& dist) const;
      
  /*!
   * Intersection with a line which is parallel to one of the coordinate axises
   * @param c Line to intersect with
   * @param dir Coordinate axis
   * @param hit_pos Point of intersection
   */
  virtual bool intersectsParallelToDir(const line_t &l, int dir, point_t &hit_pos) const {
    int dir2, dir3;
    point_t mx, mn;

    dir2 = (dir+1)%3;
    dir3 = (dir+2)%3;

    mn[dir2] = -HUGE_VAL;
    mn[dir3] = -HUGE_VAL;
    mx[dir2] = HUGE_VAL;
    mx[dir3] = HUGE_VAL;

    for (int i = 0; i < 3; i++) {
      mn[dir2] = min(mn[dir2], m_corners[i][dir2]);
      mn[dir3] = min(mn[dir3], m_corners[i][dir3]);

      mx[dir2] = max(mx[dir2], m_corners[i][dir2]);
      mx[dir3] = max(mx[dir3], m_corners[i][dir3]);
    }

    if ((l.to[dir2] > mx[dir2]+g_geom_eps &&
	 l.to[dir3] > mx[dir3]+g_geom_eps &&
	 l.from[dir2] > mx[dir2]+g_geom_eps &&
	 l.from[dir3] > mx[dir3]+g_geom_eps) ||
	(l.to[dir2] < mn[dir2]-g_geom_eps &&
	 l.to[dir3] < mn[dir3]-g_geom_eps &&
	 l.from[dir2] < mn[dir2]-g_geom_eps &&
	 l.from[dir3] < mn[dir3]-g_geom_eps))
      return false;

    return intersects(l, hit_pos);
  }

  /*!
   * Intersection with a line
   * @param a Start of the line
   * @param b End of the line
   * @param hit_pos Point of intersection
   */
  virtual bool intersects(const point_t &a, const point_t &b, point_t &hit_pos) const {
    line_t l;
		
    l.from = a;
    l.to = b;
		
    return intersects(l, hit_pos);
  }

  /*!
   * Intersection with a line which is parallel to one of the coordinate axises
   * @param a Start of the line
   * @param b End of the line
   * @param dir Coordinate axis
   * @param hit_pos Point of intersection
   */
  virtual bool intersectsParallelToDir(const point_t &a, const point_t &b, int dir, point_t &hit_pos) const {
    line_t l;
		
    l.from = a;
    l.to = b;
		
    return intersectsParallelToDir(l, dir, hit_pos);
  }
	
  /*!
   * Intersection with a cuboid
   * @param c Cuboid to intersect with
   */
  virtual bool intersects(const cuboid_t &c) const;

//  /*!
//   * Calculates the integrator solveHitTimeEquation function
//   * @param p Solves for this particle
//   * @param force Force on this particle
//   * @param results Possible hit positions-vector
//   * @param integratorP Pointer to the IntegratorPosition used
  //   */ 
//  virtual void solveIntegratorEquation(const Particle *p, const point_t &force, vector<double>* results, IntegratorPosition* integratorP);
  
//  /*!
//   * Calculates the integrator solveHitTimeEquation function
//   * @param p Solves for this particle
//   * @param force Force on this particle
//   * @param integratorP Pointer to the IntegratorPosition used
//   * @param hit_pos The position where the particle hits the wall
  //   */
//  virtual void integratorHitPos(const Particle *p, const point_t &force, double dt, IntegratorPosition* integratorP, point_t &hit_pos);
  
  /*!
   * Check whether a particle will hit this wall during this time step
   * @param p Check for this particle
   * @param force Force on this particle
   * @param t_traveled The time the particle travelled until the hit occured
   * @param hit_pos The position where the particle hits the wall
   */   
  virtual bool hit(const Particle *p, const point_t &force, double &t_traveled, point_t &hit_pos, IntegratorPosition* integratorP){
    vector<double> results;
    
//     solveIntegratorEquation(p, force, &results, integratorP);

  integratorP->solveHitTimeEquation(this, p, force, &results);

  if (results.empty())
  {
     return false;
  }
  
  else
  {
    for (size_t i = 0; i < results.size(); ++i)
    {
      if(results[i] > p->dt)
        return false;
      
      if (results[i] > c_wt_time_eps)
      {
//         integratorHitPos(p, force, results[i], integratorP, hit_pos);
        
        integratorP->hitPos(/*this, */results[i], p, hit_pos, force);

        #ifdef TRACK_PARTICLE

	if (p->mySlot == TRACK_PARTICLE)
	MSG_DEBUG("WallTriangle::hit", "hit_pos = " << hit_pos);

	#endif

	#ifdef TRACK_PARTICLE

	bool inp = reallyInPlane(hit_pos);

	if (p->mySlot == TRACK_PARTICLE)
        {
	  MSG_DEBUG("WallTriangle::hit", "reallyInPlane = " << inp);

	  cout << toString() << endl;

	  for (i = 0; i < 3; i++)
	  {
	    MSG_DEBUG("WallTriangle::hit", "m_side_normals[i] * (hit_pos - m_corners[i]) = " << m_side_normals[i] * (hit_pos - m_corners[i]) << 
"with c_wt_dist_eps = " << c_wt_dist_eps);

	    MSG_DEBUG("WallTriangle::hit", "m_surface_normal * (hit_pos - m_corners[i]) = " << m_surface_normal * (hit_pos - m_corners[i]));
	  }
	}
	
	if (inp)
	{
	  t_traveled = results[i];
	  return inp;
	}

	#else

	if (reallyInPlane(hit_pos))
	{
// 	  bool inp = reallyInPlane(hit_pos);
	  t_traveled = results[i];
	  return true/*inp*/;
	}

	#endif

      }
    }
    return false;
  }
} 
    

    /*!
    * distance to plane, the wall is in
    */
    virtual double distToPlane(const point_t& p) const
    {
      return abs((p - m_corners[0])*m_surface_normal);
    }

    
  /*!
   * Is this wall close than \a range to position \a p?
   * WARNING: currently, this function is only used during initialisation; if it is 
   * meant to be used during simulation time, it has to be optimised concerning 
   * performance.
   * @param range Distance from this wall
   * @param p Position to check distance from
   */
    virtual bool isInRange(const double& range, const point_t& p) const
  {
    // distance to plane, the wall is in
    double dist = abs((p - m_corners[0])*m_surface_normal);
			
//     if(dist <=1 && p.x>285.458 && p.x < 285.460 && p.y>234.711 && p.y < 234.713 && p.z>273.468 && p.z < 273.470) MSG_DEBUG("WallTriangle::isInRange", "dist =" << dist << ", m_corners[0]=" << m_corners[0] << ", m_corners[1]=" << m_corners[1] << ", m_corners[2]=" << m_corners[2] << ", m_surface_normal=" << m_surface_normal << ", p=" << p);
			
    if(dist > range) return false;
    else

      // NEW STYLE
      {
	// the perpendicular projection of p to the plane
	point_t proj = p + dist*m_surface_normal;
	// see below for what the next one is needed
	bool tempBool = true;
	// loop over the side(normal)s
	for (int i = 0; i < 3; ++i)
	  {
	    // the distance of the projection to the current side
	    // <0 is outside the triangle, >0 is inside
	    // (because m_side_normals[i] always points inside the triangle)
	    double projDist = m_side_normals[i] * (proj - m_corners[i]);
	    // is projDist on the 'backside' of this side?
	    // If the normal product of the side normal with 'proj'
	    // becomes smaller than zero the point is outside.
	    // comparison to zero should be OK because, no matter if it's maybe a very small
	    // positive or negative number, the function will end up with 'true' anyhow
// 	    if(dist <=1 && p.x>285.458 && p.x < 285.460 && p.y>234.711 && p.y < 234.713 && p.z>273.468 && p.z < 273.470) MSG_DEBUG("WallTriangle::isInRange", "p=" << p << ", projDist(" << i << ")=" << projDist << ", side=" << m_sides[i]); 
	    if(projDist < 0)
	      {
		// now we need the normal vector of the 'region-boarder';
		// since the boarder is perp. to the side, the normal vector can be
		// obtained directly from the side itself
		// by construction of m_sides[i] it always points from 
		// m_corners[i] to m_corbers[(i+1)%SPACE_DIMS]
		point_t boarderNormal = m_sides[i] / m_sides[i].abs();
		// distance of 'proj' to the boarder at m_corners[i]
		double dist2BorderFrom = boarderNormal * (proj - m_corners[i]);
		// is it in the 'corner zone' ?
		// I don't see any danger to compare with zero here
		if(dist2BorderFrom < 0)
		  {
// 		    if(dist <=1 && p.x>285.458 && p.x < 285.460 && p.y>234.711 && p.y < 234.713 && p.z>273.468 && p.z < 273.470) MSG_DEBUG("WallTriangle::isInRange", "CORNER ZONE: p=" << p << ", side=" << m_sides[i]);
		    // So the final distance we need is the one from the point to the current corner
		    point_t finalDistVec = p - m_corners[i];
		    // point_t corner2Proj = proj - m_corners[i];
		    // we may not return 'false' directly because for angles>90deg. it may
		    // happen that 'proj' does not really lie in the corner zone but
		    // in the neighbouring centre zone 
		    if(/*corner2Proj*/finalDistVec.abs() > range) 
		      {
// if(dist <=1 && p.x>285.458 && p.x < 285.460 && p.y>234.711 && p.y < 234.713 && p.z>273.468 && p.z < 273.470) MSG_DEBUG("WallTriangle::isInRange", "FIRST CORNER ZONE FALSE: p=" << p << ", side=" << m_sides[i]);
			tempBool = false;
		      }
		    else 
		      {
// 			if(dist <=1 && p.x>285.458 && p.x < 285.460 && p.y>234.711 && p.y < 234.713 && p.z>273.468 && p.z < 273.470) MSG_DEBUG("WallTriangle::isInRange", "FIRST CORNER ZONE TRUE: p=" << p << ", side=" << m_sides[i]);
			// So we found at least one point of the triangle which is in range, that's enough
			return true;
		      }
		  }
		else
		  {
		    // distance of 'proj' to the boarder at m_corners[(i+1)%SPACE_DIMS]
		    double dist2Border2 = boarderNormal * (proj - m_corners[(i+1)%SPACE_DIMS]);
		    // is it in the second 'corner zone' ?
		    // I don't see any danger to compare with zero here
		    // !!! must be the opposite comparison (besides the ==) to the previous one,
		    // since we use the same 'boarderNormal' !!!
		    if(dist2Border2 > 0)
		      {
// 			if(dist <=1 && p.x>285.458 && p.x < 285.460 && p.y>234.711 && p.y < 234.713 && p.z>273.468 && p.z < 273.470) MSG_DEBUG("WallTriangle::isInRange", "SECOND CORNER ZONE: p=" << p << ", side=" << m_sides[i]);
			// So the final distance we need is the one from the point to the current corner
			point_t finalDistVec = p - m_corners[(i+1)%SPACE_DIMS];
			// point_t corner2Proj = proj - m_corners[(i+1)%SPACE_DIMS];
			if(/*corner2Proj*/finalDistVec.abs() > range) tempBool = false;
			else return true;
		      }
		    else // so it is in the 'centre zone'
		      {
			// 					MSG_DEBUG("WallTriangle::isInRange", "centre zone: p=" << p << ", projDist(" << i << ")=" << projDist << ", side=" << m_sides[i]); 
			// take care, projDist is negative !!!					
			// if((projDist+range) < 0) return false;
			if(sqrt(projDist*projDist+dist*dist) > range) tempBool = false;
			else return true;
		      }
		  }
	      } // end: if(projDist < 0)
	  } // end: loop over the side(normal)s
	// 				MSG_DEBUG("WallTriangle::isInRange", "p=" << p << ", returning tmpBool=" << tempBool);
	return tempBool;
      }		
  }


  /*!
   * Return the position of a corners
   * @param i Number of the corner
   */
  const point_t &corner(int i) const {
    return m_corners[i];
  }

  /*!
   * Return the vertex index for a corner
   * @param i Number of the corner
   */
  int cornerVertex(int i) const {
    return m_corners_v[i];
  }
	
  /*!
   * Return the normal vector
   */
  virtual const point_t& normal() const {
    return m_surface_normal;
  }

    /*!
   * Return the normal vector
   */
  virtual double nDotR() const {
    return m_ndotr;
  }
  
  /*!
   * Return a vector that lies in the plane of the surface
   */
  virtual const point_t& inPlane() const {
    return m_side_normals[0];
  }

  /*!
   * Return the vector on the sides of the triangle
   * @param i Number of the side
   */
  const point_t &side(int i) const {
    return m_sides[i];
  }

  /*!
   * Return a string identifier for this wall piece
   */
  virtual string toString() const;
	
  /*!
   * Distance of a point to the plane, the triangle is in
   * @param p Compute distance to this point
   */
  virtual double distanceTo(const point_t& p) const {
    // uses the formula (r_vec - r1_vec)*nUnit_vec,
    // where r_vec is the point, r1_vec is a point on the plane, 
    // nUnit_vec is the unit normal vector of the plane
    return (p - m_corners[0])*m_surface_normal;
  }

  /*!
   * Write this wall to VTK
   * @param s Stream to write to
   */
  virtual void toVTK(ostream &s);

  /*!
  * Check, whether the two walls are the same
  */
  virtual bool operator==(const Wall& wall) const;
  
  /*!
   * Calculate side normals, the face normal, ...
   */
  void initHelpers();
  
};

#endif
